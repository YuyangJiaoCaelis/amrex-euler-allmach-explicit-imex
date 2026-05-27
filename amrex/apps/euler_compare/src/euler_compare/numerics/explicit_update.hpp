// Explicit finite-volume update with MUSCL reconstruction and Riemann fluxes.
int count_hllc_degenerate_star_states(const amrex::MultiFab& state, const RunConfig& cfg)
{
  if (cfg.riemann != RiemannKind::Hllc && cfg.riemann != RiemannKind::LowMachHllc &&
      cfg.riemann != RiemannKind::XieAmHllcP) {
    return 0;
  }

  amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<int> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const amrex::Real gamma = cfg.gamma;
  const int spatial_order = cfg.spatial_order;
  const SlopeLimiterKind slope_limiter = cfg.slope_limiter;
  const RiemannKind riemann = cfg.riemann;
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const src = state.const_array(mfi);
    reduce_op.eval(box, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
      const auto star_state_degenerate =
          [=] AMREX_GPU_DEVICE(const PrimitiveState& left, const PrimitiveState& right, int axis) noexcept {
        if (riemann == RiemannKind::LowMachHllc) {
          PrimitiveState corrected_left;
          PrimitiveState corrected_right;
          low_mach_correct_velocity_states(left, right, gamma, corrected_left, corrected_right);
          return hllc_star_state_is_degenerate(corrected_left, corrected_right, axis, gamma);
        }
        return hllc_star_state_is_degenerate(left, right, axis, gamma);
      };
      int count = 0;
      count += star_state_degenerate(reconstruct_state(src, i - 1, j, k, 0, +1, spatial_order, gamma, slope_limiter),
                                     reconstruct_state(src, i, j, k, 0, -1, spatial_order, gamma, slope_limiter), 0)
                   ? 1
                   : 0;
      count += star_state_degenerate(reconstruct_state(src, i, j, k, 0, +1, spatial_order, gamma, slope_limiter),
                                     reconstruct_state(src, i + 1, j, k, 0, -1, spatial_order, gamma, slope_limiter), 0)
                   ? 1
                   : 0;
      count += star_state_degenerate(reconstruct_state(src, i, j - 1, k, 1, +1, spatial_order, gamma, slope_limiter),
                                     reconstruct_state(src, i, j, k, 1, -1, spatial_order, gamma, slope_limiter), 1)
                   ? 1
                   : 0;
      count += star_state_degenerate(reconstruct_state(src, i, j, k, 1, +1, spatial_order, gamma, slope_limiter),
                                     reconstruct_state(src, i, j + 1, k, 1, -1, spatial_order, gamma, slope_limiter), 1)
                   ? 1
                   : 0;
      return {count};
    });
  }

  ReduceTuple values = reduce_data.value();
  int count = amrex::get<0>(values);
  amrex::ParallelDescriptor::ReduceIntSum(count);
  return count;
}

int euler_fv_one_step(amrex::MultiFab& state, const amrex::Geometry& geom,
                      const RunConfig& cfg, amrex::Real dt,
                      amrex::Real time = amrex::Real(0.0))
{
  BL_PROFILE("euler_fv_one_step()");

  amrex::MultiFab next(state.boxArray(), state.DistributionMap(), NCons, 0);
  fill_problem_ghosts(state, geom, cfg, time);
  const int hllc_degenerate_count = count_hllc_degenerate_star_states(state, cfg);

  const auto dx = geom.CellSizeArray();
  const amrex::Real inv_dx = amrex::Real(1.0) / dx[0];
  const amrex::Real inv_dy = amrex::Real(1.0) / dx[1];
  const amrex::Real gamma = cfg.gamma;
  const RiemannKind riemann = cfg.riemann;
  const int spatial_order = cfg.spatial_order;
  const SlopeLimiterKind slope_limiter = cfg.slope_limiter;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const src = state.const_array(mfi);
    auto const dst = next.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const ConservedState fx_l =
          reconstructed_interface_flux(src, i, j, k, 0, spatial_order, gamma, slope_limiter, riemann);
      const ConservedState fx_r =
          reconstructed_interface_flux(src, i + 1, j, k, 0, spatial_order, gamma, slope_limiter, riemann);
      const ConservedState fy_b =
          reconstructed_interface_flux(src, i, j, k, 1, spatial_order, gamma, slope_limiter, riemann);
      const ConservedState fy_t =
          reconstructed_interface_flux(src, i, j + 1, k, 1, spatial_order, gamma, slope_limiter, riemann);

      const ConservedState old = load_cons(src, i, j, k);
      ConservedState updated;
      updated.rho = old.rho - dt * ((fx_r.rho - fx_l.rho) * inv_dx + (fy_t.rho - fy_b.rho) * inv_dy);
      updated.mx = old.mx - dt * ((fx_r.mx - fx_l.mx) * inv_dx + (fy_t.mx - fy_b.mx) * inv_dy);
      updated.my = old.my - dt * ((fx_r.my - fx_l.my) * inv_dx + (fy_t.my - fy_b.my) * inv_dy);
      updated.E = old.E - dt * ((fx_r.E - fx_l.E) * inv_dx + (fy_t.E - fy_b.E) * inv_dy);
      store_cons(dst, i, j, k, updated);
    });
  }

  amrex::MultiFab::Copy(state, next, 0, 0, NCons, 0);
  return hllc_degenerate_count;
}

void populate_trial_fields(ImexStepDiagnostics& diagnostics,
                               const AdmissibilityDiagnostics& trial)
{
  diagnostics.trial_density_min = trial.rho_min;
  diagnostics.trial_pressure_min = trial.pressure_min;
  diagnostics.trial_internal_energy_min = trial.internal_energy_min;
  diagnostics.trial_nonfinite_count = trial.nonfinite_count;
  diagnostics.trial_positivity_count =
      trial.density_failure_count + trial.pressure_failure_count + trial.internal_energy_failure_count;
  diagnostics.rho_min = trial.rho_min;
  diagnostics.pressure_min = trial.pressure_min;
  diagnostics.internal_energy_min = trial.internal_energy_min;
  diagnostics.positivity_failure_count = diagnostics.trial_positivity_count;
  diagnostics.nonfinite_count = trial.nonfinite_count;
}

void configure_pressure_energy_check_diagnostics(ImexStepDiagnostics& diagnostics,
                                                 const RunConfig&,
                                                 const std::string& status,
                                                 const std::string& pressure_source)
{
  diagnostics.pressure_energy_check_mode = "report_solver_consistency";
  diagnostics.pressure_energy_check_update_policy =
      "diagnostic_output_only_state_update_unchanged";
  diagnostics.pressure_energy_check_status = status;
  diagnostics.pressure_energy_check_solved_pressure_source = pressure_source;
  diagnostics.pressure_energy_check_kinetic_proxy =
      "source_map_lagged_corrected_momentum";
}
