// Source-map pressure system, Picard step, and final closeout.
struct SourceMapPressureSystem {
  HostSparseMatrix matrix;
  Eigen::VectorXd rhs;
  int nonfinite_count = 0;
};

SourceMapPressureSystem source_map_build_pressure_system(const CoupledHostGrid& grid,
                                                         const std::vector<amrex::Real>& p,
                                                         const std::vector<SourceMapMomentum>& lagged_momentum,
                                                         const RunConfig& cfg,
                                                         const amrex::Geometry& geom,
                                                         amrex::Real dt)
{
  const int cell_count = grid.nx * grid.ny;
  SourceMapPressureSystem system;
  system.rhs = Eigen::VectorXd::Zero(cell_count);
  std::vector<Eigen::Triplet<amrex::Real>> triplets;
  triplets.reserve(static_cast<std::size_t>(cell_count * 11));
  for (int j = 0; j < grid.ny; ++j) {
    for (int i = 0; i < grid.nx; ++i) {
      const int row_index = grid.row(i, j);
      const SourceMapPressureRowData row =
          source_map_pressure_row_data(grid, p, lagged_momentum, i, j, cfg, geom, dt);
      system.nonfinite_count += row.nonfinite_count;
      system.rhs(row_index) = row.rhs;
      for (const SourceMapRowTerm& term : row.lhs_terms) {
        if (term.coeff != amrex::Real(0.0)) {
          triplets.emplace_back(row_index, grid.row(term.col_i, term.col_j), term.coeff);
        }
      }
    }
  }
  system.matrix = HostSparseMatrix(cell_count, cell_count);
  system.matrix.setFromTriplets(triplets.begin(), triplets.end());
  system.matrix.makeCompressed();
  return system;
}

struct SourceMapFinalEval {
  std::vector<amrex::Real> pressure;
  std::vector<ConservedState> state;
  amrex::Real residual_l1 = amrex::Real(0.0);
  amrex::Real residual_linf = amrex::Real(0.0);
  amrex::Real residual_relative_linf = amrex::Real(0.0);
  amrex::Real pressure_mismatch_l1 = amrex::Real(0.0);
  amrex::Real pressure_mismatch_linf = amrex::Real(0.0);
  amrex::Real pressure_mismatch_relative_linf = amrex::Real(0.0);
  amrex::Real density_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real solved_pressure_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real solved_pressure_max = -std::numeric_limits<amrex::Real>::infinity();
  amrex::Real eos_pressure_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real eos_pressure_max = -std::numeric_limits<amrex::Real>::infinity();
  amrex::Real internal_energy_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real internal_energy_max = -std::numeric_limits<amrex::Real>::infinity();
  amrex::Real momentum_correction_linf = amrex::Real(0.0);
  amrex::Real energy_correction_linf = amrex::Real(0.0);
  amrex::Real energy_closeout_delta_l1 = amrex::Real(0.0);
  amrex::Real energy_closeout_delta_linf = amrex::Real(0.0);
  amrex::Real energy_flux_linf = amrex::Real(0.0);
  amrex::Real energy_flux_divergence_linf = amrex::Real(0.0);
  amrex::Real eq78_divergence_linf = amrex::Real(0.0);
  amrex::Real eq78_divergence_scaled_linf = amrex::Real(0.0);
  int residual_count = 0;
  int nonfinite_count = 0;
  int positivity_count = 0;
  int density_failure_count = 0;
  int pressure_failure_count = 0;
  int internal_energy_failure_count = 0;
  int face_enthalpy_guard_count = 0;
  int face_enthalpy_nonfinite_count = 0;
  int energy_flux_nonfinite_count = 0;
};

SourceMapFinalEval evaluate_source_map_final(const CoupledHostGrid& grid,
                                             const std::vector<amrex::Real>& p,
                                             const RunConfig& cfg,
                                             const amrex::Geometry& geom,
                                             amrex::Real dt)
{
  SourceMapFinalEval eval;
  const int cell_count = grid.nx * grid.ny;
  eval.pressure = p;
  eval.state.assign(static_cast<std::size_t>(cell_count), ConservedState{});
  const auto dx = geom.CellSizeArray();
  const amrex::Real inv_dx = amrex::Real(1.0) / dx[0];
  const amrex::Real inv_dy = amrex::Real(1.0) / dx[1];
  const amrex::Real gamma_minus_one = cfg.gamma - amrex::Real(1.0);
  const std::vector<SourceMapMomentum> corrected_momentum =
      source_map_eq54_corrected_momentum_field(grid, p, dt, inv_dx, inv_dy, &cfg);
  amrex::Real residual_sum = amrex::Real(0.0);
  amrex::Real mismatch_sum = amrex::Real(0.0);
  amrex::Real energy_delta_sum = amrex::Real(0.0);

  for (int j = 0; j < grid.ny; ++j) {
    for (int i = 0; i < grid.nx; ++i) {
      CoupledFaceStats face_stats;
      const amrex::Real flux_e =
          source_map_face_energy_flux(grid, p, corrected_momentum, i, j, i + 1, j, 0,
                                      cfg, dt, inv_dx, inv_dy, &face_stats);
      const amrex::Real flux_w =
          source_map_face_energy_flux(grid, p, corrected_momentum, i - 1, j, i, j, 0,
                                      cfg, dt, inv_dx, inv_dy, &face_stats);
      const amrex::Real flux_n =
          source_map_face_energy_flux(grid, p, corrected_momentum, i, j, i, j + 1, 1,
                                      cfg, dt, inv_dx, inv_dy, &face_stats);
      const amrex::Real flux_s =
          source_map_face_energy_flux(grid, p, corrected_momentum, i, j - 1, i, j, 1,
                                      cfg, dt, inv_dx, inv_dy, &face_stats);
      eval.face_enthalpy_guard_count += face_stats.enthalpy_guard_count;
      eval.face_enthalpy_nonfinite_count += face_stats.enthalpy_nonfinite_count;
      eval.energy_flux_nonfinite_count += face_stats.flux_nonfinite_count;
      eval.energy_flux_linf = std::max(eval.energy_flux_linf, face_stats.flux_linf);

      const int row = grid.row(i, j);
      const ConservedState& star = grid.at(i, j);
      const SourceMapMomentum momentum =
          corrected_momentum[static_cast<std::size_t>(grid.row(i, j))];
      const amrex::Real div_flux = (flux_e - flux_w) * inv_dx + (flux_n - flux_s) * inv_dy;
      const amrex::Real energy = star.E - dt * div_flux;
      const amrex::Real kinetic =
          amrex::Real(0.5) * (momentum[0] * momentum[0] + momentum[1] * momentum[1]) / star.rho;
      const amrex::Real internal = energy - kinetic;
      const amrex::Real eos_pressure = gamma_minus_one * internal;
      const amrex::Real residual = p[static_cast<std::size_t>(row)] / gamma_minus_one - energy + kinetic;
      const amrex::Real mismatch = p[static_cast<std::size_t>(row)] - eos_pressure;
      const amrex::Real residual_scale =
          std::max({std::abs(p[static_cast<std::size_t>(row)]) / gamma_minus_one,
                    std::abs(energy), std::abs(kinetic), amrex::Real(1.0)});
      const amrex::Real mismatch_scale =
          std::max({std::abs(p[static_cast<std::size_t>(row)]), std::abs(eos_pressure),
                    std::abs(kinetic), amrex::Real(1.0)});

      ConservedState trial;
      trial.rho = star.rho;
      trial.mx = momentum[0];
      trial.my = momentum[1];
      trial.E = energy;
      eval.state[static_cast<std::size_t>(row)] = trial;

      const bool nonfinite = !std::isfinite(trial.rho) || !std::isfinite(trial.mx) ||
                             !std::isfinite(trial.my) || !std::isfinite(trial.E) ||
                             !std::isfinite(p[static_cast<std::size_t>(row)]) ||
                             !std::isfinite(kinetic) || !std::isfinite(internal) ||
                             !std::isfinite(eos_pressure) || !std::isfinite(residual) ||
                             !std::isfinite(div_flux);
      if (nonfinite) {
        ++eval.nonfinite_count;
        continue;
      }

      eval.density_min = std::min(eval.density_min, trial.rho);
      eval.solved_pressure_min = std::min(eval.solved_pressure_min, p[static_cast<std::size_t>(row)]);
      eval.solved_pressure_max = std::max(eval.solved_pressure_max, p[static_cast<std::size_t>(row)]);
      eval.eos_pressure_min = std::min(eval.eos_pressure_min, eos_pressure);
      eval.eos_pressure_max = std::max(eval.eos_pressure_max, eos_pressure);
      eval.internal_energy_min = std::min(eval.internal_energy_min, internal);
      eval.internal_energy_max = std::max(eval.internal_energy_max, internal);
      eval.momentum_correction_linf =
          std::max(eval.momentum_correction_linf,
                   std::sqrt((trial.mx - star.mx) * (trial.mx - star.mx) +
                             (trial.my - star.my) * (trial.my - star.my)));
      const amrex::Real energy_delta = std::abs(trial.E - star.E);
      eval.energy_correction_linf = std::max(eval.energy_correction_linf, energy_delta);
      eval.energy_closeout_delta_linf = std::max(eval.energy_closeout_delta_linf, energy_delta);
      eval.energy_flux_divergence_linf = std::max(eval.energy_flux_divergence_linf, std::abs(div_flux));

      if (trial.rho <= amrex::Real(0.0)) {
        ++eval.density_failure_count;
      }
      if (p[static_cast<std::size_t>(row)] <= amrex::Real(0.0) || eos_pressure <= amrex::Real(0.0)) {
        ++eval.pressure_failure_count;
      }
      if (internal <= amrex::Real(0.0)) {
        ++eval.internal_energy_failure_count;
      }
      residual_sum += std::abs(residual);
      mismatch_sum += std::abs(mismatch);
      energy_delta_sum += energy_delta;
      eval.residual_linf = std::max(eval.residual_linf, std::abs(residual));
      eval.residual_relative_linf =
          std::max(eval.residual_relative_linf, std::abs(residual) / residual_scale);
      eval.pressure_mismatch_linf = std::max(eval.pressure_mismatch_linf, std::abs(mismatch));
      eval.pressure_mismatch_relative_linf =
          std::max(eval.pressure_mismatch_relative_linf, std::abs(mismatch) / mismatch_scale);
      ++eval.residual_count;
    }
  }

  int eq78_count = 0;
  auto final_velocity_at = [&](int i, int j, int axis) {
    const SourceMapCellIndex mapped = source_map_pressure_index(grid, i, j);
    const ConservedState& state =
        eval.state[static_cast<std::size_t>(grid.row(mapped.i, mapped.j))];
    const amrex::Real momentum = axis == 0 ? state.mx : state.my;
    if (!std::isfinite(state.rho) || state.rho == amrex::Real(0.0) ||
        !std::isfinite(momentum)) {
      return std::numeric_limits<amrex::Real>::quiet_NaN();
    }
    return momentum / state.rho;
  };

  for (int j = 0; j < grid.ny; ++j) {
    for (int i = 0; i < grid.nx; ++i) {
      const amrex::Real u_e = final_velocity_at(i + 1, j, 0);
      const amrex::Real u_w = final_velocity_at(i - 1, j, 0);
      const amrex::Real v_n = final_velocity_at(i, j + 1, 1);
      const amrex::Real v_s = final_velocity_at(i, j - 1, 1);
      const amrex::Real div =
          amrex::Real(0.5) * (u_e - u_w) * inv_dx + amrex::Real(0.5) * (v_n - v_s) * inv_dy;
      const amrex::Real velocity_scale =
          std::max({std::abs(u_e), std::abs(u_w), std::abs(v_n), std::abs(v_s),
                    amrex::Real(1.0)});
      const amrex::Real scaled_denominator = velocity_scale * std::max(inv_dx, inv_dy);
      if (!std::isfinite(div) || !std::isfinite(scaled_denominator) ||
          scaled_denominator <= amrex::Real(0.0)) {
        continue;
      }
      eval.eq78_divergence_linf = std::max(eval.eq78_divergence_linf, std::abs(div));
      eval.eq78_divergence_scaled_linf =
          std::max(eval.eq78_divergence_scaled_linf, std::abs(div) / scaled_denominator);
      ++eq78_count;
    }
  }

  eval.positivity_count =
      eval.density_failure_count + eval.pressure_failure_count + eval.internal_energy_failure_count;
  if (eval.nonfinite_count > 0 || eval.residual_count == 0) {
    eval.residual_l1 = std::numeric_limits<amrex::Real>::infinity();
    eval.residual_linf = std::numeric_limits<amrex::Real>::infinity();
    eval.residual_relative_linf = std::numeric_limits<amrex::Real>::infinity();
    eval.pressure_mismatch_l1 = std::numeric_limits<amrex::Real>::infinity();
    eval.pressure_mismatch_linf = std::numeric_limits<amrex::Real>::infinity();
    eval.pressure_mismatch_relative_linf = std::numeric_limits<amrex::Real>::infinity();
    eval.energy_closeout_delta_l1 = std::numeric_limits<amrex::Real>::infinity();
    eval.eq78_divergence_linf = std::numeric_limits<amrex::Real>::infinity();
    eval.eq78_divergence_scaled_linf = std::numeric_limits<amrex::Real>::infinity();
  } else if (eq78_count == 0) {
    eval.eq78_divergence_linf = std::numeric_limits<amrex::Real>::infinity();
    eval.eq78_divergence_scaled_linf = std::numeric_limits<amrex::Real>::infinity();
  } else {
    const amrex::Real inv_count = amrex::Real(1.0) / static_cast<amrex::Real>(eval.residual_count);
    eval.residual_l1 = residual_sum * inv_count;
    eval.pressure_mismatch_l1 = mismatch_sum * inv_count;
    eval.energy_closeout_delta_l1 = energy_delta_sum * inv_count;
  }
  const auto finite_or_nan = [](amrex::Real value) {
    return std::isfinite(value) ? value : std::numeric_limits<amrex::Real>::quiet_NaN();
  };
  eval.density_min = finite_or_nan(eval.density_min);
  eval.solved_pressure_min = finite_or_nan(eval.solved_pressure_min);
  eval.solved_pressure_max = finite_or_nan(eval.solved_pressure_max);
  eval.eos_pressure_min = finite_or_nan(eval.eos_pressure_min);
  eval.eos_pressure_max = finite_or_nan(eval.eos_pressure_max);
  eval.internal_energy_min = finite_or_nan(eval.internal_energy_min);
  eval.internal_energy_max = finite_or_nan(eval.internal_energy_max);
  return eval;
}

void copy_source_map_eval_to_state(const SourceMapFinalEval& eval, amrex::MultiFab& state,
                                   const CoupledHostGrid& grid, const amrex::Geometry& geom,
                                   const RunConfig& cfg)
{
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        store_cons(arr, i, j, 0, eval.state[static_cast<std::size_t>(grid.row(i, j))]);
      }
    }
  }
  fill_problem_ghosts(state, geom, cfg);
}

amrex::Real source_map_star_pressure_mismatch_linf(const CoupledHostGrid& grid,
                                                   const std::vector<amrex::Real>& p,
                                                   amrex::Real gamma)
{
  const amrex::Real gamma_minus_one = gamma - amrex::Real(1.0);
  amrex::Real linf = amrex::Real(0.0);
  for (int j = 0; j < grid.ny; ++j) {
    for (int i = 0; i < grid.nx; ++i) {
      const int row = grid.row(i, j);
      const ConservedState& star = grid.at(i, j);
      const amrex::Real kinetic =
          amrex::Real(0.5) * (star.mx * star.mx + star.my * star.my) / star.rho;
      const amrex::Real p_eos = gamma_minus_one * (star.E - kinetic);
      if (!std::isfinite(p_eos) || !std::isfinite(p[static_cast<std::size_t>(row)])) {
        return std::numeric_limits<amrex::Real>::infinity();
      }
      linf = std::max(linf, std::abs(p[static_cast<std::size_t>(row)] - p_eos));
    }
  }
  return linf;
}

ImexStepDiagnostics euler_imex_bdltv20_o1_source_map_picard_one_step(
    amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg, amrex::Real dt,
    amrex::MultiFab* solved_pressure_check = nullptr)
{
  BL_PROFILE("euler_imex_bdltv20_o1_source_map_picard_one_step()");

  const bool t1s2 = is_bdltv20_t1_s2_source_map_form(cfg.imex_form);
  const std::string pressure_check_source =
      t1s2 ? "bdltv20_t1_s2_source_map_picard_host_solve"
           : "bdltv20_o1_source_map_picard_host_solve";
  ImexStepDiagnostics diagnostics;
  diagnostics.solver = "host_gmres_source_map_pressure_row_paper_gmres_aligned";
  diagnostics.solver_status = "not_run";
  diagnostics.solver_convergence_criterion = "linear_pressure_row_absolute_or_relative_per_picard_fixed_Rp";
  diagnostics.pressure_solver_path = "host_single_level_gmres";
  diagnostics.pressure_operator_label =
      "bdltv20_o1_source_picard_eq49_50_52_eq76_77_zdrop_2d_pressure_row_no_eq51_claim";
  diagnostics.pressure_row_stencil =
      "eq49_76_lhs_enthalpy_eq50_77_rhs_star_energy_lagged_kinetic_star_momentum_eq52_host_sparse_zdrop_2d_no_eq51_claim";
  diagnostics.pressure_energy_form = to_string(cfg.imex_form);
  diagnostics.pressure_energy_closeout_mode =
      "bdltv20_eq45_eq74_zdrop_2d_conservative_enthalpy_flux_closeout";
  diagnostics.pressure_energy_closeout_conservation_policy = "finite_volume_conservative_enthalpy_flux";
  diagnostics.density_denominator_floor_policy = "fail_if_predictor_density_le_floor";
  if (cfg.problem == ProblemKind::ShockDensityBubble2D) {
    diagnostics.pressure_boundary_policy =
        "x_left=fixed_post_shock_pressure_dirichlet;"
        "x_right=outflow_pressure_neumann;"
        "y_bottom=reflective_symmetry_pressure_neumann;"
        "y_top=outflow_pressure_neumann";
    diagnostics.boundary_flux_residual_policy =
        "shock_density_bubble_physical_boundary_states_for_pressure_row_face_enthalpy_and_momentum";
  } else {
    diagnostics.pressure_boundary_policy =
        std::string("x=") + (geom.isPeriodic(0) ? "periodic_row_wrap"
                                                : "bdltv20_even_neumann_pressure_reflection") +
        ";y=" + (geom.isPeriodic(1) ? "periodic_row_wrap"
                                    : "bdltv20_even_neumann_pressure_reflection");
    diagnostics.boundary_flux_residual_policy =
        "source_map_pressure_periodic_wrap_or_bdltv20_even_neumann_reflection";
  }
  diagnostics.pressure_retry_rule = "none";
  diagnostics.timestep_rescue_rule = "none";
  diagnostics.kinetic_lag_mode = "fixed_Rp_eq54_eq71_72_corrected_momentum_picard_refresh";
  diagnostics.kinetic_lag_source =
      "K0_eq53_star_momentum_Kr_eq54_eq71_72_centered_pressure_gradient_zdrop_2d_epsilon1";
  diagnostics.pressure_energy_coupling_mode =
      "bdltv20_eq49_50_52_eq76_77_zdrop_2d_pressure_row_coupling_no_eq51_claim";
  diagnostics.kinetic_jacobian_mode = "source_picard_without_local_kinetic_jacobian_tuning";
  diagnostics.face_enthalpy_policy =
      "bdltv20_eq46_corrected_momentum_weighted_face_enthalpy_with_tolerance_guard";
  diagnostics.face_momentum_policy =
      "centered_eq54_eq71_72_corrected_normal_momentum";
  diagnostics.pressure_row_support_mode = "source_map_pressure_row_terms_and_face_fluxes";
  diagnostics.pressure_row_support_radius = 2;
  diagnostics.failed_row_face_flux_entries_count = 4;
  diagnostics.source_map_version =
      t1s2 ? "bdltv20_t1_s2_section3p2_eq43_to_55_section3p4_zdrop_eq69_to_79_section4p2_eq87_to_89_conserved_minmod"
           : "bdltv20_o1_section3p2_eq43_to_55_section3p4_zdrop_eq69_to_79";
  diagnostics.source_primary_paper = "bdltv20_2020_jcp_415_109486";
  diagnostics.source_equation_set =
      t1s2 ? "eq25_28_split_eq43_55_source_map_eq69_79_zdrop_2d_plus_eq87_89_conserved_minmod_explicit_advection"
           : "eq25_28_split_eq43_55_source_map_eq69_79_zdrop_2d";
  diagnostics.source_time_order = "first_order_time_t1";
  diagnostics.source_space_order = t1s2 ? "second_order_space_s2" : "first_order_space_s1";
  diagnostics.source_advection_reconstruction =
      t1s2 ? "bdltv20_eq87_89_conserved_variable_minmod_explicit_toro_vazquez_star_update"
           : "piecewise_constant_first_order_source_map_star_update";
  diagnostics.source_scaling = "epsilon_1_dimensional_no_formal_ap_claim";
  diagnostics.source_pressure_row_normalization = "normalized_by_cell_measure";
  diagnostics.source_pressure_row_boundary_policy =
      cfg.problem == ProblemKind::ShockDensityBubble2D
          ? "shock_density_bubble_left_postshock_dirichlet_right_top_outflow_bottom_reflective"
          : "periodic_wrap_or_bdltv20_even_neumann_reflection";
  diagnostics.source_picard_count_source = "euler.imex_picard_iterations";
  diagnostics.source_picard_initialization = "eq53_star_momentum_old_pressure_initial_guess";
  diagnostics.source_picard_iterations_config = cfg.imex_picard_iterations;
  diagnostics.pressure_picard_convergence_policy =
      "fixed_Rp_eq49_55_pressure_solve_then_eq54_momentum_refresh_then_"
      "eq55_accept_final_iterate";
  diagnostics.eq78_divergence_diagnostic =
      "bdltv20_eq78_zdrop_2d_centered_final_corrected_velocity_diagnostic_only_no_ap_proof_no_eq51_claim";
  diagnostics.eq78_divergence_boundary_policy =
      "periodic_wrap_or_clamped_state_extension_approx_no_eq51_claim";
  diagnostics.eq78_divergence_scaled_normalization =
      "abs_divergence_over_max_neighbor_velocity_component_unit_floor_times_max_inverse_cell_width";
  diagnostics.pressure_stabilization = "off";
  diagnostics.pressure_stabilization_trigger_count = 0;
  diagnostics.density_denominator_floor_value = amrex::Real(1.0e-14);
  diagnostics.picard_extra_iterations_requested = cfg.imex_picard_extra_iterations;
  diagnostics.picard_extra_iterations_schema =
      "noncontrolling_for_source_map_fixed_Rp_uses_"
      "euler_imex_picard_iterations";
  diagnostics.solver_requested_relative_tol = cfg.imex_solver_tol;
  diagnostics.solver_requested_absolute_tol = cfg.imex_solver_tol;
  configure_pressure_energy_check_diagnostics(diagnostics, cfg, "not_captured_before_pressure_solve",
                                              pressure_check_source);

  if (amrex::ParallelDescriptor::NProcs() != 1) {
    diagnostics.solver_status = "unsupported_mpi_for_host_single_level_gmres";
    return diagnostics;
  }

  const amrex::BoxArray& ba = state.boxArray();
  const amrex::DistributionMapping& dm = state.DistributionMap();
  if (cfg.problem == ProblemKind::ShockDensityBubble2D) {
    fill_problem_ghosts(state, geom, cfg);
  } else {
    fill_physical_outflow_ghosts(state, geom);
  }
  amrex::MultiFab star(ba, dm, NCons, 2);
  compute_bdltv20_pressure_split_star(state, star, geom, cfg, dt);
  const CoupledHostGrid grid = load_coupled_host_grid(star, geom, cfg);
  const int cell_count = grid.nx * grid.ny;
  const auto dx = geom.CellSizeArray();
  const amrex::Real inv_dx = amrex::Real(1.0) / dx[0];
  const amrex::Real inv_dy = amrex::Real(1.0) / dx[1];
  const amrex::Real density_floor = diagnostics.density_denominator_floor_value;
  std::vector<SourceMapMomentum> lagged_momentum = source_map_eq53_star_momentum_field(grid);
  diagnostics.kinetic_lag_nonfinite_count =
      source_map_momentum_nonfinite_count(lagged_momentum);

  std::vector<amrex::Real> pressure(static_cast<std::size_t>(cell_count), amrex::Real(0.0));
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        pressure[static_cast<std::size_t>(grid.row(i, j))] =
            to_primitive(load_cons(arr, i, j, 0), cfg.gamma).p;
      }
    }
  }

  diagnostics.trial_density_min = std::numeric_limits<amrex::Real>::infinity();
  for (const ConservedState& cell : grid.cell) {
    if (!std::isfinite(cell.rho) || cell.rho <= density_floor) {
      ++diagnostics.density_denominator_floor_count;
    }
    if (std::isfinite(cell.rho)) {
      diagnostics.trial_density_min = std::min(diagnostics.trial_density_min, cell.rho);
    }
  }
  diagnostics.floor_used_count = diagnostics.density_denominator_floor_count;
  if (diagnostics.density_denominator_floor_count > 0) {
    diagnostics.solver_status = "density_denominator_floor_used";
    diagnostics.trial_positivity_count = diagnostics.density_denominator_floor_count;
    diagnostics.positivity_failure_count = diagnostics.density_denominator_floor_count;
    return diagnostics;
  }

  diagnostics.pressure_energy_residual_before_linf =
      source_map_star_pressure_mismatch_linf(grid, pressure, cfg.gamma);
  auto copy_pressure_check = [&]() {
    if (solved_pressure_check != nullptr) {
      copy_pressure_vector_to_multifab(pressure, *solved_pressure_check, grid, geom);
      configure_pressure_energy_check_diagnostics(diagnostics, cfg, "captured",
                                                  pressure_check_source);
    }
  };

  for (int picard = 0; picard < cfg.imex_picard_iterations; ++picard) {
    const SourceMapPressureSystem system =
        source_map_build_pressure_system(grid, pressure, lagged_momentum, cfg, geom, dt);
    diagnostics.pressure_matrix_rows = system.matrix.rows();
    diagnostics.pressure_matrix_nonzeros = system.matrix.nonZeros();
    diagnostics.pressure_matrix_max_row_width = 0;
    diagnostics.kinetic_jacobian_abs_max = amrex::Real(0.0);
    for (int outer = 0; outer < system.matrix.outerSize(); ++outer) {
      int row_width = 0;
      for (HostSparseMatrix::InnerIterator it(system.matrix, outer); it; ++it) {
        ++row_width;
        diagnostics.kinetic_jacobian_abs_max =
            std::max(diagnostics.kinetic_jacobian_abs_max, std::abs(it.value()));
      }
      diagnostics.pressure_matrix_max_row_width =
          std::max(diagnostics.pressure_matrix_max_row_width, row_width);
    }
    diagnostics.failed_row_sparse_entries_count =
        std::max(diagnostics.failed_row_sparse_entries_count, diagnostics.pressure_matrix_max_row_width);
    const HostSparseMatrix symmetry_diff = system.matrix - HostSparseMatrix(system.matrix.transpose());
    amrex::Real symmetry_linf = amrex::Real(0.0);
    for (int outer = 0; outer < symmetry_diff.outerSize(); ++outer) {
      for (HostSparseMatrix::InnerIterator it(symmetry_diff, outer); it; ++it) {
        symmetry_linf = std::max(symmetry_linf, std::abs(it.value()));
      }
    }
    diagnostics.pressure_matrix_symmetry_linf = symmetry_linf;

    if (system.nonfinite_count > 0 || !system.rhs.allFinite()) {
      diagnostics.solver_status = "source_map_pressure_row_nonfinite";
      diagnostics.kinetic_lag_nonfinite_count = system.nonfinite_count;
      diagnostics.nonfinite_count = system.nonfinite_count;
      copy_pressure_check();
      return diagnostics;
    }

    Eigen::VectorXd guess(cell_count);
    for (int row = 0; row < cell_count; ++row) {
      guess(row) = pressure[static_cast<std::size_t>(row)];
    }
    const Eigen::VectorXd initial_linear_residual = system.matrix * guess - system.rhs;
    diagnostics.solver_initial_residual = initial_linear_residual.lpNorm<Eigen::Infinity>();
    diagnostics.solver_initial_rhs_norm =
        std::max(system.rhs.lpNorm<Eigen::Infinity>(), amrex::Real(1.0));

    const HostGmresResult gmres =
        solve_host_sparse_gmres(system.matrix, system.rhs, guess, cfg.imex_solver_max_iter,
                                cfg.imex_solver_tol);
    const Eigen::VectorXd solved = gmres.solution;
    diagnostics.solver_iterations += gmres.iterations;
    const Eigen::VectorXd linear_residual = system.matrix * solved - system.rhs;
    diagnostics.pressure_linear_residual_linf = linear_residual.lpNorm<Eigen::Infinity>();
    diagnostics.pressure_linear_residual_relative_linf =
        diagnostics.pressure_linear_residual_linf / diagnostics.solver_initial_rhs_norm;
    diagnostics.solver_residual = diagnostics.pressure_linear_residual_linf;
    diagnostics.solver_relative_residual = diagnostics.pressure_linear_residual_relative_linf;
    const bool abs_ok = diagnostics.pressure_linear_residual_linf <= cfg.imex_solver_tol;
    const bool rel_ok = diagnostics.pressure_linear_residual_relative_linf <= cfg.imex_solver_tol;
    diagnostics.linear_solver_status_last =
        abs_ok ? "converged_absolute" : (rel_ok ? "converged_relative" : "linear_solve_not_converged");
    if (!std::isfinite(diagnostics.pressure_linear_residual_linf) || !solved.allFinite()) {
      diagnostics.solver_status = "linear_residual_nonfinite";
      diagnostics.nonfinite_count = 1;
      copy_pressure_check();
      return diagnostics;
    }
    if (!abs_ok && !rel_ok) {
      diagnostics.solver_status =
          gmres.breakdown ? "gmres_breakdown_not_converged" : "linear_solve_not_converged";
      copy_pressure_check();
      return diagnostics;
    }

    amrex::Real delta_linf = amrex::Real(0.0);
    amrex::Real delta_l2_sum = amrex::Real(0.0);
    amrex::Real pmax = amrex::Real(1.0);
    for (int row = 0; row < cell_count; ++row) {
      const amrex::Real old_p = pressure[static_cast<std::size_t>(row)];
      const amrex::Real new_p = solved(row);
      const amrex::Real delta = new_p - old_p;
      delta_linf = std::max(delta_linf, std::abs(delta));
      delta_l2_sum += delta * delta;
      pmax = std::max({pmax, std::abs(old_p), std::abs(new_p)});
      pressure[static_cast<std::size_t>(row)] = new_p;
    }
    lagged_momentum =
        source_map_eq54_corrected_momentum_field(grid, pressure, dt, inv_dx, inv_dy, &cfg);
    diagnostics.kinetic_lag_nonfinite_count =
        source_map_momentum_nonfinite_count(lagged_momentum);
    if (diagnostics.kinetic_lag_nonfinite_count > 0) {
      diagnostics.solver_status = "source_map_lagged_momentum_nonfinite";
      diagnostics.nonfinite_count = diagnostics.kinetic_lag_nonfinite_count;
      copy_pressure_check();
      return diagnostics;
    }
    ++diagnostics.kinetic_lag_update_count;
    diagnostics.pressure_picard_delta_linf = delta_linf;
    diagnostics.source_picard_delta_linf = delta_linf;
    diagnostics.source_picard_delta_l2_relative_to_pmax =
        std::sqrt(delta_l2_sum / static_cast<amrex::Real>(cell_count)) / pmax;
    ++diagnostics.source_picard_iterations_used;
    ++diagnostics.picard_pressure_solves_used;
  }

  diagnostics.solver_status = diagnostics.linear_solver_status_last;
  const SourceMapFinalEval final_eval =
      evaluate_source_map_final(grid, pressure, cfg, geom, dt);
  diagnostics.coupled_residual_l1_final = final_eval.residual_l1;
  diagnostics.coupled_residual_linf_final = final_eval.residual_linf;
  diagnostics.coupled_residual_relative_linf_final = final_eval.residual_relative_linf;
  diagnostics.pressure_energy_residual_l1 = final_eval.pressure_mismatch_l1;
  diagnostics.pressure_energy_residual_linf = final_eval.pressure_mismatch_linf;
  diagnostics.pressure_energy_residual_relative_linf = final_eval.pressure_mismatch_relative_linf;
  diagnostics.pressure_energy_residual_after_linf = final_eval.pressure_mismatch_linf;
  diagnostics.pressure_energy_residual_run_max_linf = final_eval.pressure_mismatch_linf;
  diagnostics.pressure_energy_residual_reduction_ratio =
      std::isfinite(diagnostics.pressure_energy_residual_before_linf) &&
              diagnostics.pressure_energy_residual_before_linf > amrex::Real(0.0)
          ? final_eval.pressure_mismatch_linf / diagnostics.pressure_energy_residual_before_linf
          : std::numeric_limits<amrex::Real>::quiet_NaN();
  diagnostics.eos_pressure_residual_l1 = final_eval.pressure_mismatch_l1;
  diagnostics.eos_pressure_residual_linf = final_eval.pressure_mismatch_linf;
  diagnostics.rho_min = final_eval.density_min;
  diagnostics.pressure_min = final_eval.eos_pressure_min;
  diagnostics.internal_energy_min = final_eval.internal_energy_min;
  diagnostics.trial_density_min = final_eval.density_min;
  diagnostics.trial_pressure_min = final_eval.eos_pressure_min;
  diagnostics.trial_internal_energy_min = final_eval.internal_energy_min;
  diagnostics.solved_pressure_min = final_eval.solved_pressure_min;
  diagnostics.solved_pressure_max = final_eval.solved_pressure_max;
  diagnostics.eos_pressure_min = final_eval.eos_pressure_min;
  diagnostics.eos_pressure_max = final_eval.eos_pressure_max;
  diagnostics.internal_energy_max = final_eval.internal_energy_max;
  diagnostics.trial_nonfinite_count = final_eval.nonfinite_count;
  diagnostics.trial_positivity_count = final_eval.positivity_count;
  diagnostics.nonfinite_count = final_eval.nonfinite_count;
  diagnostics.positivity_failure_count = final_eval.positivity_count;
  diagnostics.momentum_correction_linf = final_eval.momentum_correction_linf;
  diagnostics.energy_correction_linf = final_eval.energy_correction_linf;
  diagnostics.energy_flux_linf = final_eval.energy_flux_linf;
  diagnostics.energy_flux_divergence_linf = final_eval.energy_flux_divergence_linf;
  diagnostics.eq78_divergence_linf = final_eval.eq78_divergence_linf;
  diagnostics.eq78_divergence_scaled_linf = final_eval.eq78_divergence_scaled_linf;
  diagnostics.energy_flux_closeout_consistency_linf = amrex::Real(0.0);
  diagnostics.energy_closeout_delta_l1 = final_eval.energy_closeout_delta_l1;
  diagnostics.energy_closeout_delta_linf = final_eval.energy_closeout_delta_linf;
  diagnostics.face_enthalpy_denominator_guard_count = final_eval.face_enthalpy_guard_count;
  diagnostics.face_enthalpy_denominator_guard_count_max = final_eval.face_enthalpy_guard_count;
  diagnostics.face_enthalpy_nonfinite_count = final_eval.face_enthalpy_nonfinite_count;
  diagnostics.energy_flux_nonfinite_count = final_eval.energy_flux_nonfinite_count;
  copy_pressure_check();
  copy_source_map_eval_to_state(final_eval, state, grid, geom, cfg);
  return diagnostics;
}

amrex::MultiFab make_plot_state(const amrex::MultiFab& state, const RunConfig& cfg)
{
  amrex::MultiFab plot(state.boxArray(), state.DistributionMap(), 5, 0);
  const amrex::Real gamma = cfg.gamma;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const src = state.const_array(mfi);
    auto const dst = plot.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real rho = src(i, j, k, Rho);
      const amrex::Real inv_rho = amrex::Real(1.0) / rho;
      const amrex::Real ux = src(i, j, k, Mx) * inv_rho;
      const amrex::Real uy = src(i, j, k, My) * inv_rho;
      const amrex::Real kinetic = amrex::Real(0.5) * rho * (ux * ux + uy * uy);
      const amrex::Real pressure = (gamma - amrex::Real(1.0)) * (src(i, j, k, E) - kinetic);
      dst(i, j, k, 0) = rho;
      dst(i, j, k, 1) = src(i, j, k, Mx);
      dst(i, j, k, 2) = src(i, j, k, My);
      dst(i, j, k, 3) = src(i, j, k, E);
      dst(i, j, k, 4) = pressure;
    });
  }

  return plot;
}
