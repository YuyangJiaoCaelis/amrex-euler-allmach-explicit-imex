// Diagnostics, timestep selection, and problem metrics.
amrex::Real compute_blob_dt(const amrex::Geometry& geom, const RunConfig& cfg, amrex::Real time)
{
  const auto dx = geom.CellSizeArray();
  amrex::Real dt = std::numeric_limits<amrex::Real>::max();
  if (std::abs(cfg.velocity_x) > amrex::Real(0.0)) {
    dt = std::min(dt, cfg.cfl * dx[0] / std::abs(cfg.velocity_x));
  }
  if (std::abs(cfg.velocity_y) > amrex::Real(0.0)) {
    dt = std::min(dt, cfg.cfl * dx[1] / std::abs(cfg.velocity_y));
  }
  if (dt == std::numeric_limits<amrex::Real>::max()) {
    dt = cfg.stop_time - time;
  }
  return std::min(dt, cfg.stop_time - time);
}

Diagnostics compute_diagnostics(const amrex::MultiFab& state, const RunConfig& cfg)
{
  amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin, amrex::ReduceOpSum> reduce_op;
  amrex::ReduceData<amrex::Real, amrex::Real, int> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const amrex::Real gamma = cfg.gamma;
  const amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    reduce_op.eval(box, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
      const ConservedState u = load_cons(arr, i, j, k);
      const PrimitiveState q = to_primitive(u, gamma);
      const bool bad = !std::isfinite(q.rho) || !std::isfinite(q.u) || !std::isfinite(q.v) || !std::isfinite(q.p) ||
                       !std::isfinite(u.E);
      const amrex::Real rho_min = bad ? -inf : q.rho;
      const amrex::Real p_min = bad ? -inf : q.p;
      return {rho_min, p_min, bad ? 1 : 0};
    });
  }

  ReduceTuple values = reduce_data.value();
  Diagnostics diagnostics{amrex::get<0>(values), amrex::get<1>(values), amrex::get<2>(values)};
  amrex::ParallelDescriptor::ReduceRealMin(diagnostics.rho_min);
  amrex::ParallelDescriptor::ReduceRealMin(diagnostics.pressure_min);
  amrex::ParallelDescriptor::ReduceIntSum(diagnostics.nonfinite_count);
  return diagnostics;
}

AdmissibilityDiagnostics compute_admissibility_diagnostics(const amrex::MultiFab& state, const RunConfig& cfg)
{
  amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin, amrex::ReduceOpMin, amrex::ReduceOpSum,
                   amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum>
      reduce_op;
  amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, int, int, int, int> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const amrex::Real gamma = cfg.gamma;
  const amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    reduce_op.eval(box, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
      const ConservedState u = load_cons(arr, i, j, k);
      const bool conserved_nonfinite = !std::isfinite(u.rho) || !std::isfinite(u.mx) ||
                                       !std::isfinite(u.my) || !std::isfinite(u.E);
      const bool density_bad = !conserved_nonfinite && u.rho <= amrex::Real(0.0);

      amrex::Real internal_energy = -inf;
      amrex::Real pressure = -inf;
      bool thermodynamic_nonfinite = false;
      if (!conserved_nonfinite && !density_bad) {
        const amrex::Real kinetic = amrex::Real(0.5) * (u.mx * u.mx + u.my * u.my) / u.rho;
        internal_energy = u.E - kinetic;
        pressure = (gamma - amrex::Real(1.0)) * internal_energy;
        thermodynamic_nonfinite = !std::isfinite(kinetic) || !std::isfinite(internal_energy) ||
                                  !std::isfinite(pressure);
      }

      const bool nonfinite = conserved_nonfinite || thermodynamic_nonfinite;
      const bool pressure_bad = !nonfinite && !density_bad && pressure <= amrex::Real(0.0);
      const bool internal_energy_bad =
          !nonfinite && !density_bad && internal_energy <= amrex::Real(0.0);

      return {nonfinite ? -inf : u.rho,
              (nonfinite || density_bad) ? -inf : pressure,
              (nonfinite || density_bad) ? -inf : internal_energy,
              nonfinite ? 1 : 0,
              density_bad ? 1 : 0,
              pressure_bad ? 1 : 0,
              internal_energy_bad ? 1 : 0};
    });
  }

  ReduceTuple values = reduce_data.value();
  AdmissibilityDiagnostics diagnostics{amrex::get<0>(values),
                                       amrex::get<1>(values),
                                       amrex::get<2>(values),
                                       amrex::get<3>(values),
                                       amrex::get<4>(values),
                                       amrex::get<5>(values),
                                       amrex::get<6>(values)};
  amrex::ParallelDescriptor::ReduceRealMin(diagnostics.rho_min);
  amrex::ParallelDescriptor::ReduceRealMin(diagnostics.pressure_min);
  amrex::ParallelDescriptor::ReduceRealMin(diagnostics.internal_energy_min);
  amrex::ParallelDescriptor::ReduceIntSum(diagnostics.nonfinite_count);
  amrex::ParallelDescriptor::ReduceIntSum(diagnostics.density_failure_count);
  amrex::ParallelDescriptor::ReduceIntSum(diagnostics.pressure_failure_count);
  amrex::ParallelDescriptor::ReduceIntSum(diagnostics.internal_energy_failure_count);
  return diagnostics;
}

bool admissible(const AdmissibilityDiagnostics& diagnostics)
{
  return diagnostics.nonfinite_count == 0 && diagnostics.density_failure_count == 0 &&
         diagnostics.pressure_failure_count == 0 && diagnostics.internal_energy_failure_count == 0 &&
         diagnostics.rho_min > amrex::Real(0.0) && diagnostics.pressure_min > amrex::Real(0.0) &&
         diagnostics.internal_energy_min > amrex::Real(0.0);
}

std::string failure_category_from_admissibility(const AdmissibilityDiagnostics& diagnostics)
{
  if (diagnostics.nonfinite_count > 0) {
    return "nonfinite_state";
  }
  if (diagnostics.density_failure_count > 0 || diagnostics.rho_min <= amrex::Real(0.0)) {
    return "nonpositive_density";
  }
  if (diagnostics.pressure_failure_count > 0 || diagnostics.pressure_min <= amrex::Real(0.0)) {
    return "nonpositive_pressure";
  }
  if (diagnostics.internal_energy_failure_count > 0 ||
      diagnostics.internal_energy_min <= amrex::Real(0.0)) {
    return "nonpositive_internal_energy";
  }
  return "none";
}

void update_high_trial_run(HighTrialRunDiagnostics& run, const AdmissibilityDiagnostics& step)
{
  run.rho_min = std::min(run.rho_min, step.rho_min);
  run.pressure_min = std::min(run.pressure_min, step.pressure_min);
  run.internal_energy_min = std::min(run.internal_energy_min, step.internal_energy_min);
  run.nonfinite_step_count += step.nonfinite_count > 0 ? 1 : 0;
  run.density_failure_step_count += step.density_failure_count > 0 ? 1 : 0;
  run.pressure_failure_step_count += step.pressure_failure_count > 0 ? 1 : 0;
  run.internal_energy_failure_step_count += step.internal_energy_failure_count > 0 ? 1 : 0;
}

GreshoMetrics compute_gresho_metrics(const amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg)
{
  if (cfg.problem != ProblemKind::GreshoVortex) {
    return GreshoMetrics{};
  }

  amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
                   amrex::ReduceOpSum, amrex::ReduceOpSum>
      reduce_op;
  amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real>
      reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real p0 = cfg.pressure;
  const amrex::Real cell_vol = cell_volume(geom);

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    reduce_op.eval(box, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
      const amrex::Real x = plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
      const amrex::Real y = plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
      const PrimitiveState numerical = to_primitive(load_cons(arr, i, j, k), gamma);
      const PrimitiveState exact = gresho_state(x, y, cfg);
      const amrex::Real du = numerical.u - exact.u;
      const amrex::Real dv = numerical.v - exact.v;
      const amrex::Real pressure_perturbation_error = std::abs((numerical.p - p0) - (exact.p - p0));
      return {std::abs(numerical.rho - exact.rho),
              std::sqrt(du * du + dv * dv),
              std::abs(numerical.p - exact.p),
              pressure_perturbation_error,
              std::abs(exact.p - p0),
              kinetic_energy_density(exact) * cell_vol,
              kinetic_energy_density(numerical) * cell_vol};
    });
  }

  ReduceTuple values = reduce_data.value();
  amrex::Real density_l1 = amrex::get<0>(values);
  amrex::Real velocity_l1 = amrex::get<1>(values);
  amrex::Real pressure_l1 = amrex::get<2>(values);
  amrex::Real pressure_perturbation_l1 = amrex::get<3>(values);
  amrex::Real pressure_perturbation_scale = amrex::get<4>(values);
  amrex::Real kinetic_initial = amrex::get<5>(values);
  amrex::Real kinetic_final = amrex::get<6>(values);
  amrex::ParallelDescriptor::ReduceRealSum(density_l1);
  amrex::ParallelDescriptor::ReduceRealSum(velocity_l1);
  amrex::ParallelDescriptor::ReduceRealSum(pressure_l1);
  amrex::ParallelDescriptor::ReduceRealSum(pressure_perturbation_l1);
  amrex::ParallelDescriptor::ReduceRealSum(pressure_perturbation_scale);
  amrex::ParallelDescriptor::ReduceRealSum(kinetic_initial);
  amrex::ParallelDescriptor::ReduceRealSum(kinetic_final);

  const amrex::Real inv_count = amrex::Real(1.0) /
                                static_cast<amrex::Real>(cfg.n_cell[0] * cfg.n_cell[1]);
  GreshoMetrics metrics;
  metrics.density_l1_error = density_l1 * inv_count;
  metrics.velocity_l1_error = velocity_l1 * inv_count;
  metrics.pressure_l1_error = pressure_l1 * inv_count;
  metrics.pressure_perturbation_l1_error = pressure_perturbation_l1 * inv_count;
  metrics.pressure_perturbation_l1_relative_error =
      pressure_perturbation_scale > amrex::Real(1.0e-30)
          ? pressure_perturbation_l1 / pressure_perturbation_scale
          : std::numeric_limits<amrex::Real>::quiet_NaN();
  metrics.kinetic_energy_initial = kinetic_initial;
  metrics.kinetic_energy_final = kinetic_final;
  metrics.kinetic_energy_ratio =
      kinetic_initial > amrex::Real(0.0) ? kinetic_final / kinetic_initial
                                         : std::numeric_limits<amrex::Real>::quiet_NaN();
  metrics.reference_sound_speed = reference_sound_speed(cfg);
  metrics.velocity_scale = velocity_scale_from_mach(cfg);
  return metrics;
}

amrex::Real compute_euler_dt(const amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg,
                             amrex::Real time)
{
  amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
  amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const auto dx = geom.CellSizeArray();
  const amrex::Real gamma = cfg.gamma;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    reduce_op.eval(box, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
      const PrimitiveState q = to_primitive(load_cons(arr, i, j, k), gamma);
      const amrex::Real c = sound_speed(q, gamma);
      const amrex::Real spectral_rate = (std::abs(q.u) + c) / dx[0] + (std::abs(q.v) + c) / dx[1];
      return {spectral_rate};
    });
  }

  ReduceTuple values = reduce_data.value();
  amrex::Real max_rate = amrex::get<0>(values);
  amrex::ParallelDescriptor::ReduceRealMax(max_rate);

  amrex::Real dt = cfg.stop_time - time;
  if (max_rate > amrex::Real(0.0)) {
    dt = cfg.cfl / max_rate;
  }
  if (!std::isfinite(dt) || dt <= amrex::Real(0.0)) {
    dt = cfg.stop_time - time;
  }
  return std::min(dt, cfg.stop_time - time);
}

bool use_bdltv20_eq79_timestep(const RunConfig& cfg)
{
  if (has_bdltv20_o1_eq79_selector_controls(cfg)) {
    return true;
  }
  return cfg.method == MethodKind::Imex &&
         is_bdltv20_t1_s2_source_map_form(cfg.imex_form) &&
         has_base_bdltv20_paper_controls(cfg);
}

amrex::Real compute_imex_dt(const amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg,
                            amrex::Real time, ImexTimeStepDiagnostics& diagnostics)
{
  amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_op;
  amrex::ReduceData<amrex::Real, amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const auto dx = geom.CellSizeArray();
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real min_dx = std::min(dx[0], dx[1]);
  const bool bdltv20_eq79_timestep = use_bdltv20_eq79_timestep(cfg);

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    reduce_op.eval(box, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
      const PrimitiveState q = to_primitive(load_cons(arr, i, j, k), gamma);
      const amrex::Real c = sound_speed(q, gamma);
      const amrex::Real speed_norm = std::sqrt(q.u * q.u + q.v * q.v);
      const amrex::Real material_rate =
          bdltv20_eq79_timestep ? (amrex::Real(2.0) * speed_norm / min_dx)
                                : (std::abs(q.u) / dx[0] + std::abs(q.v) / dx[1]);
      const amrex::Real acoustic_rate = (std::abs(q.u) + c) / dx[0] + (std::abs(q.v) + c) / dx[1];
      return {material_rate, acoustic_rate};
    });
  }

  ReduceTuple values = reduce_data.value();
  diagnostics.material_rate_max = amrex::get<0>(values);
  diagnostics.acoustic_rate_max = amrex::get<1>(values);
  amrex::ParallelDescriptor::ReduceRealMax(diagnostics.material_rate_max);
  amrex::ParallelDescriptor::ReduceRealMax(diagnostics.acoustic_rate_max);

  amrex::Real max_rate = diagnostics.material_rate_max;
  if (cfg.imex_acoustic_startup != 0 && diagnostics.material_rate_max <= amrex::Real(1.0e-12) &&
      diagnostics.acoustic_rate_max > amrex::Real(0.0)) {
    max_rate = diagnostics.acoustic_rate_max;
    diagnostics.acoustic_startup_used = 1;
  }
  amrex::Real dt = cfg.stop_time - time;
  if (max_rate > amrex::Real(0.0)) {
    dt = cfg.imex_cfl / max_rate;
  }
  if (cfg.imex_acoustic_cfl_cap > amrex::Real(0.0) && diagnostics.acoustic_rate_max > amrex::Real(0.0)) {
    const amrex::Real acoustic_limited_dt = cfg.imex_acoustic_cfl_cap / diagnostics.acoustic_rate_max;
    if (acoustic_limited_dt < dt) {
      dt = acoustic_limited_dt;
      diagnostics.acoustic_cap_used = 1;
    }
  }
  if (!std::isfinite(dt) || dt <= amrex::Real(0.0)) {
    dt = cfg.stop_time - time;
  }
  return std::min(dt, cfg.stop_time - time);
}
