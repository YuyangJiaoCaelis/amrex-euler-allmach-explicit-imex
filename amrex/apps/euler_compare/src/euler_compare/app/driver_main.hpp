// Application driver, time integration loop, and report-facing status output.
namespace euler_compare {

int run_application(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  int return_code = 0;
  {
    BL_PROFILE("main()");
    const RunConfig cfg = read_config();

    if (bdltv20_paper_t1_s2_requested(cfg)) {
      return_code = run_bdltv20_paper_t1_s2(cfg);
      goto finalize_amrex;
    }

    const amrex::Geometry geom = make_geometry(cfg);
    amrex::BoxArray ba(geom.Domain());
    ba.maxSize(cfg.max_grid_size);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab state(ba, dm, NCons, 2);

    initialize_state(state, geom, cfg);
    amrex::Gpu::streamSynchronize();

    const amrex::Real volume = cell_volume(geom);
    const amrex::Real initial_mass = state.sum(Rho) * volume;
    const amrex::Real initial_energy = state.sum(E) * volume;

    int step = 0;
    int hllc_degenerate_star_count = 0;
    int imex_step_count = 0;
    int imex_acoustic_startup_steps = 0;
    int imex_acoustic_cap_steps = 0;
    int steps_rejected = 0;
    amrex::Real time = 0.0;
    amrex::Real imex_timestep_material_rate_max = amrex::Real(0.0);
    amrex::Real imex_timestep_acoustic_rate_max = amrex::Real(0.0);
    amrex::Real dt_min_run = std::numeric_limits<amrex::Real>::infinity();
    amrex::Real dt_max_run = amrex::Real(0.0);
    amrex::Real dt_sum_run = amrex::Real(0.0);
    int dt_count_run = 0;

    ImexStepDiagnostics last_imex_diagnostics;
    ImexStepDiagnostics max_imex_diagnostics;
    AdmissibilityDiagnostics imex_admissibility_last;
    HighTrialRunDiagnostics imex_admissibility_run;

    const auto max_real = [](amrex::Real current, amrex::Real trial) {
      if (!std::isfinite(current)) {
        return trial;
      }
      if (!std::isfinite(trial)) {
        return current;
      }
      return std::max(current, trial);
    };
    const auto update_imex_maxima = [&](const ImexStepDiagnostics& current) {
      max_imex_diagnostics.solver_iterations =
          std::max(max_imex_diagnostics.solver_iterations, current.solver_iterations);
      max_imex_diagnostics.solver_residual =
          max_real(max_imex_diagnostics.solver_residual, current.solver_residual);
      max_imex_diagnostics.solver_initial_residual =
          max_real(max_imex_diagnostics.solver_initial_residual,
                   current.solver_initial_residual);
      max_imex_diagnostics.solver_initial_rhs_norm =
          max_real(max_imex_diagnostics.solver_initial_rhs_norm,
                   current.solver_initial_rhs_norm);
      max_imex_diagnostics.solver_relative_residual =
          max_real(max_imex_diagnostics.solver_relative_residual,
                   current.solver_relative_residual);
      max_imex_diagnostics.pressure_matrix_rows =
          std::max(max_imex_diagnostics.pressure_matrix_rows,
                   current.pressure_matrix_rows);
      max_imex_diagnostics.pressure_matrix_nonzeros =
          std::max(max_imex_diagnostics.pressure_matrix_nonzeros,
                   current.pressure_matrix_nonzeros);
      max_imex_diagnostics.pressure_matrix_max_row_width =
          std::max(max_imex_diagnostics.pressure_matrix_max_row_width,
                   current.pressure_matrix_max_row_width);
      max_imex_diagnostics.pressure_matrix_symmetry_linf =
          max_real(max_imex_diagnostics.pressure_matrix_symmetry_linf,
                   current.pressure_matrix_symmetry_linf);
      max_imex_diagnostics.pressure_linear_residual_linf =
          max_real(max_imex_diagnostics.pressure_linear_residual_linf,
                   current.pressure_linear_residual_linf);
      max_imex_diagnostics.pressure_linear_residual_relative_linf =
          max_real(max_imex_diagnostics.pressure_linear_residual_relative_linf,
                   current.pressure_linear_residual_relative_linf);
      max_imex_diagnostics.picard_pressure_solves_used =
          std::max(max_imex_diagnostics.picard_pressure_solves_used,
                   current.picard_pressure_solves_used);
      max_imex_diagnostics.kinetic_lag_update_count =
          std::max(max_imex_diagnostics.kinetic_lag_update_count,
                   current.kinetic_lag_update_count);
      max_imex_diagnostics.kinetic_lag_nonfinite_count =
          std::max(max_imex_diagnostics.kinetic_lag_nonfinite_count,
                   current.kinetic_lag_nonfinite_count);
      max_imex_diagnostics.kinetic_lag_delta_linf =
          max_real(max_imex_diagnostics.kinetic_lag_delta_linf,
                   current.kinetic_lag_delta_linf);
      max_imex_diagnostics.kinetic_picard_residual_linf =
          max_real(max_imex_diagnostics.kinetic_picard_residual_linf,
                   current.kinetic_picard_residual_linf);
      max_imex_diagnostics.pressure_picard_delta_linf =
          max_real(max_imex_diagnostics.pressure_picard_delta_linf,
                   current.pressure_picard_delta_linf);
      max_imex_diagnostics.predictor_material_speed_max =
          std::max(max_imex_diagnostics.predictor_material_speed_max,
                   current.predictor_material_speed_max);
      max_imex_diagnostics.predictor_dissipation_l1 =
          std::max(max_imex_diagnostics.predictor_dissipation_l1,
                   current.predictor_dissipation_l1);
      max_imex_diagnostics.momentum_correction_linf =
          std::max(max_imex_diagnostics.momentum_correction_linf,
                   current.momentum_correction_linf);
      max_imex_diagnostics.energy_correction_linf =
          std::max(max_imex_diagnostics.energy_correction_linf,
                   current.energy_correction_linf);
      max_imex_diagnostics.energy_flux_linf =
          std::max(max_imex_diagnostics.energy_flux_linf, current.energy_flux_linf);
      max_imex_diagnostics.energy_flux_divergence_linf =
          std::max(max_imex_diagnostics.energy_flux_divergence_linf,
                   current.energy_flux_divergence_linf);
      max_imex_diagnostics.eq78_divergence_linf =
          max_real(max_imex_diagnostics.eq78_divergence_linf,
                   current.eq78_divergence_linf);
      max_imex_diagnostics.eq78_divergence_scaled_linf =
          max_real(max_imex_diagnostics.eq78_divergence_scaled_linf,
                   current.eq78_divergence_scaled_linf);
      max_imex_diagnostics.energy_flux_closeout_consistency_linf =
          max_real(max_imex_diagnostics.energy_flux_closeout_consistency_linf,
                   current.energy_flux_closeout_consistency_linf);
      max_imex_diagnostics.face_enthalpy_denominator_guard_count_max =
          std::max(max_imex_diagnostics.face_enthalpy_denominator_guard_count_max,
                   current.face_enthalpy_denominator_guard_count);
      max_imex_diagnostics.face_enthalpy_nonfinite_count =
          std::max(max_imex_diagnostics.face_enthalpy_nonfinite_count,
                   current.face_enthalpy_nonfinite_count);
      max_imex_diagnostics.energy_flux_nonfinite_count =
          std::max(max_imex_diagnostics.energy_flux_nonfinite_count,
                   current.energy_flux_nonfinite_count);
      max_imex_diagnostics.energy_closeout_delta_l1 =
          std::max(max_imex_diagnostics.energy_closeout_delta_l1,
                   current.energy_closeout_delta_l1);
      max_imex_diagnostics.energy_closeout_delta_linf =
          std::max(max_imex_diagnostics.energy_closeout_delta_linf,
                   current.energy_closeout_delta_linf);
      max_imex_diagnostics.kinetic_pressure_contribution_l1 =
          std::max(max_imex_diagnostics.kinetic_pressure_contribution_l1,
                   current.kinetic_pressure_contribution_l1);
      max_imex_diagnostics.kinetic_pressure_contribution_linf =
          std::max(max_imex_diagnostics.kinetic_pressure_contribution_linf,
                   current.kinetic_pressure_contribution_linf);
      max_imex_diagnostics.eos_pressure_residual_l1 =
          std::max(max_imex_diagnostics.eos_pressure_residual_l1,
                   current.eos_pressure_residual_l1);
      max_imex_diagnostics.eos_pressure_residual_linf =
          std::max(max_imex_diagnostics.eos_pressure_residual_linf,
                   current.eos_pressure_residual_linf);
      max_imex_diagnostics.pressure_energy_residual_l1 =
          max_real(max_imex_diagnostics.pressure_energy_residual_l1,
                   current.pressure_energy_residual_l1);
      max_imex_diagnostics.pressure_energy_residual_linf =
          max_real(max_imex_diagnostics.pressure_energy_residual_linf,
                   current.pressure_energy_residual_linf);
      max_imex_diagnostics.pressure_energy_residual_relative_linf =
          max_real(max_imex_diagnostics.pressure_energy_residual_relative_linf,
                   current.pressure_energy_residual_relative_linf);
      max_imex_diagnostics.source_picard_iterations_used =
          std::max(max_imex_diagnostics.source_picard_iterations_used,
                   current.source_picard_iterations_used);
      max_imex_diagnostics.source_picard_delta_l2_relative_to_pmax =
          max_real(max_imex_diagnostics.source_picard_delta_l2_relative_to_pmax,
                   current.source_picard_delta_l2_relative_to_pmax);
      max_imex_diagnostics.source_picard_delta_linf =
          max_real(max_imex_diagnostics.source_picard_delta_linf,
                   current.source_picard_delta_linf);
    };

    std::string status = "ok";
    std::string failure_category = "none";
    const bool shock_density_bubble_case =
        is_shock_density_bubble_case(cfg.problem);
    std::vector<amrex::Real> shock_density_snapshot_times;
    int shock_density_next_snapshot = 0;
    int shock_density_snapshots_written = 0;
    bool shock_density_snapshot_failed = false;
    std::string shock_density_snapshot_status =
        cfg.shock_density_bubble_snapshot_dir.empty() ? "not_requested" : "pending";

    if (shock_density_bubble_case) {
      shock_density_snapshot_times =
          parse_shock_density_bubble_snapshot_times(
              cfg.shock_density_bubble_snapshot_times, cfg.stop_time);
      fill_problem_ghosts(state, geom, cfg, time);
      while (shock_density_next_snapshot <
                 static_cast<int>(shock_density_snapshot_times.size()) &&
             shock_density_snapshot_times[shock_density_next_snapshot] <=
                 time + amrex::Real(1.0e-13)) {
        if (!write_shock_density_bubble_snapshot(
                state, geom, cfg, step, time, shock_density_next_snapshot)) {
          shock_density_snapshot_failed = true;
          status = "failed";
          failure_category = "snapshot_write_failed";
          break;
        }
        ++shock_density_snapshots_written;
        ++shock_density_next_snapshot;
      }
    }

    write_plotfile(state, geom, cfg, step, time);
    const auto start_time = std::chrono::steady_clock::now();

    Diagnostics diagnostics = compute_diagnostics(state, cfg);
    if (diagnostics.nonfinite_count > 0 || diagnostics.rho_min <= amrex::Real(0.0) ||
        diagnostics.pressure_min <= amrex::Real(0.0)) {
      status = "failed";
      failure_category = shock_density_bubble_case ? "stage_nonphysical"
                                                   : "initial_nonphysical_state";
    }

    while (status == "ok" && step < cfg.max_step &&
           time < cfg.stop_time - amrex::Real(1.0e-14)) {
      amrex::Real dt = cfg.stop_time - time;
      if (cfg.method == MethodKind::Imex) {
        ImexTimeStepDiagnostics timestep_diagnostics;
        dt = compute_imex_dt(state, geom, cfg, time, timestep_diagnostics);
        imex_acoustic_startup_steps += timestep_diagnostics.acoustic_startup_used;
        imex_acoustic_cap_steps += timestep_diagnostics.acoustic_cap_used;
        imex_timestep_material_rate_max =
            std::max(imex_timestep_material_rate_max,
                     timestep_diagnostics.material_rate_max);
        imex_timestep_acoustic_rate_max =
            std::max(imex_timestep_acoustic_rate_max,
                     timestep_diagnostics.acoustic_rate_max);
      } else {
        dt = compute_euler_dt(state, geom, cfg, time);
      }

      if (!std::isfinite(dt) || dt <= amrex::Real(0.0)) {
        status = "failed";
        failure_category = shock_density_bubble_case ? "invalid_cfl_dt"
                                                     : "invalid_time_step";
        break;
      }

      if (shock_density_bubble_case &&
          shock_density_next_snapshot <
              static_cast<int>(shock_density_snapshot_times.size())) {
        const amrex::Real next_snapshot_time =
            shock_density_snapshot_times[shock_density_next_snapshot];
        if (next_snapshot_time > time + amrex::Real(1.0e-13) &&
            time + dt > next_snapshot_time + amrex::Real(1.0e-13)) {
          dt = next_snapshot_time - time;
        }
      }

      if (cfg.method == MethodKind::Imex) {
        last_imex_diagnostics =
            euler_imex_bdltv20_o1_source_map_picard_one_step(
                state, geom, cfg, dt, nullptr);
        ++imex_step_count;
        update_imex_maxima(last_imex_diagnostics);

        const bool imex_solver_ok =
            last_imex_diagnostics.solver_status == "converged_absolute" ||
            last_imex_diagnostics.solver_status == "converged_relative" ||
            last_imex_diagnostics.solver_status == "returned_finite_before_max_iter";
        if (!imex_solver_ok) {
          status = "failed";
          failure_category = shock_density_bubble_case
                                 ? "imex_pressure_solve_nonconverged"
                                 : "imex_linear_solver_" +
                                       last_imex_diagnostics.solver_status;
        }
        if (status == "ok") {
          imex_admissibility_last =
              compute_admissibility_diagnostics(state, cfg);
          update_high_trial_run(imex_admissibility_run, imex_admissibility_last);
        }
      } else {
        hllc_degenerate_star_count += euler_fv_one_step(state, geom, cfg, dt, time);
      }

      if (status != "ok") {
        break;
      }

      apply_shock_density_bubble_cylindrical_source(state, geom, cfg, dt);
      amrex::Gpu::streamSynchronize();
      dt_min_run = std::min(dt_min_run, dt);
      dt_max_run = std::max(dt_max_run, dt);
      dt_sum_run += dt;
      ++dt_count_run;
      time += dt;
      ++step;

      diagnostics = compute_diagnostics(state, cfg);
      if (diagnostics.nonfinite_count > 0) {
        status = "failed";
        failure_category = "nonfinite_state";
      } else if (diagnostics.rho_min <= amrex::Real(0.0)) {
        status = "failed";
        failure_category = "nonpositive_density";
      } else if (diagnostics.pressure_min <= amrex::Real(0.0)) {
        status = "failed";
        failure_category = "nonpositive_pressure";
      }

      if (shock_density_bubble_case && status == "ok" &&
          !shock_density_snapshot_failed) {
        fill_problem_ghosts(state, geom, cfg, time);
        while (shock_density_next_snapshot <
                   static_cast<int>(shock_density_snapshot_times.size()) &&
               shock_density_snapshot_times[shock_density_next_snapshot] <=
                   time + amrex::Real(1.0e-13)) {
          if (!write_shock_density_bubble_snapshot(
                  state, geom, cfg, step,
                  shock_density_snapshot_times[shock_density_next_snapshot],
                  shock_density_next_snapshot)) {
            shock_density_snapshot_failed = true;
            status = "failed";
            failure_category = "snapshot_write_failed";
            break;
          }
          ++shock_density_snapshots_written;
          ++shock_density_next_snapshot;
        }
      }
      write_plotfile(state, geom, cfg, step, time);
    }

    if (status == "ok" && step >= cfg.max_step &&
        time < cfg.stop_time - amrex::Real(1.0e-14)) {
      status = shock_density_bubble_case ? "failed" : "capped_max_step";
      failure_category = shock_density_bubble_case ? "max_steps_exceeded"
                                                   : "max_step_before_stop_time";
    }
    if (shock_density_bubble_case && status == "ok") {
      failure_category = "ok";
    }

    const bool final_plot_already_written =
        cfg.plot_int >= 0 && (step == 0 || step % cfg.plot_int == 0);
    if (!final_plot_already_written) {
      write_plotfile(state, geom, cfg, step, time, true);
    }

    HighTrialRunDiagnostics field_high_trial_run = imex_admissibility_run;
    if (cfg.method != MethodKind::Imex || imex_step_count == 0) {
      field_high_trial_run.rho_min = std::numeric_limits<amrex::Real>::quiet_NaN();
      field_high_trial_run.pressure_min =
          std::numeric_limits<amrex::Real>::quiet_NaN();
      field_high_trial_run.internal_energy_min =
          std::numeric_limits<amrex::Real>::quiet_NaN();
    }

    write_gresho_final_csv(state, geom, cfg, step, time, field_high_trial_run);
    write_toro_final_csv(state, geom, cfg);
    write_riemann_quadrant_final_csv(state, geom, cfg, step, time);
    write_advection_blob_final_csv(state, geom, cfg, step, time);

    const amrex::Real final_mass = state.sum(Rho) * volume;
    const amrex::Real final_energy = state.sum(E) * volume;
    const amrex::Real mass_drift = final_mass - initial_mass;
    const amrex::Real energy_drift = final_energy - initial_energy;
    const auto stop_time = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = stop_time - start_time;
    const GreshoMetrics gresho_metrics = compute_gresho_metrics(state, geom, cfg);

    ShockDensityBubbleMetrics shock_density_metrics;
    if (shock_density_bubble_case) {
      fill_problem_ghosts(state, geom, cfg, time);
      shock_density_metrics =
          compute_shock_density_bubble_metrics(state, geom, cfg);
      if (shock_density_snapshot_failed) {
        shock_density_snapshot_status = "snapshot_write_failed";
      } else if (cfg.shock_density_bubble_snapshot_dir.empty()) {
        shock_density_snapshot_status = "not_requested";
      } else if (shock_density_snapshots_written ==
                 static_cast<int>(shock_density_snapshot_times.size())) {
        shock_density_snapshot_status = "ok";
      } else {
        shock_density_snapshot_status = "missing_requested_snapshot";
        if (status == "ok") {
          status = "failed";
          failure_category = "snapshot_write_failed";
        }
      }
      const amrex::Real dt_mean_run =
          dt_count_run > 0 ? dt_sum_run / static_cast<amrex::Real>(dt_count_run)
                           : std::numeric_limits<amrex::Real>::quiet_NaN();
      const amrex::Real dt_min_output =
          dt_count_run > 0 ? dt_min_run
                           : std::numeric_limits<amrex::Real>::quiet_NaN();
      const amrex::Real dt_max_output =
          dt_count_run > 0 ? dt_max_run
                           : std::numeric_limits<amrex::Real>::quiet_NaN();
      write_shock_density_bubble_summary_csv(
          state, geom, cfg, step, time, status, failure_category, steps_rejected,
          elapsed.count(), dt_min_output, dt_max_output, dt_mean_run,
          shock_density_snapshot_status, shock_density_metrics);
    }

    amrex::Print() << "project_amrex_euler_compare_status=" << status << "\n";
    amrex::Print() << "failure_category=" << failure_category << "\n";
    amrex::Print() << "problem=" << to_string(cfg.problem) << "\n";
    amrex::Print() << "method=" << to_string(cfg.method) << "\n";
    amrex::Print() << "riemann="
                   << (cfg.method == MethodKind::Imex ? "not_applicable"
                                                       : to_string(cfg.riemann))
                   << "\n";
    amrex::Print() << "imex_form="
                   << (cfg.method == MethodKind::Imex ? to_string(cfg.imex_form)
                                                       : "not_applicable")
                   << "\n";
    amrex::Print() << "imex_predictor_flux="
                   << (cfg.method == MethodKind::Imex
                           ? imex_predictor_flux_label(cfg.imex_form)
                           : "not_applicable")
                   << "\n";
    amrex::Print() << "imex_predictor_dissipation="
                   << (cfg.method == MethodKind::Imex
                           ? to_string(cfg.imex_predictor_dissipation)
                           : "not_applicable")
                   << "\n";
    amrex::Print() << "slope_limiter=" << to_string(cfg.slope_limiter) << "\n";
    amrex::Print() << "spatial_order=" << cfg.spatial_order << "\n";
    amrex::Print() << "toro_test=" << cfg.toro_test << "\n";
    amrex::Print() << "field_boundary=" << cfg.field_boundary << "\n";
    amrex::Print() << "geometry_is_periodic=" << cfg.is_periodic[0] << " "
                   << cfg.is_periodic[1] << "\n";
    amrex::Print() << "completed_steps=" << step << "\n";
    amrex::Print() << "completed_time=" << time << "\n";
    amrex::Print() << "target_time=" << cfg.stop_time << "\n";
    amrex::Print() << "max_step=" << cfg.max_step << "\n";
    amrex::Print() << "wall_time_sec=" << elapsed.count() << "\n";
    amrex::Print() << "initial_mass=" << initial_mass << "\n";
    amrex::Print() << "final_mass=" << final_mass << "\n";
    amrex::Print() << "mass_drift=" << mass_drift << "\n";
    amrex::Print() << "initial_energy=" << initial_energy << "\n";
    amrex::Print() << "final_energy=" << final_energy << "\n";
    amrex::Print() << "energy_drift=" << energy_drift << "\n";
    amrex::Print() << "rho_min=" << diagnostics.rho_min << "\n";
    amrex::Print() << "pressure_min=" << diagnostics.pressure_min << "\n";
    amrex::Print() << "nonfinite_count=" << diagnostics.nonfinite_count << "\n";
    amrex::Print() << "hllc_degenerate_star_count="
                   << hllc_degenerate_star_count << "\n";
    amrex::Print() << "imex_step_count=" << imex_step_count << "\n";

    if (cfg.method == MethodKind::Imex) {
      amrex::Print() << "imex_cfl=" << cfg.imex_cfl << "\n";
      amrex::Print() << "imex_acoustic_cfl_cap=" << cfg.imex_acoustic_cfl_cap << "\n";
      amrex::Print() << "imex_timestep_rule="
                     << (use_bdltv20_eq79_timestep(cfg)
                             ? ((cfg.imex_acoustic_startup != 0 ||
                                 cfg.imex_acoustic_cfl_cap > amrex::Real(0.0))
                                    ? "material_speed_eq79_zdrop_2d_with_acoustic_startup_or_cap"
                                    : "material_speed_eq79_zdrop_2d")
                             : "material_speed_additive_directional_rate")
                     << "\n";
      amrex::Print() << "imex_acoustic_startup_enabled="
                     << cfg.imex_acoustic_startup << "\n";
      amrex::Print() << "imex_acoustic_startup_steps="
                     << imex_acoustic_startup_steps << "\n";
      amrex::Print() << "imex_acoustic_cap_steps=" << imex_acoustic_cap_steps
                     << "\n";
      amrex::Print() << "imex_timestep_material_rate_max="
                     << imex_timestep_material_rate_max << "\n";
      amrex::Print() << "imex_timestep_acoustic_rate_max="
                     << imex_timestep_acoustic_rate_max << "\n";
      amrex::Print() << "imex_solver=" << last_imex_diagnostics.solver << "\n";
      amrex::Print() << "imex_solver_status="
                     << last_imex_diagnostics.solver_status << "\n";
      amrex::Print() << "imex_solver_iterations_last="
                     << last_imex_diagnostics.solver_iterations << "\n";
      amrex::Print() << "imex_solver_iterations_max="
                     << max_imex_diagnostics.solver_iterations << "\n";
      amrex::Print() << "imex_solver_residual_last="
                     << last_imex_diagnostics.solver_residual << "\n";
      amrex::Print() << "imex_solver_residual_max="
                     << max_imex_diagnostics.solver_residual << "\n";
      amrex::Print() << "imex_solver_relative_residual_last="
                     << last_imex_diagnostics.solver_relative_residual << "\n";
      amrex::Print() << "imex_solver_relative_residual_max="
                     << max_imex_diagnostics.solver_relative_residual << "\n";
      amrex::Print() << "imex_pressure_solver_path="
                     << last_imex_diagnostics.pressure_solver_path << "\n";
      amrex::Print() << "imex_pressure_operator_label="
                     << last_imex_diagnostics.pressure_operator_label << "\n";
      amrex::Print() << "imex_pressure_row_stencil="
                     << last_imex_diagnostics.pressure_row_stencil << "\n";
      amrex::Print() << "imex_pressure_energy_form="
                     << last_imex_diagnostics.pressure_energy_form << "\n";
      amrex::Print() << "imex_source_primary_paper="
                     << last_imex_diagnostics.source_primary_paper << "\n";
      amrex::Print() << "imex_source_equation_set="
                     << last_imex_diagnostics.source_equation_set << "\n";
      amrex::Print() << "imex_source_time_order="
                     << last_imex_diagnostics.source_time_order << "\n";
      amrex::Print() << "imex_source_space_order="
                     << last_imex_diagnostics.source_space_order << "\n";
      amrex::Print() << "imex_source_advection_reconstruction="
                     << last_imex_diagnostics.source_advection_reconstruction
                     << "\n";
      amrex::Print() << "imex_picard_pressure_solves_used="
                     << last_imex_diagnostics.picard_pressure_solves_used << "\n";
      amrex::Print() << "imex_picard_pressure_solves_used_max="
                     << max_imex_diagnostics.picard_pressure_solves_used << "\n";
      amrex::Print() << "imex_kinetic_lag_update_count="
                     << last_imex_diagnostics.kinetic_lag_update_count << "\n";
      amrex::Print() << "imex_kinetic_lag_update_count_max="
                     << max_imex_diagnostics.kinetic_lag_update_count << "\n";
      amrex::Print() << "imex_momentum_correction_linf_last="
                     << last_imex_diagnostics.momentum_correction_linf << "\n";
      amrex::Print() << "imex_momentum_correction_linf_max="
                     << max_imex_diagnostics.momentum_correction_linf << "\n";
      amrex::Print() << "imex_energy_correction_linf_last="
                     << last_imex_diagnostics.energy_correction_linf << "\n";
      amrex::Print() << "imex_energy_correction_linf_max="
                     << max_imex_diagnostics.energy_correction_linf << "\n";
      amrex::Print() << "imex_face_enthalpy_policy="
                     << last_imex_diagnostics.face_enthalpy_policy << "\n";
      amrex::Print() << "imex_face_momentum_policy="
                     << last_imex_diagnostics.face_momentum_policy << "\n";
      amrex::Print() << "imex_eq78_divergence_linf_last="
                     << last_imex_diagnostics.eq78_divergence_linf << "\n";
      amrex::Print() << "imex_eq78_divergence_linf_max="
                     << max_imex_diagnostics.eq78_divergence_linf << "\n";
      amrex::Print() << "imex_energy_closeout_delta_linf_last="
                     << last_imex_diagnostics.energy_closeout_delta_linf << "\n";
      amrex::Print() << "imex_energy_closeout_delta_linf_max="
                     << max_imex_diagnostics.energy_closeout_delta_linf << "\n";
      amrex::Print() << "imex_solved_pressure_min="
                     << last_imex_diagnostics.solved_pressure_min << "\n";
      amrex::Print() << "imex_solved_pressure_max="
                     << last_imex_diagnostics.solved_pressure_max << "\n";
      amrex::Print() << "imex_trial_density_min="
                     << last_imex_diagnostics.trial_density_min << "\n";
      amrex::Print() << "imex_trial_pressure_min="
                     << last_imex_diagnostics.trial_pressure_min << "\n";
      amrex::Print() << "imex_trial_internal_energy_min="
                     << last_imex_diagnostics.trial_internal_energy_min << "\n";
      amrex::Print() << "imex_nonfinite_count_last="
                     << last_imex_diagnostics.nonfinite_count << "\n";
    }

    if (cfg.problem == ProblemKind::GreshoVortex) {
      amrex::Print() << "density_l1_error=" << gresho_metrics.density_l1_error
                     << "\n";
      amrex::Print() << "velocity_l1_error=" << gresho_metrics.velocity_l1_error
                     << "\n";
      amrex::Print() << "pressure_l1_error=" << gresho_metrics.pressure_l1_error
                     << "\n";
      amrex::Print() << "pressure_perturbation_l1_error="
                     << gresho_metrics.pressure_perturbation_l1_error << "\n";
      amrex::Print() << "pressure_perturbation_l1_relative_error="
                     << gresho_metrics.pressure_perturbation_l1_relative_error
                     << "\n";
      amrex::Print() << "kinetic_energy_initial="
                     << gresho_metrics.kinetic_energy_initial << "\n";
      amrex::Print() << "kinetic_energy_final="
                     << gresho_metrics.kinetic_energy_final << "\n";
      amrex::Print() << "kinetic_energy_ratio="
                     << gresho_metrics.kinetic_energy_ratio << "\n";
      amrex::Print() << "reference_sound_speed="
                     << gresho_metrics.reference_sound_speed << "\n";
      amrex::Print() << "velocity_scale=" << gresho_metrics.velocity_scale << "\n";
    }

    if (shock_density_bubble_case) {
      amrex::Print() << "source_case_id="
                     << shock_density_bubble_source_case_id(cfg) << "\n";
      amrex::Print() << "scheme_name=" << shock_density_bubble_scheme_name(cfg)
                     << "\n";
      amrex::Print() << "order=" << shock_density_bubble_order_label(cfg)
                     << "\n";
      amrex::Print() << "gamma=" << cfg.gamma << "\n";
      amrex::Print() << "grid_nx=" << cfg.n_cell[0] << "\n";
      amrex::Print() << "grid_ny=" << cfg.n_cell[1] << "\n";
      amrex::Print() << "shock_initial_x=" << shock_density_bubble_shock_x()
                     << "\n";
      amrex::Print() << "post_shock_density="
                     << shock_density_bubble_post_shock_state().rho << "\n";
      amrex::Print() << "post_shock_velocity="
                     << shock_density_bubble_post_shock_state().u << "\n";
      amrex::Print() << "post_shock_pressure="
                     << shock_density_bubble_post_shock_state().p << "\n";
      amrex::Print() << "bubble_density=0.1\n";
      amrex::Print() << "bubble_center_x=" << shock_density_bubble_center_x()
                     << "\n";
      amrex::Print() << "bubble_center_y=" << shock_density_bubble_center_y()
                     << "\n";
      amrex::Print() << "bubble_radius=" << shock_density_bubble_radius() << "\n";
      amrex::Print() << "boundary_left=fixed_post_shock_inflow\n";
      amrex::Print() << "boundary_right=outflow_zero_gradient\n";
      amrex::Print() << "boundary_bottom=reflective_symmetry\n";
      amrex::Print() << "boundary_top=outflow_zero_gradient\n";
      amrex::Print() << "snapshot_times="
                     << cfg.shock_density_bubble_snapshot_times << "\n";
      amrex::Print() << "snapshot_status=" << shock_density_snapshot_status
                     << "\n";
      amrex::Print() << "snapshot_dir="
                     << cfg.shock_density_bubble_snapshot_dir << "\n";
      amrex::Print() << "snapshots_written=" << shock_density_snapshots_written
                     << "\n";
      amrex::Print() << "density_bubble_centroid_x="
                     << shock_density_metrics.density_bubble_centroid_x << "\n";
      amrex::Print() << "density_bubble_centroid_y="
                     << shock_density_metrics.density_bubble_centroid_y << "\n";
      amrex::Print() << "density_bubble_area_threshold="
                     << shock_density_metrics.density_bubble_area_threshold
                     << "\n";
      amrex::Print() << "schlieren_max=" << shock_density_metrics.schlieren_max
                     << "\n";
      amrex::Print() << "internal_energy_min="
                     << shock_density_metrics.internal_energy_min << "\n";
      amrex::Print() << "validation_claim="
                     << shock_density_bubble_validation_claim(cfg) << "\n";
      amrex::Print() << "gfm_used=false\n";
      amrex::Print() << "level_set_used=false\n";
      amrex::Print() << "material_count=1\n";
      amrex::Print() << "cylindrical_source_used="
                     << (shock_density_bubble_uses_cylindrical_source(cfg)
                             ? "true"
                             : "false")
                     << "\n";
      amrex::Print() << "geometric_source_form="
                     << (shock_density_bubble_uses_cylindrical_source(cfg)
                             ? "radial_symmetry_fractional_source_step"
                             : "none")
                     << "\n";
    }

    if (status != "ok" && status != "capped_max_step") {
      return_code = 2;
    }
  }

finalize_amrex:
  amrex::Finalize();
  return return_code;
}

}  // namespace euler_compare
