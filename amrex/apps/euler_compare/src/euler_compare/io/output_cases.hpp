// CSV, plotfile, snapshot, and error-output routines.
std::string gather_text_to_io_processor(const std::string& local_text)
{
  const int root = amrex::ParallelDescriptor::IOProcessorNumber();
  if (amrex::ParallelDescriptor::NProcs() == 1) {
    return local_text;
  }

  if (local_text.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    amrex::Abort("CSV output from one rank is too large for MPI Gatherv.");
  }

  const int local_count = static_cast<int>(local_text.size());
  const std::vector<int> counts = amrex::ParallelDescriptor::Gather(local_count, root);
  std::vector<int> offsets;
  std::vector<char> gathered;
  if (amrex::ParallelDescriptor::IOProcessor()) {
    offsets.resize(counts.size(), 0);
    int total_count = 0;
    for (std::size_t n = 0; n < counts.size(); ++n) {
      offsets[n] = total_count;
      total_count += counts[n];
    }
    gathered.resize(static_cast<std::size_t>(total_count));
  }

  amrex::ParallelDescriptor::Gatherv(
      local_text.data(), local_count, gathered.data(), counts, offsets, root);
  if (!amrex::ParallelDescriptor::IOProcessor()) {
    return std::string();
  }
  return std::string(gathered.begin(), gathered.end());
}

bool write_gathered_csv_file(const std::string& filename,
                             const std::string& header,
                             const std::string& local_body)
{
  const std::string gathered_body = gather_text_to_io_processor(local_body);
  if (!amrex::ParallelDescriptor::IOProcessor()) {
    return true;
  }
  std::ofstream out(filename);
  if (!out) {
    return false;
  }
  out << header << gathered_body;
  return static_cast<bool>(out);
}

void write_plotfile(const amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg, int step,
                    amrex::Real time, bool force = false)
{
  if (cfg.plot_int < 0) {
    return;
  }
  if (!force && step != 0 && step % cfg.plot_int != 0) {
    return;
  }
  amrex::MultiFab plot = make_plot_state(state, cfg);
  const amrex::Vector<std::string> names{"rho", "momx", "momy", "E", "pressure"};
  amrex::WriteSingleLevelPlotfile(cfg.plot_file + amrex::Concatenate("", step, 5), plot, names, geom, time, step);
}

void write_gresho_final_csv(const amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg,
                            int step, amrex::Real time,
                            const HighTrialRunDiagnostics& high_trial_run)
{
  if (cfg.problem != ProblemKind::GreshoVortex || cfg.final_csv.empty()) {
    return;
  }

  const std::string header =
      "x,y,rho,u,v,pressure,exact_rho,exact_u,exact_v,exact_pressure,"
      "density_error,velocity_error,pressure_error,pressure_perturbation_error,"
      "time,step,mach,method,riemann,imex_form,imex_predictor_flux,imex_predictor_dissipation,"
      "spatial_order,slope_limiter,"
      "imex_trial_density_run_min,"
      "imex_high_trial_pressure_run_min,imex_high_trial_internal_energy_run_min\n";
  std::ostringstream rows;
  rows << std::setprecision(17);

  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real p0 = cfg.pressure;
  const char* riemann_label = cfg.method == MethodKind::Imex ? "not_applicable" : to_string(cfg.riemann);
  const char* imex_form_label = cfg.method == MethodKind::Imex ? to_string(cfg.imex_form) : "not_applicable";
  const char* imex_predictor_flux =
      cfg.method == MethodKind::Imex ? imex_predictor_flux_label(cfg.imex_form) : "not_applicable";
  const char* imex_predictor_dissipation =
      cfg.method == MethodKind::Imex ? to_string(cfg.imex_predictor_dissipation) : "not_applicable";
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        const amrex::Real x = plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
        const amrex::Real y = plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
        const PrimitiveState numerical = to_primitive(load_cons(arr, i, j, 0), cfg.gamma);
        const PrimitiveState exact = gresho_state(x, y, cfg);
        const amrex::Real du = numerical.u - exact.u;
        const amrex::Real dv = numerical.v - exact.v;
        rows << x << ',' << y << ',' << numerical.rho << ',' << numerical.u << ',' << numerical.v << ','
             << numerical.p << ',' << exact.rho << ',' << exact.u << ',' << exact.v << ',' << exact.p << ','
             << numerical.rho - exact.rho << ',' << std::sqrt(du * du + dv * dv) << ','
             << numerical.p - exact.p << ',' << (numerical.p - p0) - (exact.p - p0) << ',' << time << ','
             << step << ',' << cfg.mach << ',' << to_string(cfg.method) << ',' << riemann_label << ','
             << imex_form_label << ',' << imex_predictor_flux << ',' << imex_predictor_dissipation << ','
             << cfg.spatial_order << ',' << to_string(cfg.slope_limiter) << ','
             << high_trial_run.rho_min << ',' << high_trial_run.pressure_min << ','
             << high_trial_run.internal_energy_min << '\n';
      }
    }
  }
  write_gathered_csv_file(cfg.final_csv, header, rows.str());
}

void write_toro_final_csv(const amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg)
{
  if (cfg.problem != ProblemKind::Toro1 || cfg.final_csv.empty()) {
    return;
  }

  const std::string header =
      "i,j,x,y,rho,u,v,p,exact_rho,exact_u,exact_v,exact_p,"
      "rho_error,u_error,v_error,p_error\n";
  std::ostringstream rows;
  rows << std::setprecision(17);

  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
      for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
        const amrex::Real x =
            plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
        const amrex::Real y =
            plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
        const ConservedState c = load_cons(arr, i, j, 0);
        const PrimitiveState prim = to_primitive(c, cfg.gamma);
        rows << i << ',' << j << ',' << x << ',' << y << ',' << prim.rho << ','
             << prim.u << ',' << prim.v << ',' << prim.p << ',' << nan << ','
             << nan << ',' << nan << ',' << nan << ',' << nan << ',' << nan << ','
             << nan << ',' << nan << '\n';
      }
    }
  }
  write_gathered_csv_file(cfg.final_csv, header, rows.str());
}

void write_riemann_quadrant_final_csv(const amrex::MultiFab& state,
                                      const amrex::Geometry& geom,
                                      const RunConfig& cfg,
                                      int step,
                                      amrex::Real time)
{
  if (cfg.problem != ProblemKind::RiemannQuadrant || cfg.final_csv.empty()) {
    return;
  }

  const std::string header =
      "i,j,x,y,rho,u,v,p,initial_rho,initial_u,initial_v,initial_p,"
      "time,step,method,riemann,spatial_order,slope_limiter,"
      "field_boundary,geometry_is_periodic_x,geometry_is_periodic_y,"
      "quadrant_case,claim_limit\n";
  std::ostringstream rows;
  rows << std::setprecision(17);

  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const char* riemann_label =
      cfg.method == MethodKind::Imex ? "imex_pressure_split" : to_string(cfg.riemann);
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
      for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
        const amrex::Real x =
            plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
        const amrex::Real y =
            plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
        const PrimitiveState numerical =
            to_primitive(load_cons(arr, i, j, 0), cfg.gamma);
        const PrimitiveState initial = riemann_quadrant_state(x, y, cfg);
        rows << i << ',' << j << ',' << x << ',' << y << ','
             << numerical.rho << ',' << numerical.u << ',' << numerical.v << ','
             << numerical.p << ',' << initial.rho << ',' << initial.u << ','
             << initial.v << ',' << initial.p << ',' << time << ',' << step << ','
             << to_string(cfg.method) << ',' << riemann_label << ','
             << cfg.spatial_order << ',' << to_string(cfg.slope_limiter) << ','
             << cfg.field_boundary << ',' << geom.isPeriodic(0) << ','
             << geom.isPeriodic(1) << ','
             << "liska_wendroff_clawpack_quadrant,"
             << "genuine_2d_visual_finite_smoke_no_exact_solution\n";
      }
    }
  }
  write_gathered_csv_file(cfg.final_csv, header, rows.str());
}

void write_advection_blob_final_csv(const amrex::MultiFab& state,
                                    const amrex::Geometry& geom,
                                    const RunConfig& cfg,
                                    int step,
                                    amrex::Real time)
{
  if (cfg.problem != ProblemKind::AdvectionBlob || cfg.final_csv.empty()) {
    return;
  }

  const std::string header =
      "x,y,rho,u,v,pressure,exact_rho,exact_u,exact_v,exact_pressure,"
      "density_error,velocity_error,pressure_error,time,step,method,riemann,"
      "spatial_order,slope_limiter,geometry_is_periodic_x,geometry_is_periodic_y\n";
  std::ostringstream rows;
  rows << std::setprecision(17);

  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const char* riemann_label =
      cfg.method == MethodKind::Imex ? "imex_pressure_split" : to_string(cfg.riemann);
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        const amrex::Real x =
            plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
        const amrex::Real y =
            plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
        const PrimitiveState numerical = to_primitive(load_cons(arr, i, j, 0), cfg.gamma);
        const PrimitiveState exact = exact_advection_blob_state(x, y, geom, cfg, time);
        const amrex::Real du = numerical.u - exact.u;
        const amrex::Real dv = numerical.v - exact.v;
        rows << x << ',' << y << ',' << numerical.rho << ',' << numerical.u << ','
             << numerical.v << ',' << numerical.p << ',' << exact.rho << ','
             << exact.u << ',' << exact.v << ',' << exact.p << ','
             << numerical.rho - exact.rho << ',' << std::sqrt(du * du + dv * dv)
             << ',' << numerical.p - exact.p << ',' << time << ',' << step << ','
             << to_string(cfg.method) << ',' << riemann_label << ','
             << cfg.spatial_order << ',' << to_string(cfg.slope_limiter) << ','
             << geom.isPeriodic(0) << ',' << geom.isPeriodic(1) << '\n';
      }
    }
  }
  write_gathered_csv_file(cfg.final_csv, header, rows.str());
}

struct LocalCellValues {
  int i = 0;
  int j = 0;
  amrex::Real x = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real y = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real rho = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real u = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real v = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real p = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real internal_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real kinetic_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real mx = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real my = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real total_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
};

struct FailureExtrema {
  int density_i = -1;
  int density_j = -1;
  int pressure_i = -1;
  int pressure_j = -1;
  int internal_i = -1;
  int internal_j = -1;
  int pressure_grad_i = -1;
  int pressure_grad_j = -1;
  int internal_grad_i = -1;
  int internal_grad_j = -1;
  amrex::Real density_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real pressure_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real internal_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real pressure_gradient_max = amrex::Real(0.0);
  amrex::Real internal_gradient_max = amrex::Real(0.0);
};

std::string csv_escape(const std::string& value)
{
  if (value.find_first_of(",\"\n") == std::string::npos) {
    return value;
  }
  std::string escaped = "\"";
  for (const char c : value) {
    if (c == '"') {
      escaped += "\"\"";
    } else {
      escaped += c;
    }
  }
  escaped += '"';
  return escaped;
}

struct ShockDensityBubbleMetrics {
  amrex::Real rho_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real rho_max = -std::numeric_limits<amrex::Real>::infinity();
  amrex::Real pressure_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real pressure_max = -std::numeric_limits<amrex::Real>::infinity();
  amrex::Real internal_energy_min = std::numeric_limits<amrex::Real>::infinity();
  int nonfinite_count = 0;
  amrex::Real density_bubble_centroid_x = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real density_bubble_centroid_y = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real density_bubble_area_threshold = amrex::Real(0.55);
  amrex::Real density_bubble_area = amrex::Real(0.0);
  amrex::Real schlieren_max = amrex::Real(0.0);
};

std::vector<amrex::Real>
parse_shock_density_bubble_snapshot_times(const std::string& text,
                                          amrex::Real final_time)
{
  std::vector<amrex::Real> times;
  std::stringstream input(text);
  std::string token;
  while (std::getline(input, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(),
                               [](unsigned char c) { return std::isspace(c); }),
                token.end());
    if (token.empty()) {
      continue;
    }
    const amrex::Real t = static_cast<amrex::Real>(std::stod(token));
    if (!std::isfinite(t) || t < amrex::Real(-1.0e-14) ||
        t > final_time + amrex::Real(1.0e-14)) {
      amrex::Abort("shock_density_bubble_2d snapshot_times must lie in [0, final_time].");
    }
    times.push_back(std::max(amrex::Real(0.0), std::min(t, final_time)));
  }
  std::sort(times.begin(), times.end());
  times.erase(std::unique(times.begin(), times.end(),
                          [](amrex::Real a, amrex::Real b) {
                            return std::abs(a - b) < amrex::Real(1.0e-13);
                          }),
              times.end());
  return times;
}

std::string shock_density_bubble_scheme_name(const RunConfig& cfg)
{
  if (cfg.method == MethodKind::Imex) {
    return "imex_t1s2_bdltv20_ideal_gas";
  }
  if (cfg.riemann == RiemannKind::XieAmHllcP) {
    return "explicit_o2_low_mach_corrected_hllc_p";
  }
  if (cfg.riemann == RiemannKind::Hllc) {
    return "explicit_o2_hllc";
  }
  return std::string("explicit_o") + std::to_string(cfg.spatial_order) + "_" +
         to_string(cfg.riemann);
}

std::string shock_density_bubble_order_label(const RunConfig& cfg)
{
  if (cfg.method == MethodKind::Imex) {
    return "T1/S2";
  }
  return std::string("O") + std::to_string(cfg.spatial_order);
}

ShockDensityBubbleMetrics
compute_shock_density_bubble_metrics(const amrex::MultiFab& state,
                                     const amrex::Geometry& geom,
                                     const RunConfig& cfg)
{
  ShockDensityBubbleMetrics metrics;
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real cell_vol = cell_volume(geom);
  const amrex::Real threshold = metrics.density_bubble_area_threshold;
  amrex::Real centroid_x_sum = amrex::Real(0.0);
  amrex::Real centroid_y_sum = amrex::Real(0.0);

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
      for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
        const ConservedState u = load_cons(arr, i, j, 0);
        const PrimitiveState q = to_primitive(u, gamma);
        const amrex::Real kinetic =
            u.rho > amrex::Real(0.0)
                ? amrex::Real(0.5) * (u.mx * u.mx + u.my * u.my) / u.rho
                : std::numeric_limits<amrex::Real>::quiet_NaN();
        const amrex::Real internal = u.E - kinetic;
        const bool bad = !std::isfinite(u.rho) || !std::isfinite(u.mx) ||
                         !std::isfinite(u.my) || !std::isfinite(u.E) ||
                         !std::isfinite(q.p) || !std::isfinite(internal);
        if (bad) {
          ++metrics.nonfinite_count;
          continue;
        }
        metrics.rho_min = std::min(metrics.rho_min, q.rho);
        metrics.rho_max = std::max(metrics.rho_max, q.rho);
        metrics.pressure_min = std::min(metrics.pressure_min, q.p);
        metrics.pressure_max = std::max(metrics.pressure_max, q.p);
        metrics.internal_energy_min = std::min(metrics.internal_energy_min, internal);

        const amrex::Real rho_x =
            (arr(i + 1, j, 0, Rho) - arr(i - 1, j, 0, Rho)) /
            (amrex::Real(2.0) * dx[0]);
        const amrex::Real rho_y =
            (arr(i, j + 1, 0, Rho) - arr(i, j - 1, 0, Rho)) /
            (amrex::Real(2.0) * dx[1]);
        const amrex::Real grad = std::sqrt(rho_x * rho_x + rho_y * rho_y);
        if (std::isfinite(grad)) {
          metrics.schlieren_max = std::max(metrics.schlieren_max, grad);
        }

        if (q.rho < threshold) {
          const amrex::Real x =
              plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
          const amrex::Real y =
              plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
          metrics.density_bubble_area += cell_vol;
          centroid_x_sum += x * cell_vol;
          centroid_y_sum += y * cell_vol;
        }
      }
    }
  }

  amrex::ParallelDescriptor::ReduceRealMin(metrics.rho_min);
  amrex::ParallelDescriptor::ReduceRealMax(metrics.rho_max);
  amrex::ParallelDescriptor::ReduceRealMin(metrics.pressure_min);
  amrex::ParallelDescriptor::ReduceRealMax(metrics.pressure_max);
  amrex::ParallelDescriptor::ReduceRealMin(metrics.internal_energy_min);
  amrex::ParallelDescriptor::ReduceIntSum(metrics.nonfinite_count);
  amrex::ParallelDescriptor::ReduceRealSum(metrics.density_bubble_area);
  amrex::ParallelDescriptor::ReduceRealSum(centroid_x_sum);
  amrex::ParallelDescriptor::ReduceRealSum(centroid_y_sum);
  amrex::ParallelDescriptor::ReduceRealMax(metrics.schlieren_max);
  if (metrics.density_bubble_area > amrex::Real(0.0)) {
    metrics.density_bubble_centroid_x = centroid_x_sum / metrics.density_bubble_area;
    metrics.density_bubble_centroid_y = centroid_y_sum / metrics.density_bubble_area;
  }
  return metrics;
}

std::string shock_density_bubble_snapshot_filename(const std::string& dir,
                                                   int index,
                                                   amrex::Real time,
                                                   int step)
{
  std::ostringstream tstream;
  tstream << std::fixed << std::setprecision(6) << time;
  std::string time_label = tstream.str();
  std::replace(time_label.begin(), time_label.end(), '.', 'p');

  std::ostringstream name;
  name << dir << "/snapshot_" << std::setw(2) << std::setfill('0') << index
       << "_t" << time_label << "_step" << std::setw(6) << std::setfill('0')
       << step << ".csv";
  return name.str();
}

bool write_shock_density_bubble_snapshot(const amrex::MultiFab& state,
                                         const amrex::Geometry& geom,
                                         const RunConfig& cfg,
                                         int step,
                                         amrex::Real time,
                                         int snapshot_index)
{
  if (!is_shock_density_bubble_case(cfg.problem) ||
      cfg.shock_density_bubble_snapshot_dir.empty()) {
    return true;
  }

  std::error_code ec;
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::filesystem::create_directories(cfg.shock_density_bubble_snapshot_dir, ec);
  }

  const std::string filename =
      shock_density_bubble_snapshot_filename(cfg.shock_density_bubble_snapshot_dir,
                                             snapshot_index, time, step);

  const std::string header =
      "source_case_id,step,time,i,j,x,y,rho,u,v,pressure,internal_energy,"
      "schlieren,density_bubble_indicator,density_bubble_area_threshold,gamma,"
      "method,scheme_name,riemann,spatial_order,validation_claim,gfm_used,"
      "level_set_used,material_count,cylindrical_source_used,"
      "geometric_source_form\n";
  std::ostringstream rows;
  rows << std::setprecision(17);

  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real threshold = amrex::Real(0.55);
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
      for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
        const amrex::Real x =
            plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
        const amrex::Real y =
            plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
        const ConservedState u = load_cons(arr, i, j, 0);
        const PrimitiveState q = to_primitive(u, gamma);
        const amrex::Real kinetic =
            amrex::Real(0.5) * (u.mx * u.mx + u.my * u.my) / u.rho;
        const amrex::Real internal = u.E - kinetic;
        const amrex::Real rho_x =
            (arr(i + 1, j, 0, Rho) - arr(i - 1, j, 0, Rho)) /
            (amrex::Real(2.0) * dx[0]);
        const amrex::Real rho_y =
            (arr(i, j + 1, 0, Rho) - arr(i, j - 1, 0, Rho)) /
            (amrex::Real(2.0) * dx[1]);
        const amrex::Real schlieren = std::sqrt(rho_x * rho_x + rho_y * rho_y);
        rows << shock_density_bubble_source_case_id(cfg) << ','
             << step << ',' << time << ',' << i << ',' << j << ',' << x << ','
             << y << ',' << q.rho << ',' << q.u << ',' << q.v << ',' << q.p
             << ',' << internal << ',' << schlieren << ','
             << (q.rho < threshold ? 1 : 0) << ',' << threshold << ','
             << gamma << ',' << to_string(cfg.method) << ','
             << shock_density_bubble_scheme_name(cfg) << ','
             << (cfg.method == MethodKind::Imex ? "imex_pressure_split" : to_string(cfg.riemann))
             << ',' << cfg.spatial_order << ','
             << shock_density_bubble_validation_claim(cfg)
             << ",false,false,1,"
             << (shock_density_bubble_uses_cylindrical_source(cfg) ? "true" : "false")
             << ','
             << (shock_density_bubble_uses_cylindrical_source(cfg)
                     ? "clawpack_radial_symmetry_fractional_source_step"
                     : "none")
             << "\n";
      }
    }
  }
  return write_gathered_csv_file(filename, header, rows.str());
}

void write_shock_density_bubble_summary_csv(
    const amrex::MultiFab& state,
    const amrex::Geometry& geom,
    const RunConfig& cfg,
    int step,
    amrex::Real time,
    const std::string& status,
    const std::string& failure_category,
    int steps_rejected,
    amrex::Real wall_time_sec,
    amrex::Real dt_min,
    amrex::Real dt_max,
    amrex::Real dt_mean,
    const std::string& snapshot_status,
    const ShockDensityBubbleMetrics& metrics)
{
  if (!is_shock_density_bubble_case(cfg.problem) || cfg.final_csv.empty()) {
    return;
  }
  if (!amrex::ParallelDescriptor::IOProcessor()) {
    return;
  }

  std::ofstream out(cfg.final_csv);
  if (!out) {
    return;
  }
  out << std::setprecision(17);
  out << "source_case_id,run_id,backend,scheme_family,scheme_name,order,riemann,"
         "gamma,grid_nx,grid_ny,x_lower,x_upper,y_lower,y_upper,shock_initial_x,"
         "post_shock_density,post_shock_velocity,post_shock_pressure,bubble_density,"
         "bubble_center_x,bubble_center_y,bubble_radius,boundary_left,boundary_right,"
         "boundary_bottom,boundary_top,final_time,snapshot_times,status,"
         "failure_category,steps_accepted,steps_rejected,wall_time_sec,dt_min,dt_max,"
         "dt_mean,rho_min,rho_max,pressure_min,pressure_max,internal_energy_min,"
         "nonfinite_count,density_bubble_centroid_x,density_bubble_centroid_y,"
         "density_bubble_area_threshold,schlieren_max,snapshot_status,snapshot_dir,"
         "validation_claim,gfm_used,level_set_used,material_count,"
         "cylindrical_source_used,geometric_source_form\n";

  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const PrimitiveState post = shock_density_bubble_post_shock_state();
  const std::string run_id = shock_density_bubble_scheme_name(cfg) + "_" +
                             std::to_string(cfg.n_cell[0]) + "x" +
                             std::to_string(cfg.n_cell[1]);
  out << shock_density_bubble_source_case_id(cfg) << ','
      << csv_escape(run_id) << ",amrex_single_level,"
      << to_string(cfg.method) << ','
      << csv_escape(shock_density_bubble_scheme_name(cfg)) << ','
      << csv_escape(shock_density_bubble_order_label(cfg)) << ','
      << (cfg.method == MethodKind::Imex ? "imex_pressure_split" : to_string(cfg.riemann))
      << ',' << cfg.gamma << ',' << cfg.n_cell[0] << ',' << cfg.n_cell[1]
      << ',' << plo[0] << ',' << phi[0] << ',' << plo[1] << ',' << phi[1]
      << ',' << shock_density_bubble_shock_x() << ',' << post.rho << ','
      << post.u << ',' << post.p << ',' << amrex::Real(0.1) << ','
      << shock_density_bubble_center_x() << ',' << shock_density_bubble_center_y()
      << ',' << shock_density_bubble_radius()
      << ",fixed_post_shock_inflow,outflow_zero_gradient,reflective_symmetry,"
      << "outflow_zero_gradient," << cfg.stop_time << ','
      << csv_escape(cfg.shock_density_bubble_snapshot_times) << ','
      << status << ',' << failure_category << ',' << step << ','
      << steps_rejected << ',' << wall_time_sec << ',' << dt_min << ','
      << dt_max << ',' << dt_mean << ',' << metrics.rho_min << ','
      << metrics.rho_max << ',' << metrics.pressure_min << ','
      << metrics.pressure_max << ',' << metrics.internal_energy_min << ','
      << metrics.nonfinite_count << ',' << metrics.density_bubble_centroid_x
      << ',' << metrics.density_bubble_centroid_y << ','
      << metrics.density_bubble_area_threshold << ',' << metrics.schlieren_max
      << ',' << snapshot_status << ','
      << csv_escape(cfg.shock_density_bubble_snapshot_dir)
      << ',' << shock_density_bubble_validation_claim(cfg)
      << ",false,false,1,"
      << (shock_density_bubble_uses_cylindrical_source(cfg) ? "true" : "false")
      << ','
      << (shock_density_bubble_uses_cylindrical_source(cfg)
              ? "clawpack_radial_symmetry_fractional_source_step"
              : "none")
      << "\n";
  (void)state;
}
