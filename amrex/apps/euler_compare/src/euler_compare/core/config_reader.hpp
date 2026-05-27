// ParmParse reader and geometry construction for AMReX runtime options.
#include "config_parse_helpers.hpp"

RunConfig read_config()
{
  RunConfig cfg;

  amrex::ParmParse pp;
  const bool max_step_set = pp.query("max_step", cfg.max_step);
  pp.query("stop_time", cfg.stop_time);

  amrex::ParmParse pp_amr("amr");
  amrex::Vector<int> n_cell(AMREX_SPACEDIM);
  const bool amr_n_cell_set =
      pp_amr.countval("n_cell") >= static_cast<int>(AMREX_SPACEDIM);
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    n_cell[d] = cfg.n_cell[d];
  }
  pp_amr.queryarr("n_cell", n_cell, 0, AMREX_SPACEDIM);
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    cfg.n_cell[d] = n_cell[d];
  }
  int nx_alias = cfg.n_cell[0];
  int ny_alias = cfg.n_cell[1];
  const bool nx_alias_set = pp.query("nx", nx_alias);
  if (nx_alias_set) {
    cfg.n_cell[0] = nx_alias;
  }
  const bool ny_alias_set = pp.query("ny", ny_alias);
  if (ny_alias_set) {
    cfg.n_cell[1] = ny_alias;
  }
  const bool grid_size_set = amr_n_cell_set || nx_alias_set || ny_alias_set;
  pp_amr.query("max_grid_size", cfg.max_grid_size);
  pp_amr.query("plot_file", cfg.plot_file);
  pp_amr.query("plot_int", cfg.plot_int);

  amrex::ParmParse pp_geom("geometry");
  amrex::Vector<amrex::Real> prob_lo(AMREX_SPACEDIM);
  amrex::Vector<amrex::Real> prob_hi(AMREX_SPACEDIM);
  amrex::Vector<int> is_periodic(AMREX_SPACEDIM);
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    prob_lo[d] = cfg.prob_lo[d];
    prob_hi[d] = cfg.prob_hi[d];
    is_periodic[d] = cfg.is_periodic[d];
  }
  pp_geom.queryarr("prob_lo", prob_lo, 0, AMREX_SPACEDIM);
  pp_geom.queryarr("prob_hi", prob_hi, 0, AMREX_SPACEDIM);
  pp_geom.queryarr("is_periodic", is_periodic, 0, AMREX_SPACEDIM);
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    cfg.prob_lo[d] = prob_lo[d];
    cfg.prob_hi[d] = prob_hi[d];
    cfg.is_periodic[d] = is_periodic[d];
  }
  pp.query("x_lower", cfg.prob_lo[0]);
  pp.query("x_upper", cfg.prob_hi[0]);
  pp.query("y_lower", cfg.prob_lo[1]);
  pp.query("y_upper", cfg.prob_hi[1]);

  amrex::ParmParse pp_euler("euler");
  pp_euler.query("gamma", cfg.gamma);
  pp_euler.query("cfl", cfg.cfl);
  pp_euler.query("imex_cfl", cfg.imex_cfl);
  pp_euler.query("imex_acoustic_cfl_cap", cfg.imex_acoustic_cfl_cap);
  pp_euler.query("imex_solver_tol", cfg.imex_solver_tol);
  pp_euler.query("imex_solver_max_iter", cfg.imex_solver_max_iter);
  pp_euler.query("imex_solver_verbose", cfg.imex_solver_verbose);
  pp_euler.query("imex_acoustic_startup", cfg.imex_acoustic_startup);
  pp_euler.query("imex_picard_extra_iterations", cfg.imex_picard_extra_iterations);
  pp_euler.query("imex_picard_iterations", cfg.imex_picard_iterations);

  std::string problem = to_string(cfg.problem);
  std::string riemann = to_string(cfg.riemann);
  std::string method = to_string(cfg.method);
  std::string imex_form = to_string(cfg.imex_form);
  std::string slope_limiter = to_string(cfg.slope_limiter);
  std::string imex_predictor_dissipation = to_string(cfg.imex_predictor_dissipation);
  pp_euler.query("problem", problem);
  pp_euler.query("riemann", riemann);
  pp_euler.query("method", method);
  pp_euler.query("imex_form", imex_form);
  pp_euler.query("slope_limiter", slope_limiter);
  pp_euler.query("imex_predictor_dissipation", imex_predictor_dissipation);
  pp_euler.query("imex_pressure_solver_path", cfg.imex_pressure_solver_path);
  pp_euler.query("field_boundary", cfg.field_boundary);
  pp_euler.query("imex_pressure_stabilization", cfg.imex_pressure_stabilization);
  pp_euler.query("bdltv20_paper_t1_s2", cfg.bdltv20_paper_t1_s2);
  pp_euler.query("bdltv20_paper_pressure_solver",
                 cfg.bdltv20_paper_pressure_solver);
  pp_euler.query("bdltv20_paper_t1_s2_dt", cfg.bdltv20_paper_t1_s2_dt);
  pp_euler.query("bdltv20_paper_epsilon", cfg.bdltv20_paper_epsilon);
  pp_euler.query("bdltv20_paper_t1_s2_max_steps",
                 cfg.bdltv20_paper_t1_s2_max_steps);
  pp_euler.query("shock_density_bubble_snapshot_dir",
                 cfg.shock_density_bubble_snapshot_dir);
  pp_euler.query("shock_density_bubble_snapshot_times",
                 cfg.shock_density_bubble_snapshot_times);
  pp_euler.query("spatial_order", cfg.spatial_order);

  pp.query("case", problem);
  pp.query("problem", problem);
  pp.query("method", method);
  pp.query("riemann", riemann);
  pp.query("order", cfg.spatial_order);
  pp.query("spatial_order", cfg.spatial_order);
  pp.query("gamma", cfg.gamma);
  pp.query("cfl", cfg.cfl);
  pp.query("final_time", cfg.stop_time);
  pp.query("snapshot_dir", cfg.shock_density_bubble_snapshot_dir);
  pp.query("snapshot_times", cfg.shock_density_bubble_snapshot_times);

  cfg.problem = parse_problem(problem);
  cfg.riemann = parse_riemann(riemann);
  cfg.method = parse_method(method);
  cfg.imex_form = parse_imex_form(imex_form);
  cfg.slope_limiter = parse_slope_limiter(slope_limiter);
  cfg.imex_predictor_dissipation = parse_imex_predictor_dissipation(imex_predictor_dissipation);
  cfg.bdltv20_paper_t1_s2 = lower_string(cfg.bdltv20_paper_t1_s2);
  cfg.bdltv20_paper_pressure_solver =
      lower_string(cfg.bdltv20_paper_pressure_solver);
  cfg.field_boundary = lower_string(cfg.field_boundary);

  if (cfg.field_boundary != "outflow" && cfg.field_boundary != "exact_dirichlet") {
    amrex::Abort("euler.field_boundary must be outflow or exact_dirichlet.");
  }
  if (cfg.spatial_order != 1 && cfg.spatial_order != 2) {
    amrex::Abort("euler.spatial_order must be 1 or 2.");
  }
  if (cfg.imex_picard_extra_iterations < 0) {
    amrex::Abort("euler.imex_picard_extra_iterations must be non-negative.");
  }
  if (cfg.imex_picard_iterations < 1) {
    amrex::Abort("euler.imex_picard_iterations must be at least one.");
  }

  if (cfg.bdltv20_paper_t1_s2 != "off") {
    if (cfg.bdltv20_paper_t1_s2 != "toro_x_exact_dirichlet_y_periodic_2d" &&
        cfg.bdltv20_paper_t1_s2 != "toro_xy_exact_dirichlet_2d" &&
        cfg.bdltv20_paper_t1_s2 != "isentropic_vortex_periodic_2d" &&
        cfg.bdltv20_paper_t1_s2 != "gresho_exact_dirichlet_2d" &&
        cfg.bdltv20_paper_t1_s2 != "advection_blob_periodic_2d") {
      amrex::Abort("Unsupported euler.bdltv20_paper_t1_s2 setting.");
    }
    if (cfg.method != MethodKind::Imex) {
      amrex::Abort("euler.bdltv20_paper_t1_s2 requires euler.method=imex.");
    }
    if (cfg.spatial_order != 2) {
      amrex::Abort("euler.bdltv20_paper_t1_s2 requires euler.spatial_order=2.");
    }
    if (cfg.slope_limiter != SlopeLimiterKind::Minmod) {
      amrex::Abort("euler.bdltv20_paper_t1_s2 requires euler.slope_limiter=minmod.");
    }
    if (!(cfg.gamma > amrex::Real(1.0)) || !std::isfinite(cfg.gamma)) {
      amrex::Abort("euler.bdltv20_paper_t1_s2 requires finite gamma > 1.");
    }
    if (!(cfg.imex_solver_tol > amrex::Real(0.0)) ||
        !std::isfinite(cfg.imex_solver_tol) || cfg.imex_solver_max_iter < 0) {
      amrex::Abort("euler.bdltv20_paper_t1_s2 requires positive finite GMRES tolerance and nonnegative max_iter.");
    }
    if (!(cfg.bdltv20_paper_t1_s2_dt < amrex::Real(0.0)) &&
        (!(cfg.bdltv20_paper_t1_s2_dt > amrex::Real(0.0)) ||
         !std::isfinite(cfg.bdltv20_paper_t1_s2_dt))) {
      amrex::Abort("euler.bdltv20_paper_t1_s2_dt must be negative for CFL mode or finite positive.");
    }
    if (cfg.bdltv20_paper_pressure_solver != "gmres") {
      amrex::Abort("euler.bdltv20_paper_pressure_solver must be gmres.");
    }
    if (!(cfg.bdltv20_paper_epsilon > amrex::Real(0.0)) ||
        !std::isfinite(cfg.bdltv20_paper_epsilon)) {
      amrex::Abort("euler.bdltv20_paper_epsilon must be finite and positive.");
    }
    if (cfg.bdltv20_paper_t1_s2_max_steps < 1) {
      amrex::Abort("euler.bdltv20_paper_t1_s2_max_steps must be at least one.");
    }
  }

  if (cfg.method == MethodKind::Imex) {
    if (cfg.spatial_order != 2) {
      amrex::Abort("euler.method=imex requires euler.spatial_order=2 in this code pack.");
    }
    if (cfg.imex_pressure_solver_path != "host_single_level_gmres") {
      amrex::Abort("IMEX rows require euler.imex_pressure_solver_path=host_single_level_gmres.");
    }
    if (cfg.imex_pressure_stabilization != "off") {
      amrex::Abort("IMEX rows in this code pack use euler.imex_pressure_stabilization=off.");
    }
    if (cfg.imex_predictor_dissipation != ImexPredictorDissipationKind::Material) {
      amrex::Abort("IMEX rows in this code pack use euler.imex_predictor_dissipation=material.");
    }
    if (cfg.slope_limiter != SlopeLimiterKind::Minmod) {
      amrex::Abort("IMEX rows in this code pack use euler.slope_limiter=minmod.");
    }
    if (!has_base_bdltv20_paper_controls(cfg)) {
      amrex::Abort("IMEX rows require material dissipation, acoustic startup off, acoustic cap zero, and no pressure stabilization.");
    }
  }

  amrex::Vector<amrex::Real> velocity{cfg.velocity_x, cfg.velocity_y};
  amrex::Vector<amrex::Real> center{cfg.blob_cx, cfg.blob_cy};
  amrex::Vector<amrex::Real> vortex_center{cfg.vortex_cx, cfg.vortex_cy};
  pp_euler.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
  pp_euler.queryarr("blob_center", center, 0, AMREX_SPACEDIM);
  pp_euler.queryarr("vortex_center", vortex_center, 0, AMREX_SPACEDIM);
  cfg.velocity_x = velocity[0];
  cfg.velocity_y = velocity[1];
  cfg.blob_cx = center[0];
  cfg.blob_cy = center[1];
  cfg.vortex_cx = vortex_center[0];
  cfg.vortex_cy = vortex_center[1];
  pp_euler.query("density_inner", cfg.density_inner);
  pp_euler.query("density_outer", cfg.density_outer);
  pp_euler.query("pressure", cfg.pressure);
  pp_euler.query("mach", cfg.mach);
  pp_euler.query("blob_radius", cfg.blob_radius);
  pp_euler.query("isentropic_vortex_strength",
                 cfg.isentropic_vortex_strength);
  pp_euler.query("riemann_interface_x", cfg.riemann_interface_x);
  pp_euler.query("riemann_interface_y", cfg.riemann_interface_y);
  pp_euler.query("final_csv", cfg.final_csv);
  pp.query("final_csv", cfg.final_csv);

  amrex::ParmParse pp_prob("prob");
  pp_prob.query("test", cfg.toro_test);
  pp_prob.query("gamma", cfg.gamma);
  pp_prob.query("x0", cfg.riemann_interface_x);
  pp_prob.query("y0", cfg.riemann_interface_y);

  if (cfg.bdltv20_paper_t1_s2 == "toro_x_exact_dirichlet_y_periodic_2d") {
    if (cfg.problem != ProblemKind::Toro1 || cfg.is_periodic[0] != 0 ||
        cfg.is_periodic[1] != 1 || cfg.field_boundary != "exact_dirichlet" ||
        cfg.n_cell[0] < 16 || cfg.n_cell[1] < 2) {
      amrex::Abort("BDLTV20 Toro x-exact mode requires problem=toro1, periodicity=0 1, exact Dirichlet, and at least 16x2 cells.");
    }
  }
  if (cfg.bdltv20_paper_t1_s2 == "toro_xy_exact_dirichlet_2d") {
    if (cfg.problem != ProblemKind::Toro1 || cfg.is_periodic[0] != 0 ||
        cfg.is_periodic[1] != 0 || cfg.field_boundary != "exact_dirichlet" ||
        cfg.n_cell[0] < 16 || cfg.n_cell[1] < 2) {
      amrex::Abort("BDLTV20 Toro xy-exact mode requires problem=toro1, periodicity=0 0, exact Dirichlet, and at least 16x2 cells.");
    }
  }
  if (cfg.bdltv20_paper_t1_s2 == "gresho_exact_dirichlet_2d") {
    if (cfg.problem != ProblemKind::GreshoVortex || cfg.is_periodic[0] != 0 ||
        cfg.is_periodic[1] != 0 || cfg.field_boundary != "exact_dirichlet" ||
        cfg.n_cell[0] < 8 || cfg.n_cell[1] < 8) {
      amrex::Abort("BDLTV20 Gresho mode requires problem=gresho_vortex, periodicity=0 0, exact Dirichlet, and at least 8x8 cells.");
    }
  }
  if (cfg.bdltv20_paper_t1_s2 == "isentropic_vortex_periodic_2d") {
    if (cfg.problem != ProblemKind::IsentropicVortex || cfg.is_periodic[0] != 1 ||
        cfg.is_periodic[1] != 1 || cfg.n_cell[0] < 8 || cfg.n_cell[1] < 8) {
      amrex::Abort("BDLTV20 isentropic-vortex mode requires problem=isentropic_vortex, periodicity=1 1, and at least 8x8 cells.");
    }
  }
  if (cfg.bdltv20_paper_t1_s2 == "advection_blob_periodic_2d") {
    if (cfg.problem != ProblemKind::AdvectionBlob || cfg.is_periodic[0] != 1 ||
        cfg.is_periodic[1] != 1 || cfg.n_cell[0] < 8 || cfg.n_cell[1] < 8) {
      amrex::Abort("BDLTV20 advection-blob mode requires problem=advection_blob, periodicity=1 1, and at least 8x8 cells.");
    }
  }

  if (cfg.problem == ProblemKind::ShockDensityBubble2D ||
      cfg.problem == ProblemKind::ShockDensityBubbleCylindrical) {
    if (!max_step_set && cfg.max_step == 8) {
      cfg.max_step = 200000;
    }
    if (!grid_size_set && cfg.n_cell[0] == 64 && cfg.n_cell[1] == 64) {
      cfg.n_cell[0] = 160;
      cfg.n_cell[1] = 40;
    }
    if (cfg.prob_lo[0] == amrex::Real(0.0) && cfg.prob_lo[1] == amrex::Real(0.0) &&
        cfg.prob_hi[0] == amrex::Real(1.0) && cfg.prob_hi[1] == amrex::Real(1.0)) {
      cfg.prob_hi[0] = amrex::Real(2.0);
      cfg.prob_hi[1] = amrex::Real(0.5);
    }
    if (cfg.stop_time == amrex::Real(0.05)) {
      cfg.stop_time = amrex::Real(0.3);
    }
    if (cfg.cfl == amrex::Real(0.4)) {
      cfg.cfl = amrex::Real(0.45);
    }
    if (cfg.spatial_order == 1) {
      cfg.spatial_order = 2;
    }
    cfg.is_periodic[0] = 0;
    cfg.is_periodic[1] = 0;

    const bool valid_domain =
        std::abs(cfg.prob_lo[0] - amrex::Real(0.0)) < amrex::Real(1.0e-14) &&
        std::abs(cfg.prob_hi[0] - amrex::Real(2.0)) < amrex::Real(1.0e-14) &&
        std::abs(cfg.prob_lo[1] - amrex::Real(0.0)) < amrex::Real(1.0e-14) &&
        std::abs(cfg.prob_hi[1] - amrex::Real(0.5)) < amrex::Real(1.0e-14);
    if (!valid_domain) {
      amrex::Abort("shock_density_bubble requires domain x=[0,2], y=[0,0.5].");
    }
    if (std::abs(cfg.gamma - amrex::Real(1.4)) > amrex::Real(1.0e-14)) {
      amrex::Abort("shock_density_bubble is a same-gamma ideal-gas test with gamma=1.4.");
    }
    if (cfg.n_cell[0] < 16 || cfg.n_cell[1] < 4) {
      amrex::Abort("shock_density_bubble requires at least 16x4 cells.");
    }
    if (!(cfg.stop_time > amrex::Real(0.0)) || !std::isfinite(cfg.stop_time)) {
      amrex::Abort("shock_density_bubble requires positive finite final_time/stop_time.");
    }
    if (!(cfg.cfl > amrex::Real(0.0)) || !std::isfinite(cfg.cfl)) {
      amrex::Abort("shock_density_bubble requires positive finite CFL.");
    }
    if (cfg.method == MethodKind::Explicit) {
      if (cfg.spatial_order != 2) {
        amrex::Abort("shock_density_bubble explicit rows require order/spatial_order=2.");
      }
      if (cfg.riemann != RiemannKind::Hllc &&
          cfg.riemann != RiemannKind::XieAmHllcP) {
        amrex::Abort("shock_density_bubble explicit rows support hllc or xie_am_hllc_p.");
      }
    } else {
      if (cfg.problem == ProblemKind::ShockDensityBubbleCylindrical) {
        amrex::Abort("cylindrical shock_density_bubble IMEX is not part of the cleaned Cartesian code pack.");
      }
      if (!is_bdltv20_t1_s2_source_map_form(cfg.imex_form) || cfg.spatial_order != 2) {
        amrex::Abort("shock_density_bubble IMEX rows require bdltv20_t1_s2_source_map_picard and spatial_order=2.");
      }
    }
  }

  return cfg;
}

amrex::Geometry make_geometry(const RunConfig& cfg)
{
  const amrex::IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
  const amrex::IntVect dom_hi(AMREX_D_DECL(cfg.n_cell[0] - 1, cfg.n_cell[1] - 1, 0));
  const amrex::Box domain(dom_lo, dom_hi);
  const amrex::RealBox real_box(cfg.prob_lo.data(), cfg.prob_hi.data());
  return amrex::Geometry(domain, &real_box, amrex::CoordSys::cartesian, cfg.is_periodic.data());
}
