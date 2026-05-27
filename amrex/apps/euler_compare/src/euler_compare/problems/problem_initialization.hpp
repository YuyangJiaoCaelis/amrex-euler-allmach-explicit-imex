// Initial conditions and exact-state fills used by report problems.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
shock_density_bubble_shock_x() noexcept
{
  return amrex::Real(0.2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
shock_density_bubble_center_x() noexcept
{
  return amrex::Real(0.5);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
shock_density_bubble_center_y() noexcept
{
  return amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
shock_density_bubble_radius() noexcept
{
  return amrex::Real(0.2);
}

bool is_shock_density_bubble_case(ProblemKind problem) noexcept
{
  return problem == ProblemKind::ShockDensityBubble2D ||
         problem == ProblemKind::ShockDensityBubbleCylindrical;
}

bool shock_density_bubble_uses_cylindrical_source(const RunConfig& cfg) noexcept
{
  return cfg.problem == ProblemKind::ShockDensityBubbleCylindrical;
}

const char* shock_density_bubble_source_case_id(const RunConfig& cfg) noexcept
{
  return shock_density_bubble_uses_cylindrical_source(cfg)
             ? "same_gamma_shock_density_bubble_clawpack_cylindrical_v1"
             : "same_gamma_shock_density_bubble_clawpack_style_v1";
}

const char* shock_density_bubble_validation_claim(const RunConfig& cfg) noexcept
{
  return shock_density_bubble_uses_cylindrical_source(cfg)
             ? "qualitative_project_clawpack_style_cylindrical_high_speed_2d_only"
             : "qualitative_project_high_speed_2d_only";
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState
shock_density_bubble_post_shock_state() noexcept
{
  return PrimitiveState{amrex::Real(2.8181818181818183),
                        amrex::Real(1.6064386578049976),
                        amrex::Real(0.0),
                        amrex::Real(5.0)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState
shock_density_bubble_ambient_state(amrex::Real rho) noexcept
{
  return PrimitiveState{rho, amrex::Real(0.0), amrex::Real(0.0), amrex::Real(1.0)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
shock_density_bubble_fraction(amrex::Real xc, amrex::Real yc,
                              amrex::Real dx, amrex::Real dy) noexcept
{
  amrex::Real inside = amrex::Real(0.0);
  for (int sy = 0; sy < 4; ++sy) {
    for (int sx = 0; sx < 4; ++sx) {
      const amrex::Real x =
          xc + ((static_cast<amrex::Real>(sx) + amrex::Real(0.5)) / amrex::Real(4.0) -
                amrex::Real(0.5)) * dx;
      const amrex::Real y =
          yc + ((static_cast<amrex::Real>(sy) + amrex::Real(0.5)) / amrex::Real(4.0) -
                amrex::Real(0.5)) * dy;
      const amrex::Real rx = x - shock_density_bubble_center_x();
      const amrex::Real ry = y - shock_density_bubble_center_y();
      if (x >= shock_density_bubble_shock_x() &&
          rx * rx + ry * ry <=
              shock_density_bubble_radius() * shock_density_bubble_radius()) {
        inside += amrex::Real(1.0);
      }
    }
  }
  return inside / amrex::Real(16.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState
shock_density_bubble_cell_state(amrex::Real x, amrex::Real y,
                                amrex::Real dx, amrex::Real dy) noexcept
{
  if (x < shock_density_bubble_shock_x()) {
    return shock_density_bubble_post_shock_state();
  }
  const amrex::Real f = shock_density_bubble_fraction(x, y, dx, dy);
  const amrex::Real rho = f * amrex::Real(0.1) +
                          (amrex::Real(1.0) - f) * amrex::Real(1.0);
  return shock_density_bubble_ambient_state(rho);
}

void initialize_shock_density_bubble(amrex::MultiFab& state,
                                     const amrex::Geometry& geom,
                                     const RunConfig& cfg)
{
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x =
          plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
      const amrex::Real y =
          plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
      store_cons(arr, i, j, k,
                 to_conserved(shock_density_bubble_cell_state(x, y, dx[0], dx[1]),
                              gamma));
    });
  }
}

void initialize_blob(amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg)
{
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real u = cfg.velocity_x;
  const amrex::Real v = cfg.velocity_y;
  const amrex::Real p = cfg.pressure;
  const amrex::Real rho_inner = cfg.density_inner;
  const amrex::Real rho_outer = cfg.density_outer;
  const amrex::Real cx = cfg.blob_cx;
  const amrex::Real cy = cfg.blob_cy;
  const amrex::Real radius = cfg.blob_radius;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x = plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
      const amrex::Real y = plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
      const amrex::Real r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
      const amrex::Real rho = r2 <= radius * radius ? rho_inner : rho_outer;
      arr(i, j, k, Rho) = rho;
      arr(i, j, k, Mx) = rho * u;
      arr(i, j, k, My) = rho * v;
      arr(i, j, k, E) = p / (gamma - amrex::Real(1.0)) + amrex::Real(0.5) * rho * (u * u + v * v);
    });
  }
}

amrex::Real periodic_signed_distance(amrex::Real x, amrex::Real center,
                                     amrex::Real lo, amrex::Real hi,
                                     bool periodic)
{
  amrex::Real d = x - center;
  if (periodic) {
    const amrex::Real length = hi - lo;
    if (length > amrex::Real(0.0)) {
      d -= std::round(d / length) * length;
    }
  }
  return d;
}

amrex::Real wrap_periodic_coordinate(amrex::Real x, amrex::Real lo, amrex::Real hi)
{
  const amrex::Real length = hi - lo;
  if (length <= amrex::Real(0.0)) {
    return x;
  }
  amrex::Real wrapped = lo + std::fmod(x - lo, length);
  if (wrapped < lo) {
    wrapped += length;
  }
  return wrapped;
}

PrimitiveState exact_advection_blob_state(amrex::Real x, amrex::Real y,
                                          const amrex::Geometry& geom,
                                          const RunConfig& cfg,
                                          amrex::Real time)
{
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  amrex::Real cx = cfg.blob_cx + cfg.velocity_x * time;
  amrex::Real cy = cfg.blob_cy + cfg.velocity_y * time;
  if (geom.isPeriodic(0)) {
    cx = wrap_periodic_coordinate(cx, plo[0], phi[0]);
  }
  if (geom.isPeriodic(1)) {
    cy = wrap_periodic_coordinate(cy, plo[1], phi[1]);
  }
  const amrex::Real dx =
      periodic_signed_distance(x, cx, plo[0], phi[0], geom.isPeriodic(0));
  const amrex::Real dy =
      periodic_signed_distance(y, cy, plo[1], phi[1], geom.isPeriodic(1));
  const amrex::Real r2 = dx * dx + dy * dy;
  const amrex::Real rho =
      r2 <= cfg.blob_radius * cfg.blob_radius ? cfg.density_inner : cfg.density_outer;
  return PrimitiveState{rho, cfg.velocity_x, cfg.velocity_y, cfg.pressure};
}

void fill_advection_blob_exact_dirichlet_ghosts(amrex::MultiFab& state,
                                                const amrex::Geometry& geom,
                                                const RunConfig& cfg,
                                                amrex::Real time)
{
  state.FillBoundary(geom.periodicity());

  const amrex::Box domain = geom.Domain();
  const auto dom_lo = domain.smallEnd();
  const auto dom_hi = domain.bigEnd();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const bool periodic_x = geom.isPeriodic(0);
  const bool periodic_y = geom.isPeriodic(1);
  const int nx = dom_hi[0] - dom_lo[0] + 1;
  const int ny = dom_hi[1] - dom_lo[1] + 1;
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real u = cfg.velocity_x;
  const amrex::Real v = cfg.velocity_y;
  const amrex::Real p = cfg.pressure;
  const amrex::Real rho_inner = cfg.density_inner;
  const amrex::Real rho_outer = cfg.density_outer;
  const amrex::Real cx = cfg.blob_cx + cfg.velocity_x * time;
  const amrex::Real cy = cfg.blob_cy + cfg.velocity_y * time;
  const amrex::Real radius = cfg.blob_radius;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box box = mfi.growntilebox(state.nGrowVect());
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const bool x_physical = !periodic_x && (i < dom_lo[0] || i > dom_hi[0]);
      const bool y_physical = !periodic_y && (j < dom_lo[1] || j > dom_hi[1]);
      if (x_physical || y_physical) {
        int ii = i;
        int jj = j;
        if (periodic_x) {
          const int offset = (i - dom_lo[0]) % nx;
          ii = dom_lo[0] + (offset < 0 ? offset + nx : offset);
        }
        if (periodic_y) {
          const int offset = (j - dom_lo[1]) % ny;
          jj = dom_lo[1] + (offset < 0 ? offset + ny : offset);
        }
        const amrex::Real x =
            plo[0] + (static_cast<amrex::Real>(ii) + amrex::Real(0.5)) * dx[0];
        const amrex::Real y =
            plo[1] + (static_cast<amrex::Real>(jj) + amrex::Real(0.5)) * dx[1];
        const amrex::Real r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
        const amrex::Real rho = r2 <= radius * radius ? rho_inner : rho_outer;
        store_cons(arr, i, j, k,
                   to_conserved(PrimitiveState{rho, u, v, p}, gamma));
      }
    });
  }
}

void initialize_toro(amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg)
{
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real x0 = cfg.riemann_interface_x;
  const ToroState ts = toro_state(cfg.toro_test);

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x = plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
      const bool left = x < x0;
      const PrimitiveState q{left ? ts.rho_l : ts.rho_r, left ? ts.u_l : ts.u_r, left ? ts.v_l : ts.v_r,
                             left ? ts.p_l : ts.p_r};
      store_cons(arr, i, j, k, to_conserved(q, gamma));
    });
  }
}

void initialize_riemann_quadrant(amrex::MultiFab& state,
                                 const amrex::Geometry& geom,
                                 const RunConfig& cfg)
{
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x =
          plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
      const amrex::Real y =
          plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
      store_cons(arr, i, j, k,
                 to_conserved(riemann_quadrant_state(x, y, cfg), gamma));
    });
  }
}

void fill_toro_x_exact_dirichlet_ghosts(amrex::MultiFab& state,
                                        const amrex::Geometry& geom,
                                        const RunConfig& cfg,
                                        amrex::Real time)
{
  state.FillBoundary(geom.periodicity());

  const amrex::Box domain = geom.Domain();
  const auto dom_lo = domain.smallEnd();
  const auto dom_hi = domain.bigEnd();
  const bool periodic_y = geom.isPeriodic(1);
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real interface_x = cfg.riemann_interface_x;
  const ToroState ts = toro_state(cfg.toro_test);
  const PrimitiveState left_primitive{ts.rho_l, ts.u_l, ts.v_l, ts.p_l};
  const PrimitiveState right_primitive{ts.rho_r, ts.u_r, ts.v_r, ts.p_r};

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box box = mfi.growntilebox(state.nGrowVect());
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (i < dom_lo[0] || i > dom_hi[0] ||
          (!periodic_y && (j < dom_lo[1] || j > dom_hi[1]))) {
        const amrex::Real x =
            plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
        const PrimitiveState exact =
            exact_riemann_sample(x, time, left_primitive, right_primitive, gamma,
                                 interface_x);
        store_cons(arr, i, j, k, to_conserved(exact, gamma));
      }
    });
  }
}

void fill_riemann_quadrant_exact_dirichlet_ghosts(amrex::MultiFab& state,
                                                  const amrex::Geometry& geom,
                                                  const RunConfig& cfg)
{
  state.FillBoundary(geom.periodicity());

  const amrex::Box domain = geom.Domain();
  const auto dom_lo = domain.smallEnd();
  const auto dom_hi = domain.bigEnd();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const bool periodic_x = geom.isPeriodic(0);
  const bool periodic_y = geom.isPeriodic(1);
  const amrex::Real gamma = cfg.gamma;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box box = mfi.growntilebox(state.nGrowVect());
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const bool x_physical = !periodic_x && (i < dom_lo[0] || i > dom_hi[0]);
      const bool y_physical = !periodic_y && (j < dom_lo[1] || j > dom_hi[1]);
      if (x_physical || y_physical) {
        const amrex::Real x =
            plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
        const amrex::Real y =
            plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
        store_cons(arr, i, j, k,
                   to_conserved(riemann_quadrant_state(x, y, cfg), gamma));
      }
    });
  }
}

void initialize_gresho(amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg)
{
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x = plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
      const amrex::Real y = plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
      store_cons(arr, i, j, k, to_conserved(gresho_state(x, y, cfg), gamma));
    });
  }
}

void initialize_isentropic_vortex(amrex::MultiFab& state,
                                  const amrex::Geometry& geom,
                                  const RunConfig& cfg)
{
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const amrex::Real gamma = cfg.gamma;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x =
          plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
      const amrex::Real y =
          plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
      store_cons(arr, i, j, k,
                 to_conserved(isentropic_vortex_state(
                                  x, y, cfg, plo[0], phi[0], plo[1], phi[1],
                                  amrex::Real(0.0)),
                              gamma));
    });
  }
}

void initialize_state(amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg)
{
  if (cfg.problem == ProblemKind::Toro1) {
    initialize_toro(state, geom, cfg);
  } else if (cfg.problem == ProblemKind::GreshoVortex) {
    initialize_gresho(state, geom, cfg);
  } else if (cfg.problem == ProblemKind::IsentropicVortex) {
    initialize_isentropic_vortex(state, geom, cfg);
  } else if (cfg.problem == ProblemKind::RiemannQuadrant) {
    initialize_riemann_quadrant(state, geom, cfg);
  } else if (is_shock_density_bubble_case(cfg.problem)) {
    initialize_shock_density_bubble(state, geom, cfg);
  } else {
    initialize_blob(state, geom, cfg);
  }
}
