// Physical ghost-cell fill policies.
void fill_physical_outflow_ghosts(amrex::MultiFab& state, const amrex::Geometry& geom)
{
  state.FillBoundary(geom.periodicity());

  const amrex::Box domain = geom.Domain();
  const auto dom_lo = domain.smallEnd();
  const auto dom_hi = domain.bigEnd();
  const bool periodic_x = geom.isPeriodic(0);
  const bool periodic_y = geom.isPeriodic(1);

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box box = mfi.growntilebox(state.nGrowVect());
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, static_cast<int>(NCons), [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      const bool x_physical = !periodic_x && (i < dom_lo[0] || i > dom_hi[0]);
      const bool y_physical = !periodic_y && (j < dom_lo[1] || j > dom_hi[1]);
      if (x_physical || y_physical) {
        const int ii = x_physical ? amrex::min(amrex::max(i, dom_lo[0]), dom_hi[0]) : i;
        const int jj = y_physical ? amrex::min(amrex::max(j, dom_lo[1]), dom_hi[1]) : j;
        arr(i, j, k, n) = arr(ii, jj, k, n);
      }
    });
  }
}

void fill_shock_density_bubble_ghosts(amrex::MultiFab& state,
                                       const amrex::Geometry& geom,
                                       const RunConfig& cfg)
{
  state.FillBoundary(geom.periodicity());

  const amrex::Box domain = geom.Domain();
  const auto dom_lo = domain.smallEnd();
  const auto dom_hi = domain.bigEnd();
  const amrex::Real gamma = cfg.gamma;
  const ConservedState left_state =
      to_conserved(shock_density_bubble_post_shock_state(), gamma);

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box box = mfi.growntilebox(state.nGrowVect());
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const bool outside_x = i < dom_lo[0] || i > dom_hi[0];
      const bool outside_y = j < dom_lo[1] || j > dom_hi[1];
      if (!outside_x && !outside_y) {
        return;
      }

      if (i < dom_lo[0]) {
        store_cons(arr, i, j, k, left_state);
        return;
      }

      const int ii = i > dom_hi[0] ? dom_hi[0] : i;
      int jj = j;
      bool reflect_y = false;
      if (j < dom_lo[1]) {
        jj = dom_lo[1] + (dom_lo[1] - j - 1);
        reflect_y = true;
      } else if (j > dom_hi[1]) {
        jj = dom_hi[1];
      }
      jj = amrex::min(amrex::max(jj, dom_lo[1]), dom_hi[1]);

      ConservedState copied = load_cons(arr, ii, jj, k);
      if (reflect_y) {
        copied.my = -copied.my;
      }
      store_cons(arr, i, j, k, copied);
    });
  }
}

void apply_shock_density_bubble_cylindrical_source(amrex::MultiFab& state,
                                                   const amrex::Geometry& geom,
                                                   const RunConfig& cfg,
                                                   amrex::Real dt)
{
  if (!shock_density_bubble_uses_cylindrical_source(cfg)) {
    return;
  }

  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real gm1 = gamma - amrex::Real(1.0);
  const amrex::Real half_dt = amrex::Real(0.5) * dt;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real r =
          plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
      if (!(r > amrex::Real(0.0))) {
        return;
      }

      const ConservedState q0 = load_cons(arr, i, j, k);
      const amrex::Real inv_rho0 =
          q0.rho > amrex::Real(0.0) ? amrex::Real(1.0) / q0.rho
                                    : amrex::Real(0.0);
      const amrex::Real u0 = q0.mx * inv_rho0;
      const amrex::Real v0 = q0.my * inv_rho0;
      const amrex::Real p0 =
          gm1 * (q0.E - amrex::Real(0.5) * (q0.mx * q0.mx + q0.my * q0.my) * inv_rho0);

      ConservedState qs = q0;
      qs.rho = q0.rho - half_dt * q0.my / r;
      qs.mx = q0.mx - half_dt * q0.rho * u0 * v0 / r;
      qs.my = q0.my - half_dt * q0.rho * v0 * v0 / r;
      qs.E = q0.E - half_dt * v0 * (q0.E + p0) / r;

      const amrex::Real inv_rhos =
          qs.rho > amrex::Real(0.0) ? amrex::Real(1.0) / qs.rho
                                    : amrex::Real(0.0);
      const amrex::Real us = qs.mx * inv_rhos;
      const amrex::Real vs = qs.my * inv_rhos;
      const amrex::Real ps =
          gm1 * (qs.E - amrex::Real(0.5) * (qs.mx * qs.mx + qs.my * qs.my) * inv_rhos);

      ConservedState q1 = q0;
      q1.rho = q0.rho - dt * qs.my / r;
      q1.mx = q0.mx - dt * qs.rho * us * vs / r;
      q1.my = q0.my - dt * qs.rho * vs * vs / r;
      q1.E = q0.E - dt * vs * (qs.E + ps) / r;
      store_cons(arr, i, j, k, q1);
    });
  }
}

void fill_gresho_exact_dirichlet_ghosts(amrex::MultiFab& state,
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
        store_cons(arr, i, j, k, to_conserved(gresho_state(x, y, cfg), gamma));
      }
    });
  }
}

bool use_gresho_exact_dirichlet_boundary(const RunConfig& cfg,
                                         const amrex::Geometry& geom)
{
  return cfg.problem == ProblemKind::GreshoVortex &&
         cfg.field_boundary == "exact_dirichlet" &&
         (!geom.isPeriodic(0) || !geom.isPeriodic(1));
}

bool use_toro_x_exact_dirichlet_boundary(const RunConfig& cfg,
                                         const amrex::Geometry& geom)
{
  return cfg.problem == ProblemKind::Toro1 &&
         cfg.field_boundary == "exact_dirichlet" &&
         (!geom.isPeriodic(0) || !geom.isPeriodic(1));
}

bool use_riemann_quadrant_exact_dirichlet_boundary(const RunConfig& cfg,
                                                   const amrex::Geometry& geom)
{
  return cfg.problem == ProblemKind::RiemannQuadrant &&
         cfg.field_boundary == "exact_dirichlet" &&
         (!geom.isPeriodic(0) || !geom.isPeriodic(1));
}

bool use_advection_blob_exact_dirichlet_boundary(const RunConfig& cfg,
                                                 const amrex::Geometry& geom)
{
  return cfg.problem == ProblemKind::AdvectionBlob &&
         cfg.field_boundary == "exact_dirichlet" &&
         (!geom.isPeriodic(0) || !geom.isPeriodic(1));
}

bool use_shock_density_bubble_boundary(const RunConfig& cfg)
{
  return is_shock_density_bubble_case(cfg.problem);
}

void fill_problem_ghosts(amrex::MultiFab& state, const amrex::Geometry& geom,
                         const RunConfig& cfg,
                         amrex::Real time = amrex::Real(0.0))
{
  if (use_shock_density_bubble_boundary(cfg)) {
    fill_shock_density_bubble_ghosts(state, geom, cfg);
  } else if (use_gresho_exact_dirichlet_boundary(cfg, geom)) {
    fill_gresho_exact_dirichlet_ghosts(state, geom, cfg);
  } else if (use_toro_x_exact_dirichlet_boundary(cfg, geom)) {
    fill_toro_x_exact_dirichlet_ghosts(state, geom, cfg, time);
  } else if (use_riemann_quadrant_exact_dirichlet_boundary(cfg, geom)) {
    fill_riemann_quadrant_exact_dirichlet_ghosts(state, geom, cfg);
  } else if (use_advection_blob_exact_dirichlet_boundary(cfg, geom)) {
    fill_advection_blob_exact_dirichlet_ghosts(state, geom, cfg, time);
  } else {
    fill_physical_outflow_ghosts(state, geom);
  }
}

void fill_scalar_outflow_ghosts(amrex::MultiFab& field, const amrex::Geometry& geom)
{
  field.FillBoundary(geom.periodicity());

  const amrex::Box domain = geom.Domain();
  const auto dom_lo = domain.smallEnd();
  const auto dom_hi = domain.bigEnd();
  const bool periodic_x = geom.isPeriodic(0);
  const bool periodic_y = geom.isPeriodic(1);
  const int ncomp = field.nComp();

  for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
    const amrex::Box box = mfi.growntilebox(field.nGrowVect());
    auto const arr = field.array(mfi);
    amrex::ParallelFor(box, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      const bool x_physical = !periodic_x && (i < dom_lo[0] || i > dom_hi[0]);
      const bool y_physical = !periodic_y && (j < dom_lo[1] || j > dom_hi[1]);
      if (x_physical || y_physical) {
        const int ii = x_physical ? amrex::min(amrex::max(i, dom_lo[0]), dom_hi[0]) : i;
        const int jj = y_physical ? amrex::min(amrex::max(j, dom_lo[1]), dom_hi[1]) : j;
        arr(i, j, k, n) = arr(ii, jj, k, n);
      }
    });
  }
}
