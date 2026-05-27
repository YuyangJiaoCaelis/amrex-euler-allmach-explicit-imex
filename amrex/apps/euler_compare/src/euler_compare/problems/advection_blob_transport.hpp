// Passive density-blob transport helper.
void advect_blob_one_step(amrex::MultiFab& state, const amrex::Geometry& geom, const RunConfig& cfg, amrex::Real dt)
{
  BL_PROFILE("advect_blob_one_step()");

  amrex::MultiFab next(state.boxArray(), state.DistributionMap(), NCons, 0);
  state.FillBoundary(geom.periodicity());

  const auto dx = geom.CellSizeArray();
  const amrex::Real u = cfg.velocity_x;
  const amrex::Real v = cfg.velocity_y;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const src = state.const_array(mfi);
    auto const dst = next.array(mfi);
    amrex::ParallelFor(box, static_cast<int>(NCons), [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      const amrex::Real x_grad =
          u >= amrex::Real(0.0) ? (src(i, j, k, n) - src(i - 1, j, k, n)) / dx[0]
                                : (src(i + 1, j, k, n) - src(i, j, k, n)) / dx[0];
      const amrex::Real y_grad =
          v >= amrex::Real(0.0) ? (src(i, j, k, n) - src(i, j - 1, k, n)) / dx[1]
                                : (src(i, j + 1, k, n) - src(i, j, k, n)) / dx[1];
      dst(i, j, k, n) = src(i, j, k, n) - dt * (u * x_grad + v * y_grad);
    });
  }

  amrex::MultiFab::Copy(state, next, 0, 0, NCons, 0);
}
