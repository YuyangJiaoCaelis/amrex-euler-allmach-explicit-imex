// BDLTV20/Toro-Vazquez source-map and pressure-correction helper routines.
void compute_bdltv20_pressure_split_star(const amrex::MultiFab& state, amrex::MultiFab& star,
                                         const amrex::Geometry& geom, const RunConfig& cfg, amrex::Real dt)
{
  const auto dx = geom.CellSizeArray();
  const amrex::Real inv_dx = amrex::Real(1.0) / dx[0];
  const amrex::Real inv_dy = amrex::Real(1.0) / dx[1];
  const amrex::Real gamma = cfg.gamma;
  const SlopeLimiterKind slope_limiter = cfg.slope_limiter;
  const ImexPredictorDissipationKind predictor_dissipation = cfg.imex_predictor_dissipation;
  const bool use_eq87_89_conserved_minmod =
      is_bdltv20_t1_s2_source_map_form(cfg.imex_form) && cfg.spatial_order == 2;
  const int spatial_order = cfg.spatial_order;

  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const src = state.const_array(mfi);
    auto const dst = star.array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const ConservedState fx_l =
          bdltv20_pressure_split_lf_flux(
              reconstruct_bdltv20_star_state(src, i - 1, j, k, 0, +1,
                                             use_eq87_89_conserved_minmod, spatial_order, gamma,
                                             slope_limiter),
              reconstruct_bdltv20_star_state(src, i, j, k, 0, -1,
                                             use_eq87_89_conserved_minmod, spatial_order, gamma,
                                             slope_limiter),
              0, gamma, predictor_dissipation);
      const ConservedState fx_r =
          bdltv20_pressure_split_lf_flux(
              reconstruct_bdltv20_star_state(src, i, j, k, 0, +1,
                                             use_eq87_89_conserved_minmod, spatial_order, gamma,
                                             slope_limiter),
              reconstruct_bdltv20_star_state(src, i + 1, j, k, 0, -1,
                                             use_eq87_89_conserved_minmod, spatial_order, gamma,
                                             slope_limiter),
              0, gamma, predictor_dissipation);
      const ConservedState fy_b =
          bdltv20_pressure_split_lf_flux(
              reconstruct_bdltv20_star_state(src, i, j - 1, k, 1, +1,
                                             use_eq87_89_conserved_minmod, spatial_order, gamma,
                                             slope_limiter),
              reconstruct_bdltv20_star_state(src, i, j, k, 1, -1,
                                             use_eq87_89_conserved_minmod, spatial_order, gamma,
                                             slope_limiter),
              1, gamma, predictor_dissipation);
      const ConservedState fy_t =
          bdltv20_pressure_split_lf_flux(
              reconstruct_bdltv20_star_state(src, i, j, k, 1, +1,
                                             use_eq87_89_conserved_minmod, spatial_order, gamma,
                                             slope_limiter),
              reconstruct_bdltv20_star_state(src, i, j + 1, k, 1, -1,
                                             use_eq87_89_conserved_minmod, spatial_order, gamma,
                                             slope_limiter),
              1, gamma, predictor_dissipation);

      const ConservedState old = load_cons(src, i, j, k);
      ConservedState updated;
      updated.rho = old.rho - dt * ((fx_r.rho - fx_l.rho) * inv_dx + (fy_t.rho - fy_b.rho) * inv_dy);
      updated.mx = old.mx - dt * ((fx_r.mx - fx_l.mx) * inv_dx + (fy_t.mx - fy_b.mx) * inv_dy);
      updated.my = old.my - dt * ((fx_r.my - fx_l.my) * inv_dx + (fy_t.my - fy_b.my) * inv_dy);
      updated.E = old.E - dt * ((fx_r.E - fx_l.E) * inv_dx + (fy_t.E - fy_b.E) * inv_dy);
      store_cons(dst, i, j, k, updated);
    });
  }
  amrex::Gpu::streamSynchronize();
  fill_problem_ghosts(star, geom, cfg);
}

struct CoupledHostGrid {
  int nx = 0;
  int ny = 0;
  bool periodic_x = false;
  bool periodic_y = false;
  std::vector<ConservedState> cell;

  int row(int i, int j) const { return j * nx + i; }

  int map_i(int i) const
  {
    if (periodic_x) {
      return (i % nx + nx) % nx;
    }
    return std::min(std::max(i, 0), nx - 1);
  }

  int map_j(int j) const
  {
    if (periodic_y) {
      return (j % ny + ny) % ny;
    }
    return std::min(std::max(j, 0), ny - 1);
  }

  const ConservedState& at(int i, int j) const { return cell[row(map_i(i), map_j(j))]; }
};

ConservedState source_map_state_at_physical_boundary(const CoupledHostGrid& grid,
                                                     const RunConfig& cfg,
                                                     int i, int j)
{
  if (cfg.problem == ProblemKind::ShockDensityBubble2D) {
    if (!grid.periodic_x && i < 0) {
      return to_conserved(shock_density_bubble_post_shock_state(), cfg.gamma);
    }
    const int ii = grid.map_i(i);
    if (!grid.periodic_y && j < 0) {
      ConservedState reflected = grid.at(ii, 0);
      reflected.my = -reflected.my;
      return reflected;
    }
    return grid.at(ii, grid.map_j(j));
  }
  return grid.at(i, j);
}

CoupledHostGrid load_coupled_host_grid(const amrex::MultiFab& state, const amrex::Geometry& geom,
                                       const RunConfig& cfg)
{
  CoupledHostGrid grid;
  grid.nx = cfg.n_cell[0];
  grid.ny = cfg.n_cell[1];
  grid.periodic_x = geom.isPeriodic(0);
  grid.periodic_y = geom.isPeriodic(1);
  grid.cell.assign(static_cast<std::size_t>(grid.nx * grid.ny), ConservedState{});
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        grid.cell[grid.row(i, j)] = load_cons(arr, i, j, 0);
      }
    }
  }
  return grid;
}

struct CoupledFaceStats {
  int enthalpy_guard_count = 0;
  int enthalpy_nonfinite_count = 0;
  int flux_nonfinite_count = 0;
  amrex::Real flux_linf = amrex::Real(0.0);
};

void copy_pressure_vector_to_multifab(const std::vector<amrex::Real>& p, amrex::MultiFab& pressure,
                                      const CoupledHostGrid& grid, const amrex::Geometry& geom)
{
  for (amrex::MFIter mfi(pressure); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const p_arr = pressure.array(mfi);
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        p_arr(i, j, 0) = p[grid.row(i, j)];
      }
    }
  }
  fill_scalar_outflow_ghosts(pressure, geom);
}

struct SourceMapIndex {
  int index = 0;
  std::string boundary_label;
};

SourceMapIndex source_map_index(int index, int n, bool periodic)
{
  if (periodic) {
    const int mapped = (index % n + n) % n;
    return SourceMapIndex{mapped, mapped == index ? std::string{} : std::string("boundary_periodic_wrap")};
  }
  if (index < 0) {
    const int layer = std::min(-index, 2);
    return SourceMapIndex{std::min(-index - 1, n - 1),
                          layer == 1 ? "boundary_neumann_reflection_layer1"
                                     : "boundary_neumann_reflection_layer2"};
  }
  if (index >= n) {
    const int layer = std::min(index - n + 1, 2);
    return SourceMapIndex{std::max(n - (index - n + 1), 0),
                          layer == 1 ? "boundary_neumann_reflection_layer1"
                                     : "boundary_neumann_reflection_layer2"};
  }
  return SourceMapIndex{index, std::string{}};
}

struct SourceMapCellIndex {
  int i = 0;
  int j = 0;
  std::string boundary_label;
};

SourceMapCellIndex source_map_pressure_index(const CoupledHostGrid& grid, int i, int j)
{
  const SourceMapIndex mi = source_map_index(i, grid.nx, grid.periodic_x);
  const SourceMapIndex mj = source_map_index(j, grid.ny, grid.periodic_y);
  SourceMapCellIndex out;
  out.i = mi.index;
  out.j = mj.index;
  if (!mi.boundary_label.empty()) {
    out.boundary_label = mi.boundary_label;
  }
  if (!mj.boundary_label.empty()) {
    if (!out.boundary_label.empty()) {
      out.boundary_label += ';';
    }
    out.boundary_label += mj.boundary_label;
  }
  return out;
}

amrex::Real source_map_pressure_at(const CoupledHostGrid& grid, const std::vector<amrex::Real>& p,
                                   int i, int j)
{
  const SourceMapCellIndex mapped = source_map_pressure_index(grid, i, j);
  return p[static_cast<std::size_t>(grid.row(mapped.i, mapped.j))];
}

bool source_map_fixed_pressure_boundary_value(const CoupledHostGrid& grid,
                                              const RunConfig& cfg,
                                              int i, int j,
                                              amrex::Real& value)
{
  (void)j;
  if (cfg.problem == ProblemKind::ShockDensityBubble2D && !grid.periodic_x && i < 0) {
    value = shock_density_bubble_post_shock_state().p;
    return true;
  }
  return false;
}

amrex::Real source_map_pressure_at_physical_boundary(const CoupledHostGrid& grid,
                                                     const std::vector<amrex::Real>& p,
                                                     const RunConfig& cfg,
                                                     int i, int j)
{
  amrex::Real boundary_value = std::numeric_limits<amrex::Real>::quiet_NaN();
  if (source_map_fixed_pressure_boundary_value(grid, cfg, i, j, boundary_value)) {
    return boundary_value;
  }
  return source_map_pressure_at(grid, p, i, j);
}

using SourceMapMomentum = std::array<amrex::Real, 2>;

std::array<amrex::Real, 2> source_map_corrected_momentum(const CoupledHostGrid& grid,
                                                         const std::vector<amrex::Real>& p,
                                                         const RunConfig* cfg,
                                                         int i, int j, amrex::Real dt,
                                                         amrex::Real inv_dx,
                                                         amrex::Real inv_dy)
{
  const ConservedState& star = grid.at(i, j);
  const auto pressure_at = [&](int pi, int pj) {
    return cfg == nullptr ? source_map_pressure_at(grid, p, pi, pj)
                          : source_map_pressure_at_physical_boundary(grid, p, *cfg, pi, pj);
  };
  const amrex::Real grad_x =
      amrex::Real(0.5) * (pressure_at(i + 1, j) - pressure_at(i - 1, j)) *
      inv_dx;
  const amrex::Real grad_y =
      amrex::Real(0.5) * (pressure_at(i, j + 1) - pressure_at(i, j - 1)) *
      inv_dy;
  return {star.mx - dt * grad_x, star.my - dt * grad_y};
}

struct SourceMapFaceFluxDetail {
  int axis = -1;
  int left_i = -1;
  int left_j = -1;
  int right_i = -1;
  int right_j = -1;
  amrex::Real p_left = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real p_right = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real rho_left = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real rho_right = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real q_left = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real q_right = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real denom = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real scale = std::numeric_limits<amrex::Real>::quiet_NaN();
  int guard = 0;
  amrex::Real h_left = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real h_right = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real h_hat = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real m_hat = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real energy_flux = std::numeric_limits<amrex::Real>::quiet_NaN();
  int flux_nonfinite = 0;
};

SourceMapMomentum source_map_lagged_momentum_at(const CoupledHostGrid& grid,
                                                const std::vector<SourceMapMomentum>& momentum,
                                                int i, int j)
{
  return momentum[static_cast<std::size_t>(grid.row(grid.map_i(i), grid.map_j(j)))];
}

SourceMapMomentum source_map_lagged_momentum_at_physical_boundary(
    const CoupledHostGrid& grid,
    const std::vector<SourceMapMomentum>& momentum,
    const RunConfig& cfg,
    int i, int j)
{
  if (cfg.problem == ProblemKind::ShockDensityBubble2D) {
    if (!grid.periodic_x && i < 0) {
      const ConservedState left =
          to_conserved(shock_density_bubble_post_shock_state(), cfg.gamma);
      return SourceMapMomentum{left.mx, left.my};
    }
    const int ii = grid.map_i(i);
    if (!grid.periodic_y && j < 0) {
      SourceMapMomentum reflected = source_map_lagged_momentum_at(grid, momentum, ii, 0);
      reflected[1] = -reflected[1];
      return reflected;
    }
  }
  return source_map_lagged_momentum_at(grid, momentum, i, j);
}

std::vector<SourceMapMomentum> source_map_eq53_star_momentum_field(const CoupledHostGrid& grid)
{
  std::vector<SourceMapMomentum> momentum(static_cast<std::size_t>(grid.nx * grid.ny));
  for (int j = 0; j < grid.ny; ++j) {
    for (int i = 0; i < grid.nx; ++i) {
      const ConservedState& star = grid.at(i, j);
      momentum[static_cast<std::size_t>(grid.row(i, j))] = SourceMapMomentum{star.mx, star.my};
    }
  }
  return momentum;
}

std::vector<SourceMapMomentum> source_map_eq54_corrected_momentum_field(
    const CoupledHostGrid& grid, const std::vector<amrex::Real>& p, amrex::Real dt,
    amrex::Real inv_dx, amrex::Real inv_dy, const RunConfig* cfg = nullptr)
{
  std::vector<SourceMapMomentum> momentum(static_cast<std::size_t>(grid.nx * grid.ny));
  for (int j = 0; j < grid.ny; ++j) {
    for (int i = 0; i < grid.nx; ++i) {
      momentum[static_cast<std::size_t>(grid.row(i, j))] =
          source_map_corrected_momentum(grid, p, cfg, i, j, dt, inv_dx, inv_dy);
    }
  }
  return momentum;
}

int source_map_momentum_nonfinite_count(const std::vector<SourceMapMomentum>& momentum)
{
  int count = 0;
  for (const SourceMapMomentum& value : momentum) {
    if (!std::isfinite(value[0]) || !std::isfinite(value[1])) {
      ++count;
    }
  }
  return count;
}

SourceMapFaceFluxDetail source_map_face_flux_detail(const CoupledHostGrid& grid,
                                                    const std::vector<amrex::Real>& p,
                                                    const std::vector<SourceMapMomentum>& momentum,
                                                    int il, int jl, int ir, int jr, int axis,
                                                    const RunConfig& cfg, amrex::Real dt,
                                                    amrex::Real inv_dx,
                                                    amrex::Real inv_dy)
{
  const amrex::Real gamma_minus_one = cfg.gamma - amrex::Real(1.0);
  const SourceMapCellIndex left_map = source_map_pressure_index(grid, il, jl);
  const SourceMapCellIndex right_map = source_map_pressure_index(grid, ir, jr);
  const ConservedState left = source_map_state_at_physical_boundary(grid, cfg, il, jl);
  const ConservedState right = source_map_state_at_physical_boundary(grid, cfg, ir, jr);
  const auto m_left =
      source_map_lagged_momentum_at_physical_boundary(grid, momentum, cfg, il, jl);
  const auto m_right =
      source_map_lagged_momentum_at_physical_boundary(grid, momentum, cfg, ir, jr);

  SourceMapFaceFluxDetail detail;
  detail.axis = axis;
  detail.left_i = left_map.i;
  detail.left_j = left_map.j;
  detail.right_i = right_map.i;
  detail.right_j = right_map.j;
  detail.rho_left = left.rho;
  detail.rho_right = right.rho;
  detail.p_left = source_map_pressure_at_physical_boundary(grid, p, cfg, il, jl);
  detail.p_right = source_map_pressure_at_physical_boundary(grid, p, cfg, ir, jr);
  detail.h_left = cfg.gamma * detail.p_left / (gamma_minus_one * left.rho);
  detail.h_right = cfg.gamma * detail.p_right / (gamma_minus_one * right.rho);
  detail.q_left = axis == 0 ? m_left[0] : m_left[1];
  detail.q_right = axis == 0 ? m_right[0] : m_right[1];
  detail.denom = detail.q_left + detail.q_right;
  detail.scale = std::max({std::abs(detail.q_left), std::abs(detail.q_right), amrex::Real(1.0)});
  detail.guard = std::abs(detail.denom) <= amrex::Real(1.0e-14) * detail.scale ? 1 : 0;
  detail.h_hat =
      detail.guard != 0 ? amrex::Real(0.5) * (detail.h_left + detail.h_right)
                        : (detail.h_left * detail.q_left + detail.h_right * detail.q_right) /
                              detail.denom;
  detail.m_hat = amrex::Real(0.5) * (detail.q_left + detail.q_right);
  detail.energy_flux = detail.h_hat * detail.m_hat;
  detail.flux_nonfinite = std::isfinite(detail.energy_flux) ? 0 : 1;
  return detail;
}

amrex::Real source_map_face_energy_flux(const CoupledHostGrid& grid,
                                        const std::vector<amrex::Real>& p,
                                        const std::vector<SourceMapMomentum>& momentum,
                                        int il, int jl, int ir, int jr, int axis, const RunConfig& cfg,
                                        amrex::Real dt, amrex::Real inv_dx,
                                        amrex::Real inv_dy, CoupledFaceStats* stats)
{
  const SourceMapFaceFluxDetail detail =
      source_map_face_flux_detail(grid, p, momentum, il, jl, ir, jr, axis, cfg, dt,
                                  inv_dx, inv_dy);
  if (stats != nullptr) {
    stats->enthalpy_guard_count += detail.guard;
    if (!std::isfinite(detail.h_left) || !std::isfinite(detail.h_right) ||
        !std::isfinite(detail.h_hat) || !std::isfinite(detail.q_left) ||
        !std::isfinite(detail.q_right)) {
      ++stats->enthalpy_nonfinite_count;
    }
    if (detail.flux_nonfinite != 0) {
      ++stats->flux_nonfinite_count;
    } else {
      stats->flux_linf = std::max(stats->flux_linf, std::abs(detail.energy_flux));
    }
  }
  return detail.energy_flux;
}

amrex::Real conserved_component(const ConservedState& u, int component)
{
  switch (component) {
    case Rho:
      return u.rho;
    case Mx:
      return u.mx;
    case My:
      return u.my;
    case E:
      return u.E;
  }
  return std::numeric_limits<amrex::Real>::quiet_NaN();
}

struct SourceMapRowTerm {
  int col_i = -1;
  int col_j = -1;
  amrex::Real coeff = std::numeric_limits<amrex::Real>::quiet_NaN();
  std::string row_term_label = "unset";
};

struct SourceMapRhsTerm {
  amrex::Real value = std::numeric_limits<amrex::Real>::quiet_NaN();
  std::string row_term_label = "unset";
};

struct SourceMapPressureRowData {
  std::vector<SourceMapRowTerm> lhs_terms;
  std::vector<SourceMapRhsTerm> rhs_terms;
  amrex::Real rhs = std::numeric_limits<amrex::Real>::quiet_NaN();
  int nonfinite_count = 0;
};

std::string source_map_row_term_label(const std::string& base,
                                      const SourceMapCellIndex& mapped)
{
  if (mapped.boundary_label.empty()) {
    return base;
  }
  return base + ";" + mapped.boundary_label;
}

SourceMapPressureRowData source_map_pressure_row_data(const CoupledHostGrid& grid,
                                                      const std::vector<amrex::Real>& p,
                                                      const std::vector<SourceMapMomentum>& lagged_momentum,
                                                      int i, int j, const RunConfig& cfg,
                                                      const amrex::Geometry& geom,
                                                      amrex::Real dt)
{
  const auto dx = geom.CellSizeArray();
  const amrex::Real inv_dx = amrex::Real(1.0) / dx[0];
  const amrex::Real inv_dy = amrex::Real(1.0) / dx[1];
  const amrex::Real gamma_minus_one = cfg.gamma - amrex::Real(1.0);
  const amrex::Real x_coeff = amrex::Real(0.25) * dt * dt * inv_dx * inv_dx;
  const amrex::Real y_coeff = amrex::Real(0.25) * dt * dt * inv_dy * inv_dy;

  SourceMapPressureRowData row;
  const auto add_lhs = [&](int pi, int pj, amrex::Real coeff, const std::string& label) {
    amrex::Real boundary_pressure = std::numeric_limits<amrex::Real>::quiet_NaN();
    if (source_map_fixed_pressure_boundary_value(grid, cfg, pi, pj, boundary_pressure)) {
      const amrex::Real rhs_contribution = -coeff * boundary_pressure;
      row.rhs_terms.push_back(
          SourceMapRhsTerm{rhs_contribution,
                           label + ";boundary_fixed_post_shock_dirichlet_pressure"});
      row.rhs += rhs_contribution;
      if (!std::isfinite(rhs_contribution)) {
        ++row.nonfinite_count;
      }
      return;
    }
    const SourceMapCellIndex mapped = source_map_pressure_index(grid, pi, pj);
    row.lhs_terms.push_back(SourceMapRowTerm{mapped.i, mapped.j, coeff,
                                             source_map_row_term_label(label, mapped)});
    if (!std::isfinite(coeff)) {
      ++row.nonfinite_count;
    }
  };
  const auto add_rhs = [&](amrex::Real value, const std::string& label) {
    row.rhs_terms.push_back(SourceMapRhsTerm{value, label});
    if (!std::isfinite(value)) {
      ++row.nonfinite_count;
    }
    row.rhs += value;
  };
  row.rhs = amrex::Real(0.0);
  const amrex::Real h_e =
      source_map_face_flux_detail(grid, p, lagged_momentum, i, j, i + 1, j, 0, cfg, dt,
                                  inv_dx, inv_dy).h_hat;
  const amrex::Real h_w =
      source_map_face_flux_detail(grid, p, lagged_momentum, i - 1, j, i, j, 0, cfg, dt,
                                  inv_dx, inv_dy).h_hat;
  const amrex::Real h_n =
      source_map_face_flux_detail(grid, p, lagged_momentum, i, j, i, j + 1, 1, cfg, dt,
                                  inv_dx, inv_dy).h_hat;
  const amrex::Real h_s =
      source_map_face_flux_detail(grid, p, lagged_momentum, i, j - 1, i, j, 1, cfg, dt,
                                  inv_dx, inv_dy).h_hat;
  if (!std::isfinite(h_e) || !std::isfinite(h_w) || !std::isfinite(h_n) ||
      !std::isfinite(h_s)) {
    ++row.nonfinite_count;
  }

  add_lhs(i, j, amrex::Real(1.0) / gamma_minus_one, "diag_eos_internal_energy");
  add_lhs(i + 2, j, -x_coeff * h_e, "x_p2_outer_enthalpy");
  add_lhs(i + 1, j, -x_coeff * (h_e - h_w), "x_p1_enthalpy_difference");
  add_lhs(i, j, x_coeff * (h_e + h_w), "x_diag_enthalpy_sum");
  add_lhs(i - 1, j, x_coeff * (h_e - h_w), "x_m1_enthalpy_difference");
  add_lhs(i - 2, j, -x_coeff * h_w, "x_m2_outer_enthalpy");
  add_lhs(i, j + 2, -y_coeff * h_n, "y_p2_outer_enthalpy");
  add_lhs(i, j + 1, -y_coeff * (h_n - h_s), "y_p1_enthalpy_difference");
  add_lhs(i, j, y_coeff * (h_n + h_s), "y_diag_enthalpy_sum");
  add_lhs(i, j - 1, y_coeff * (h_n - h_s), "y_m1_enthalpy_difference");
  add_lhs(i, j - 2, -y_coeff * h_s, "y_m2_outer_enthalpy");

  const ConservedState& star = grid.at(i, j);
  const SourceMapMomentum momentum = source_map_lagged_momentum_at(grid, lagged_momentum, i, j);
  const amrex::Real kinetic =
      amrex::Real(0.5) * (momentum[0] * momentum[0] + momentum[1] * momentum[1]) / star.rho;
  add_rhs(star.E, "rhs_energy_star");
  add_rhs(-kinetic, "rhs_kinetic_lag");
  add_rhs(-amrex::Real(0.5) * dt * inv_dx * h_e *
              source_map_state_at_physical_boundary(grid, cfg, i + 1, j).mx,
          "rhs_x_plus_star_momentum");
  add_rhs(-amrex::Real(0.5) * dt * inv_dx * (h_e - h_w) * grid.at(i, j).mx,
          "rhs_x_center_star_momentum");
  add_rhs(amrex::Real(0.5) * dt * inv_dx * h_w *
              source_map_state_at_physical_boundary(grid, cfg, i - 1, j).mx,
          "rhs_x_minus_star_momentum");
  add_rhs(-amrex::Real(0.5) * dt * inv_dy * h_n *
              source_map_state_at_physical_boundary(grid, cfg, i, j + 1).my,
          "rhs_y_plus_star_momentum");
  add_rhs(-amrex::Real(0.5) * dt * inv_dy * (h_n - h_s) * grid.at(i, j).my,
          "rhs_y_center_star_momentum");
  add_rhs(amrex::Real(0.5) * dt * inv_dy * h_s *
              source_map_state_at_physical_boundary(grid, cfg, i, j - 1).my,
          "rhs_y_minus_star_momentum");
  if (!std::isfinite(row.rhs)) {
    ++row.nonfinite_count;
  }
  return row;
}

using HostSparseMatrix = Eigen::SparseMatrix<amrex::Real, Eigen::RowMajor>;
