// Direct BDLTV20 T1/S2 driver used by formula-level verification rows.
// BDLTV20-family IMEX update implementations and dispatch helpers.
bool bdltv20_paper_t1_s2_requested(const RunConfig& cfg)
{
  return cfg.bdltv20_paper_t1_s2 != "off";
}

amrex::Real compute_bdltv20_paper_t1_s2_dt(const amrex::MultiFab& state,
                                           const amrex::Geometry& geom,
                                           const RunConfig& cfg,
                                           amrex::Real time)
{
  if (cfg.bdltv20_paper_t1_s2_dt > amrex::Real(0.0)) {
    return std::min(cfg.bdltv20_paper_t1_s2_dt, cfg.stop_time - time);
  }
  const auto dx = geom.CellSizeArray();
  amrex::Real material_rate = amrex::Real(0.0);
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
      for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
        const PrimitiveState q =
            bdltv20_paper_to_primitive(load_cons(arr, i, j, 0), cfg.gamma,
                                       cfg.bdltv20_paper_epsilon);
        if (std::isfinite(q.u) && std::isfinite(q.v)) {
          const amrex::Real speed =
              std::sqrt(q.u * q.u + q.v * q.v);
          const amrex::Real min_dx = std::min(dx[0], dx[1]);
          material_rate =
              std::max(material_rate, amrex::Real(2.0) * speed / min_dx);
        }
      }
    }
  }
  amrex::ParallelDescriptor::ReduceRealMax(material_rate);
  const amrex::Real remaining = cfg.stop_time - time;
  if (!(material_rate > amrex::Real(0.0))) {
    return remaining;
  }
  return std::min(cfg.imex_cfl / material_rate, remaining);
}

int run_bdltv20_paper_t1_s2(const RunConfig& cfg)
{
  const amrex::Geometry geom = make_geometry(cfg);
  const amrex::Box domain = geom.Domain();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const int nx = cfg.n_cell[0];
  const int ny = cfg.n_cell[1];
  const bool periodic_x = geom.isPeriodic(0);
  const bool periodic_y = geom.isPeriodic(1);
  const int cell_count = nx * ny;
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real gm1 = gamma - amrex::Real(1.0);
  const amrex::Real eps = cfg.bdltv20_paper_epsilon;
  const bool is_toro =
      cfg.bdltv20_paper_t1_s2 == "toro_x_exact_dirichlet_y_periodic_2d" ||
      cfg.bdltv20_paper_t1_s2 == "toro_xy_exact_dirichlet_2d";
  const bool is_gresho =
      cfg.bdltv20_paper_t1_s2 == "gresho_exact_dirichlet_2d";
  const bool is_isentropic_vortex =
      cfg.bdltv20_paper_t1_s2 == "isentropic_vortex_periodic_2d";
  const bool is_advection_blob =
      cfg.bdltv20_paper_t1_s2 == "advection_blob_periodic_2d";
  const bool use_riemann_x_only_pressure_row = is_toro;
  amrex::BoxArray ba(domain);
  ba.maxSize(cfg.max_grid_size);
  const amrex::DistributionMapping dm(ba);
  amrex::MultiFab state(ba, dm, NCons, 2);
  amrex::MultiFab star(ba, dm, NCons, 2);
  initialize_state(state, geom, cfg);
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.array(mfi);
    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
      for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
        const PrimitiveState q = to_primitive(load_cons(arr, i, j, 0), gamma);
        store_cons(arr, i, j, 0, bdltv20_paper_to_conserved(q, gamma, eps));
      }
    }
  }

  const ToroState ts = toro_state(cfg.toro_test);
  const PrimitiveState toro_left{ts.rho_l, ts.u_l, ts.v_l, ts.p_l};
  const PrimitiveState toro_right{ts.rho_r, ts.u_r, ts.v_r, ts.p_r};
  const auto row = [=](int i, int j) {
    return pressure_cell_row(i, j, nx);
  };
  const auto wrap = [](int index, int n) {
    return pressure_periodic_index(index, n);
  };
  const auto x_center = [=](int i) {
    return plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
  };
  const auto y_center = [=](int j) {
    return plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
  };
  const auto exact_boundary_state = [&](int i, int j, amrex::Real time_sample) {
    amrex::Real x = x_center(std::min(std::max(i, 0), nx - 1));
    amrex::Real y = y_center(std::min(std::max(j, 0), ny - 1));
    if (i < 0) {
      x = plo[0];
    } else if (i >= nx) {
      x = phi[0];
    }
    if (j < 0) {
      y = plo[1];
    } else if (j >= ny) {
      y = phi[1];
    }
    if (is_toro) {
      return bdltv20_paper_to_conserved(
          exact_riemann_sample(x, time_sample, toro_left, toro_right, gamma,
                               cfg.riemann_interface_x),
          gamma, eps);
    }
    if (is_isentropic_vortex) {
      return bdltv20_paper_to_conserved(
          isentropic_vortex_state(x, y, cfg, plo[0], phi[0], plo[1], phi[1],
                                  time_sample),
          gamma, eps);
    }
    if (is_advection_blob) {
      return bdltv20_paper_to_conserved(
          exact_advection_blob_state(x, y, geom, cfg, time_sample), gamma, eps);
    }
    return bdltv20_paper_to_conserved(gresho_state(x, y, cfg), gamma, eps);
  };
  const auto exact_boundary_pressure = [&](int i, int j, amrex::Real time_sample) {
    return bdltv20_paper_to_primitive(exact_boundary_state(i, j, time_sample),
                                      gamma, eps)
        .p;
  };

  int steps = 0;
  amrex::Real time = amrex::Real(0.0);
  int solver_failure_count = 0;
  int nonfinite_count = 0;
  int nonpositive_count = 0;
  int enthalpy_guard_count = 0;
  int max_gmres_iterations = 0;
  amrex::Real pressure_residual_linf_max = amrex::Real(0.0);
  amrex::Real pressure_relative_residual_linf_max = amrex::Real(0.0);
  StaggeredPressureMatrixSymbolDiagnostics last_matrix_diag;

  while (time < cfg.stop_time - amrex::Real(1.0e-14) &&
         steps < cfg.bdltv20_paper_t1_s2_max_steps) {
    const amrex::Real dt = compute_bdltv20_paper_t1_s2_dt(state, geom, cfg, time);
    if (!(dt > amrex::Real(0.0)) || !std::isfinite(dt)) {
      ++nonfinite_count;
      break;
    }

    std::vector<amrex::Real> old_pressure(static_cast<std::size_t>(cell_count),
                                          amrex::Real(0.0));
    std::vector<ConservedState> old_cells(static_cast<std::size_t>(cell_count));
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      auto const arr = state.const_array(mfi);
      for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
        for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
          const int r = row(i, j);
          const ConservedState old_u = load_cons(arr, i, j, 0);
          old_cells[static_cast<std::size_t>(r)] = old_u;
          old_pressure[static_cast<std::size_t>(r)] =
              bdltv20_paper_to_primitive(old_u, gamma, eps).p;
        }
      }
    }

    std::vector<ConservedState> star_cells(static_cast<std::size_t>(cell_count));
    std::vector<amrex::Real> star_pressure(static_cast<std::size_t>(cell_count),
                                           amrex::Real(0.0));
    const auto old_cell_from = [&](int i, int j, amrex::Real time_sample) {
      if (periodic_x) {
        i = wrap(i, nx);
      }
      if (periodic_y) {
        j = wrap(j, ny);
      }
      if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return old_cells[static_cast<std::size_t>(row(i, j))];
      }
      return exact_boundary_state(i, j, time_sample);
    };
    const auto reconstructed_old = [&](int i, int j, int axis, int side) {
      const ConservedState u_m =
          axis == 0 ? old_cell_from(i - 1, j, time)
                    : old_cell_from(i, j - 1, time);
      const ConservedState u_c = old_cell_from(i, j, time);
      const ConservedState u_p =
          axis == 0 ? old_cell_from(i + 1, j, time)
                    : old_cell_from(i, j + 1, time);
      const ConservedState backward = sub_cons(u_c, u_m);
      const ConservedState forward = sub_cons(u_p, u_c);
      const ConservedState slope{minmod(backward.rho, forward.rho),
                                 minmod(backward.mx, forward.mx),
                                 minmod(backward.my, forward.my),
                                 minmod(backward.E, forward.E)};
      return add_cons(u_c, scale_cons(amrex::Real(0.5 * side), slope));
    };
    const auto explicit_flux = [&](const ConservedState& left_u,
                                   const ConservedState& right_u, int axis) {
      const PrimitiveState left = bdltv20_paper_to_primitive(left_u, gamma, eps);
      const PrimitiveState right = bdltv20_paper_to_primitive(right_u, gamma, eps);
      const amrex::Real mx_l = left.rho * left.u;
      const amrex::Real my_l = left.rho * left.v;
      const amrex::Real mx_r = right.rho * right.u;
      const amrex::Real my_r = right.rho * right.v;
      const amrex::Real k_l =
          eps * amrex::Real(0.5) * left.rho * (left.u * left.u + left.v * left.v);
      const amrex::Real k_r =
          eps * amrex::Real(0.5) * right.rho * (right.u * right.u + right.v * right.v);
      ConservedState f_l;
      ConservedState f_r;
      if (axis == 0) {
        f_l = ConservedState{mx_l, mx_l * left.u, my_l * left.u, k_l * left.u};
        f_r = ConservedState{mx_r, mx_r * right.u, my_r * right.u, k_r * right.u};
      } else {
        f_l = ConservedState{my_l, mx_l * left.v, my_l * left.v, k_l * left.v};
        f_r = ConservedState{my_r, mx_r * right.v, my_r * right.v, k_r * right.v};
      }
      const amrex::Real smax =
          amrex::max(std::abs(normal_velocity(left, axis)),
                     std::abs(normal_velocity(right, axis)));
      ConservedState flux =
          sub_cons(scale_cons(amrex::Real(0.5), add_cons(f_l, f_r)),
                   scale_cons(amrex::Real(0.5) * smax, sub_cons(right_u, left_u)));
	      flux.E = amrex::Real(0.5) * (f_l.E + f_r.E) -
	               amrex::Real(0.5) * smax * (right_u.E - left_u.E);
      return flux;
    };
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        const ConservedState fx_l =
            explicit_flux(reconstructed_old(i - 1, j, 0, +1),
                          reconstructed_old(i, j, 0, -1), 0);
        const ConservedState fx_r =
            explicit_flux(reconstructed_old(i, j, 0, +1),
                          reconstructed_old(i + 1, j, 0, -1), 0);
        const ConservedState fy_b =
            explicit_flux(reconstructed_old(i, j - 1, 1, +1),
                          reconstructed_old(i, j, 1, -1), 1);
        const ConservedState fy_t =
            explicit_flux(reconstructed_old(i, j, 1, +1),
                          reconstructed_old(i, j + 1, 1, -1), 1);
        const ConservedState old_u = old_cells[static_cast<std::size_t>(row(i, j))];
        ConservedState updated;
        updated.rho = old_u.rho - dt * ((fx_r.rho - fx_l.rho) / dx[0] +
                                        (fy_t.rho - fy_b.rho) / dx[1]);
        updated.mx = old_u.mx - dt * ((fx_r.mx - fx_l.mx) / dx[0] +
                                      (fy_t.mx - fy_b.mx) / dx[1]);
        updated.my = old_u.my - dt * ((fx_r.my - fx_l.my) / dx[0] +
                                      (fy_t.my - fy_b.my) / dx[1]);
        updated.E = old_u.E - dt * ((fx_r.E - fx_l.E) / dx[0] +
                                    (fy_t.E - fy_b.E) / dx[1]);
        const int r = row(i, j);
        star_cells[static_cast<std::size_t>(r)] = updated;
        star_pressure[static_cast<std::size_t>(r)] =
            bdltv20_paper_to_primitive(updated, gamma, eps).p;
      }
    }

    const int x_face_count = pressure_x_face_count(nx, ny, periodic_x);
    const int y_face_count = pressure_y_face_count(nx, ny, periodic_y);
    std::vector<amrex::Real> h_x(static_cast<std::size_t>(x_face_count), amrex::Real(0.0));
    std::vector<amrex::Real> h_y(static_cast<std::size_t>(y_face_count), amrex::Real(0.0));
    std::vector<amrex::Real> m_x(static_cast<std::size_t>(x_face_count), amrex::Real(0.0));
    std::vector<amrex::Real> m_y(static_cast<std::size_t>(y_face_count), amrex::Real(0.0));
    std::vector<ConservedState> corrected_cells(static_cast<std::size_t>(cell_count));

    const auto cell_state_from = [&](int i, int j,
                                     const std::vector<ConservedState>& cells,
                                     amrex::Real time_sample) {
      if (periodic_x) {
        i = wrap(i, nx);
      }
      if (periodic_y) {
        j = wrap(j, ny);
      }
      if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return cells[static_cast<std::size_t>(row(i, j))];
      }
      return exact_boundary_state(i, j, time_sample);
    };
    const auto pressure_from = [&](int i, int j, const Eigen::VectorXd* pressure,
                                   amrex::Real time_sample) {
      if (periodic_x) {
        i = wrap(i, nx);
      }
      if (periodic_y) {
        j = wrap(j, ny);
      }
	      if (i >= 0 && i < nx && j >= 0 && j < ny) {
	        const int r = row(i, j);
	        return pressure == nullptr ? old_pressure[static_cast<std::size_t>(r)]
	                                   : (*pressure)(r);
      }
      return exact_boundary_pressure(i, j, time_sample);
    };
    const auto old_state_from = [&](int i, int j, amrex::Real time_sample) {
      if (periodic_x) {
        i = wrap(i, nx);
      }
      if (periodic_y) {
        j = wrap(j, ny);
      }
      if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return old_cells[static_cast<std::size_t>(row(i, j))];
      }
      return exact_boundary_state(i, j, time_sample);
    };
    const auto refresh_corrected_cells_and_faces =
        [&](const Eigen::VectorXd* pressure, amrex::Real time_sample) {
          int local_nonfinite = 0;
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              const int r = row(i, j);
              ConservedState u = star_cells[static_cast<std::size_t>(r)];
              if (pressure != nullptr) {
                const amrex::Real p_e = pressure_from(i + 1, j, pressure, time_sample);
                const amrex::Real p_w = pressure_from(i - 1, j, pressure, time_sample);
                const amrex::Real p_n = pressure_from(i, j + 1, pressure, time_sample);
                const amrex::Real p_s = pressure_from(i, j - 1, pressure, time_sample);
                u.mx -= dt * (p_e - p_w) /
                        (amrex::Real(2.0) * eps * dx[0]);
	                if (!use_riemann_x_only_pressure_row) {
	                  u.my -= dt * (p_n - p_s) /
	                          (amrex::Real(2.0) * eps * dx[1]);
	                }
              }
              corrected_cells[static_cast<std::size_t>(r)] = u;
            }
          }
          const int x_faces_per_row = periodic_x ? nx : nx + 1;
          for (int j = 0; j < ny; ++j) {
            for (int face_i = 0; face_i < x_faces_per_row; ++face_i) {
              const int left_i = periodic_x ? face_i : face_i - 1;
              const int right_i = periodic_x ? face_i + 1 : face_i;
              const ConservedState left =
                  cell_state_from(left_i, j, corrected_cells, time_sample);
              const ConservedState right =
                  cell_state_from(right_i, j, corrected_cells, time_sample);
              const amrex::Real p_l = pressure_from(left_i, j, pressure, time_sample);
              const amrex::Real p_r = pressure_from(right_i, j, pressure, time_sample);
              const int face = pressure_x_face_index(face_i, j, nx, periodic_x);
              const auto idx = static_cast<std::size_t>(face);
              h_x[idx] = bdltv20_eq46_face_enthalpy(left.rho, p_l, left.mx,
                                                    right.rho, p_r, right.mx,
                                                    gamma, enthalpy_guard_count);
              m_x[idx] = amrex::Real(0.5) * (left.mx + right.mx);
              if (!std::isfinite(h_x[idx]) || !std::isfinite(m_x[idx])) {
                ++local_nonfinite;
              }
            }
          }
          const int y_faces_per_col = periodic_y ? ny : ny + 1;
          for (int face_j = 0; face_j < y_faces_per_col; ++face_j) {
            for (int i = 0; i < nx; ++i) {
              const int bottom_j = periodic_y ? face_j : face_j - 1;
              const int top_j = periodic_y ? face_j + 1 : face_j;
              const ConservedState bottom =
                  cell_state_from(i, bottom_j, corrected_cells, time_sample);
              const ConservedState top =
                  cell_state_from(i, top_j, corrected_cells, time_sample);
              const amrex::Real p_b = pressure_from(i, bottom_j, pressure, time_sample);
              const amrex::Real p_t = pressure_from(i, top_j, pressure, time_sample);
              const int face = pressure_y_face_index(i, face_j, nx, ny, periodic_y);
              const auto idx = static_cast<std::size_t>(face);
              h_y[idx] = bdltv20_eq46_face_enthalpy(bottom.rho, p_b, bottom.my,
                                                    top.rho, p_t, top.my,
                                                    gamma, enthalpy_guard_count);
              m_y[idx] = amrex::Real(0.5) * (bottom.my + top.my);
              if (!std::isfinite(h_y[idx]) || !std::isfinite(m_y[idx])) {
                ++local_nonfinite;
              }
            }
          }
          return local_nonfinite;
        };

    Eigen::VectorXd pressure = Eigen::VectorXd::Zero(cell_count);
	    for (int r = 0; r < cell_count; ++r) {
	      pressure(r) = old_pressure[static_cast<std::size_t>(r)];
    }

    HostGmresResult solver;
    for (int picard = 0; picard < cfg.imex_picard_iterations; ++picard) {
      nonfinite_count += refresh_corrected_cells_and_faces(
          picard == 0 ? nullptr : &pressure, time + dt);
      last_matrix_diag = StaggeredPressureMatrixSymbolDiagnostics{};
      last_matrix_diag.rows = cell_count;
      last_matrix_diag.cols = cell_count;
      last_matrix_diag.eos_diag = eps / gm1;
	      last_matrix_diag.pressure_topology =
	          use_riemann_x_only_pressure_row
	              ? "paper_eq49_50_1d_x_riemann_row_embedded_in_2d_strip"
	              : "paper_collocated_source_map_eq76_77_2d_zdrop";
	      if (use_riemann_x_only_pressure_row) {
	        last_matrix_diag.pressure_matrix_construction =
	            "bdltv20_eq49_50_1d_x_riemann_row_embedded_in_2d_strip";
	      } else if (!periodic_x && periodic_y) {
	        last_matrix_diag.pressure_matrix_construction =
	            "bdltv20_eq76_77_x_exact_dirichlet_y_periodic_2d_zdrop";
      } else if (!periodic_x && !periodic_y) {
        last_matrix_diag.pressure_matrix_construction =
            "bdltv20_eq76_77_xy_exact_dirichlet_2d_zdrop";
      } else {
        last_matrix_diag.pressure_matrix_construction =
            "bdltv20_eq76_77_periodic_or_mixed_2d_zdrop";
      }
	      last_matrix_diag.pressure_matrix_stencil =
	          use_riemann_x_only_pressure_row
	              ? "wide_second_neighbor_collocated_eq49_1d_x_embedded"
	              : "wide_second_neighbor_collocated_eq76_2d_zdrop";
      last_matrix_diag.nearest_neighbor_control =
          "bdltv20_mean_centered_second_neighbor_pressure_row";
      last_matrix_diag.uniform_h_checkerboard_status =
          "paper_eq46_mean_centered_enthalpy_control_no_ap_proof";

      std::vector<Eigen::Triplet<amrex::Real>> triplets;
      triplets.reserve(static_cast<std::size_t>(cell_count * 9));
      Eigen::VectorXd rhs = Eigen::VectorXd::Zero(cell_count);
      const amrex::Real cx = amrex::Real(0.25) * dt * dt / (dx[0] * dx[0]);
	      const amrex::Real cy =
	          use_riemann_x_only_pressure_row
	              ? amrex::Real(0.0)
	              : amrex::Real(0.25) * dt * dt / (dx[1] * dx[1]);
      const auto inside = [&](int ci, int cj) {
        if (periodic_x) {
          ci = wrap(ci, nx);
        }
        if (periodic_y) {
          cj = wrap(cj, ny);
        }
        return ci >= 0 && ci < nx && cj >= 0 && cj < ny;
      };
      const auto mapped_row = [&](int ci, int cj) {
        if (periodic_x) {
          ci = wrap(ci, nx);
        }
        if (periodic_y) {
          cj = wrap(cj, ny);
        }
        return row(ci, cj);
      };
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          const int r = row(i, j);
          const ConservedState star_u = star_cells[static_cast<std::size_t>(r)];
          const ConservedState lag_u = corrected_cells[static_cast<std::size_t>(r)];
          const int e = pressure_x_face_east(i, j, nx, periodic_x);
          const int w = pressure_x_face_west(i, j, nx, periodic_x);
          const int n = pressure_y_face_north(i, j, nx, ny, periodic_y);
          const int s = pressure_y_face_south(i, j, nx, ny, periodic_y);
          const auto east = static_cast<std::size_t>(e);
          const auto west = static_cast<std::size_t>(w);
          const auto north = static_cast<std::size_t>(n);
          const auto south = static_cast<std::size_t>(s);
          const amrex::Real kinetic =
              eps * amrex::Real(0.5) *
              (lag_u.mx * lag_u.mx + lag_u.my * lag_u.my) /
              star_u.rho;
          const ConservedState star_e = cell_state_from(i + 1, j, star_cells, time + dt);
          const ConservedState star_w = cell_state_from(i - 1, j, star_cells, time + dt);
          const ConservedState star_n = cell_state_from(i, j + 1, star_cells, time + dt);
          const ConservedState star_s = cell_state_from(i, j - 1, star_cells, time + dt);
          const amrex::Real enthalpy_flux_div =
              (h_x[east] * star_e.mx +
               (h_x[east] - h_x[west]) * star_u.mx -
               h_x[west] * star_w.mx) /
                  (amrex::Real(2.0) * dx[0]) +
              (h_y[north] * star_n.my +
               (h_y[north] - h_y[south]) * star_u.my -
               h_y[south] * star_s.my) /
                  (amrex::Real(2.0) * dx[1]);
          amrex::Real rhs_value =
              eps * (star_u.E - kinetic - dt * enthalpy_flux_div);
          const auto add_lhs = [&](int ci, int cj, amrex::Real coeff) {
            if (!std::isfinite(coeff)) {
              ++last_matrix_diag.nonfinite_count;
            }
            if (coeff == amrex::Real(0.0)) {
              return;
            }
            if (inside(ci, cj)) {
              triplets.emplace_back(r, mapped_row(ci, cj), coeff);
            } else {
              rhs_value -= coeff * exact_boundary_pressure(ci, cj, time + dt);
            }
          };
          add_lhs(i + 2, j, -cx * h_x[east]);
          add_lhs(i + 1, j, -cx * (h_x[east] - h_x[west]));
          add_lhs(i, j, last_matrix_diag.eos_diag +
                            cx * (h_x[east] + h_x[west]) +
                            cy * (h_y[north] + h_y[south]));
          add_lhs(i - 1, j, cx * (h_x[east] - h_x[west]));
          add_lhs(i - 2, j, -cx * h_x[west]);
          add_lhs(i, j + 2, -cy * h_y[north]);
          add_lhs(i, j + 1, -cy * (h_y[north] - h_y[south]));
          add_lhs(i, j - 1, cy * (h_y[north] - h_y[south]));
          add_lhs(i, j - 2, -cy * h_y[south]);
          rhs(r) = rhs_value;
          if (!std::isfinite(rhs_value) || !(star_u.rho > amrex::Real(0.0))) {
            ++nonfinite_count;
          }
        }
      }
      HostSparseMatrix pressure_matrix(cell_count, cell_count);
      pressure_matrix.setFromTriplets(triplets.begin(), triplets.end());
      pressure_matrix.makeCompressed();
      last_matrix_diag.nonzeros = static_cast<int>(pressure_matrix.nonZeros());
      std::vector<amrex::Real> row_sum(static_cast<std::size_t>(cell_count),
                                       amrex::Real(0.0));
      std::vector<amrex::Real> row_diag(static_cast<std::size_t>(cell_count),
                                        amrex::Real(0.0));
      std::vector<amrex::Real> row_abs_offdiag(static_cast<std::size_t>(cell_count),
                                               amrex::Real(0.0));
      for (int matrix_row = 0; matrix_row < pressure_matrix.outerSize(); ++matrix_row) {
        for (HostSparseMatrix::InnerIterator it(pressure_matrix, matrix_row); it; ++it) {
          const amrex::Real value = it.value();
          if (!std::isfinite(value)) {
            ++last_matrix_diag.nonfinite_count;
          }
          row_sum[static_cast<std::size_t>(it.row())] += value;
          if (it.row() == it.col()) {
            row_diag[static_cast<std::size_t>(it.row())] += value;
          } else {
            row_abs_offdiag[static_cast<std::size_t>(it.row())] += std::abs(value);
            if (value > amrex::Real(0.0)) {
              ++last_matrix_diag.offdiag_positive_count;
            }
          }
        }
      }
      for (int matrix_row = 0; matrix_row < cell_count; ++matrix_row) {
        const auto idx = static_cast<std::size_t>(matrix_row);
        last_matrix_diag.diag_min = std::min(last_matrix_diag.diag_min, row_diag[idx]);
        last_matrix_diag.diag_max = std::max(last_matrix_diag.diag_max, row_diag[idx]);
        last_matrix_diag.row_sum_linf =
            std::max(last_matrix_diag.row_sum_linf, std::abs(row_sum[idx]));
        const amrex::Real dominance = row_diag[idx] - row_abs_offdiag[idx];
        last_matrix_diag.diag_dominance_min =
            std::min(last_matrix_diag.diag_dominance_min, dominance);
        last_matrix_diag.gershgorin_lower_bound =
            std::min(last_matrix_diag.gershgorin_lower_bound, dominance);
      }

      solver = solve_host_sparse_gmres(pressure_matrix, rhs, pressure,
                                       cfg.imex_solver_max_iter,
                                       cfg.imex_solver_tol);
      if (!solver.converged || solver.nonfinite || solver.breakdown) {
        ++solver_failure_count;
      }
      max_gmres_iterations = std::max(max_gmres_iterations, solver.iterations);
      pressure_residual_linf_max =
          std::max(pressure_residual_linf_max, solver.residual_linf);
      pressure_relative_residual_linf_max =
          std::max(pressure_relative_residual_linf_max, solver.relative_residual_linf);
      pressure = solver.solution;
    }

    nonfinite_count += refresh_corrected_cells_and_faces(&pressure, time + dt);
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      auto const dst = state.array(mfi);
      for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
        for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
          const int r = row(i, j);
          ConservedState u = corrected_cells[static_cast<std::size_t>(r)];
          const int e = pressure_x_face_east(i, j, nx, periodic_x);
          const int w = pressure_x_face_west(i, j, nx, periodic_x);
          const int n = pressure_y_face_north(i, j, nx, ny, periodic_y);
          const int s = pressure_y_face_south(i, j, nx, ny, periodic_y);
          const auto east = static_cast<std::size_t>(e);
          const auto west = static_cast<std::size_t>(w);
          const auto north = static_cast<std::size_t>(n);
          const auto south = static_cast<std::size_t>(s);
	          const amrex::Real kinetic =
	              eps * amrex::Real(0.5) * (u.mx * u.mx + u.my * u.my) / u.rho;
	          u.E = pressure(r) / gm1 + kinetic;
          const PrimitiveState q = bdltv20_paper_to_primitive(u, gamma, eps);
          if (!std::isfinite(q.rho) || !std::isfinite(q.u) ||
              !std::isfinite(q.v) || !std::isfinite(q.p)) {
            ++nonfinite_count;
          }
          if (!(q.rho > amrex::Real(0.0)) || !(q.p > amrex::Real(0.0))) {
            ++nonpositive_count;
          }
          store_cons(dst, i, j, 0, u);
        }
      }
    }
    time += dt;
    ++steps;
    fill_problem_ghosts(state, geom, cfg, time);
  }

  amrex::Real density_l1 = amrex::Real(0.0);
  amrex::Real velocity_l1 = amrex::Real(0.0);
  amrex::Real pressure_l1 = amrex::Real(0.0);
  amrex::Real density_l2_sum = amrex::Real(0.0);
  amrex::Real velocity_l2_sum = amrex::Real(0.0);
  amrex::Real pressure_l2_sum = amrex::Real(0.0);
  amrex::Real density_linf = amrex::Real(0.0);
  amrex::Real pressure_linf = amrex::Real(0.0);
  amrex::Real rho_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real pressure_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Long local_cells = 0;
  for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto const arr = state.const_array(mfi);
    for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
      for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
        const amrex::Real x = x_center(i);
        const amrex::Real y = y_center(j);
        const PrimitiveState q =
            bdltv20_paper_to_primitive(load_cons(arr, i, j, 0), gamma, eps);
        PrimitiveState exact;
        if (is_toro) {
          exact = exact_riemann_sample(x, time, toro_left, toro_right, gamma,
                                       cfg.riemann_interface_x);
        } else if (is_isentropic_vortex) {
          exact = isentropic_vortex_state(x, y, cfg, plo[0], phi[0], plo[1],
                                          phi[1], time);
        } else if (is_advection_blob) {
          exact = exact_advection_blob_state(x, y, geom, cfg, time);
        } else {
          exact = gresho_state(x, y, cfg);
        }
        const amrex::Real du =
            std::sqrt((q.u - exact.u) * (q.u - exact.u) +
                      (q.v - exact.v) * (q.v - exact.v));
        const amrex::Real drho = std::abs(q.rho - exact.rho);
        const amrex::Real dp = std::abs(q.p - exact.p);
        density_l1 += drho;
        velocity_l1 += du;
        pressure_l1 += dp;
        density_l2_sum += drho * drho;
        velocity_l2_sum += du * du;
        pressure_l2_sum += dp * dp;
        density_linf = std::max(density_linf, drho);
        pressure_linf = std::max(pressure_linf, dp);
        rho_min = std::min(rho_min, q.rho);
        pressure_min = std::min(pressure_min, q.p);
        ++local_cells;
      }
    }
  }
  amrex::ParallelDescriptor::ReduceRealSum(density_l1);
  amrex::ParallelDescriptor::ReduceRealSum(velocity_l1);
  amrex::ParallelDescriptor::ReduceRealSum(pressure_l1);
  amrex::ParallelDescriptor::ReduceRealSum(density_l2_sum);
  amrex::ParallelDescriptor::ReduceRealSum(velocity_l2_sum);
  amrex::ParallelDescriptor::ReduceRealSum(pressure_l2_sum);
  amrex::ParallelDescriptor::ReduceRealMax(density_linf);
  amrex::ParallelDescriptor::ReduceRealMax(pressure_linf);
  amrex::ParallelDescriptor::ReduceRealMin(rho_min);
  amrex::ParallelDescriptor::ReduceRealMin(pressure_min);
  amrex::ParallelDescriptor::ReduceLongSum(local_cells);
  const amrex::Real denom = static_cast<amrex::Real>(std::max<amrex::Long>(local_cells, 1));
  density_l1 /= denom;
  velocity_l1 /= denom;
  pressure_l1 /= denom;
  const amrex::Real density_l2 = std::sqrt(density_l2_sum / denom);
  const amrex::Real velocity_l2 = std::sqrt(velocity_l2_sum / denom);
  const amrex::Real pressure_l2 = std::sqrt(pressure_l2_sum / denom);

  bool final_csv_written = false;
	  if (!cfg.final_csv.empty() && amrex::ParallelDescriptor::IOProcessor()) {
	    std::ofstream out(cfg.final_csv);
	    out << std::setprecision(17);
	    out << "x,y,rho,u,v,p,exact_rho,exact_u,exact_v,exact_p\n";
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      auto const arr = state.const_array(mfi);
      for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
        for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
          const amrex::Real x = x_center(i);
          const amrex::Real y = y_center(j);
          const PrimitiveState q =
              bdltv20_paper_to_primitive(load_cons(arr, i, j, 0), gamma, eps);
          const PrimitiveState exact =
              is_toro ? exact_riemann_sample(x, time, toro_left, toro_right, gamma,
                                             cfg.riemann_interface_x)
                      : (is_isentropic_vortex
                             ? isentropic_vortex_state(x, y, cfg, plo[0], phi[0],
                                                       plo[1], phi[1], time)
                             : (is_advection_blob
                                    ? exact_advection_blob_state(x, y, geom, cfg, time)
                                    : gresho_state(x, y, cfg)));
          out << x << ',' << y << ',' << q.rho << ',' << q.u << ',' << q.v << ','
              << q.p << ',' << exact.rho << ',' << exact.u << ',' << exact.v
              << ',' << exact.p << '\n';
        }
      }
    }
    final_csv_written = static_cast<bool>(out);
  }

  const bool finished_time = time >= cfg.stop_time - amrex::Real(1.0e-12);
  const bool passed = finished_time && solver_failure_count == 0 &&
                      nonfinite_count == 0 && nonpositive_count == 0;
  amrex::Print() << std::setprecision(17);
  amrex::Print() << "bdltv20_paper_t1_s2=" << cfg.bdltv20_paper_t1_s2 << "\n";
  amrex::Print() << "bdltv20_paper_primary_source="
                 << "Boscheri_Dimarco_Loubere_Tavelli_Vignal_2020\n";
  amrex::Print() << "bdltv20_paper_dimension=2d_z_dropped\n";
  amrex::Print() << "bdltv20_paper_time_order=1\n";
  amrex::Print() << "bdltv20_paper_space_order=2\n";
  amrex::Print() << "bdltv20_paper_epsilon=" << eps << "\n";
  amrex::Print() << "bdltv20_paper_reconstruction="
                 << "section_4_2_eq87_89_conserved_minmod\n";
	  amrex::Print() << "bdltv20_paper_explicit_update="
	                 << "toro_vazquez_material_flux_lax_friedrichs_star_update_eq69_total_energy_lf_jump\n";
	  amrex::Print() << "bdltv20_paper_face_enthalpy="
	                 << "eq46_momentum_weighted_with_arithmetic_zero_denominator_guard\n";
	  amrex::Print() << "bdltv20_paper_energy_closeout="
	                 << "eq75_pressure_internal_energy_plus_corrected_kinetic_from_solved_pressure\n";
  amrex::Print() << "bdltv20_paper_pressure_solver="
                 << cfg.bdltv20_paper_pressure_solver << "\n";
	  amrex::Print() << "bdltv20_paper_picard_iterations=" << cfg.imex_picard_iterations << "\n";
	  amrex::Print() << "bdltv20_paper_picard_initialization="
	                 << "eq53_star_momentum_old_pressure_initial_guess\n";
	  amrex::Print() << "bdltv20_paper_pressure_topology="
	                 << (use_riemann_x_only_pressure_row
	                         ? "paper_eq49_50_1d_x_riemann_row_embedded_in_2d_strip"
	                         : "paper_collocated_source_map_eq76_77_2d_zdrop") << "\n";
  amrex::Print() << "bdltv20_paper_control_role="
                 << "bdltv20_t1s2_report_path\n";
  amrex::Print() << "bdltv20_paper_boundary_policy="
                 << (is_toro ? (periodic_y
                                    ? "x_exact_dirichlet_y_periodic_paper_2d_reduction"
                                    : "xy_exact_dirichlet_paper_2d_reduction")
                             : (is_isentropic_vortex
                                    ? "xy_periodic_isentropic_vortex_paper_section_5_1"
                                    : (is_advection_blob
                                           ? "xy_periodic_advection_blob_project_project"
                                           : "xy_exact_dirichlet_gresho"))) << "\n";
  amrex::Print() << "bdltv20_paper_timestep_policy="
                 << (cfg.bdltv20_paper_t1_s2_dt > amrex::Real(0.0)
                         ? "fixed_user_dt"
                         : "eq79_min_dx_over_two_velocity_magnitude_zdrop_2d") << "\n";
  amrex::Print() << "bdltv20_paper_steps=" << steps << "\n";
  amrex::Print() << "bdltv20_paper_final_time=" << time << "\n";
  amrex::Print() << "bdltv20_paper_solver_failure_count=" << solver_failure_count << "\n";
  amrex::Print() << "bdltv20_paper_gmres_iterations_max=" << max_gmres_iterations << "\n";
  amrex::Print() << "bdltv20_paper_pressure_residual_linf_max="
                 << pressure_residual_linf_max << "\n";
  amrex::Print() << "bdltv20_paper_pressure_relative_residual_linf_max="
                 << pressure_relative_residual_linf_max << "\n";
  amrex::Print() << "bdltv20_paper_matrix_stencil="
                 << last_matrix_diag.pressure_matrix_stencil << "\n";
  amrex::Print() << "bdltv20_paper_enthalpy_guard_count=" << enthalpy_guard_count << "\n";
  amrex::Print() << "bdltv20_paper_nonfinite_count=" << nonfinite_count << "\n";
  amrex::Print() << "bdltv20_paper_nonpositive_count=" << nonpositive_count << "\n";
  amrex::Print() << "bdltv20_paper_rho_min=" << rho_min << "\n";
  amrex::Print() << "bdltv20_paper_pressure_min=" << pressure_min << "\n";
  amrex::Print() << "bdltv20_paper_density_l1_error=" << density_l1 << "\n";
  amrex::Print() << "bdltv20_paper_velocity_l1_error=" << velocity_l1 << "\n";
  amrex::Print() << "bdltv20_paper_pressure_l1_error=" << pressure_l1 << "\n";
  amrex::Print() << "bdltv20_paper_density_l2_error=" << density_l2 << "\n";
  amrex::Print() << "bdltv20_paper_velocity_l2_error=" << velocity_l2 << "\n";
  amrex::Print() << "bdltv20_paper_pressure_l2_error=" << pressure_l2 << "\n";
  amrex::Print() << "bdltv20_paper_density_linf_error=" << density_linf << "\n";
  amrex::Print() << "bdltv20_paper_pressure_linf_error=" << pressure_linf << "\n";
  amrex::Print() << "bdltv20_paper_final_csv="
                 << (cfg.final_csv.empty() ? "none" : cfg.final_csv) << "\n";
  amrex::Print() << "bdltv20_paper_final_csv_written="
                 << (final_csv_written ? 1 : 0) << "\n";
  amrex::Print() << "bdltv20_paper_nonclaims="
                 << "no_full_3d_no_second_order_time_no_ap_proof_no_production_no_gpu_cpu_efficiency\n";
  amrex::Print() << "bdltv20_paper_t1_s2_status="
                 << (passed ? "passed" : "failed") << "\n";
  return passed ? 0 : 2;
}
