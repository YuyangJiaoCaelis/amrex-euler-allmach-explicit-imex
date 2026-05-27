// Host sparse pressure-solve utilities used by the BDLTV20 T1/S2 IMEX paths.
struct StaggeredPressureMatrixSymbolDiagnostics {
  int rows = 0;
  int cols = 0;
  int nonzeros = 0;
  int offdiag_positive_count = 0;
  int nonfinite_count = 0;
  amrex::Real eos_diag = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real diag_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real diag_max = -std::numeric_limits<amrex::Real>::infinity();
  amrex::Real row_sum_linf = amrex::Real(0.0);
  amrex::Real diag_dominance_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real gershgorin_lower_bound = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real uniform_h_reference = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real uniform_h_variation_linf = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real uniform_h_theta_pi_phi_pi_symbol =
      std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real uniform_h_theta_pi_over_2_phi_pi_over_2_symbol =
      std::numeric_limits<amrex::Real>::quiet_NaN();
  std::string pressure_topology = "bdltv20_pressure_row";
  std::string pressure_matrix_construction = "not_assigned";
  std::string collocated_source_map_role = "not_assigned";
  std::string dual_face_density_status = "not_assigned";
  std::string pressure_matrix_stencil = "not_assigned";
  std::string nearest_neighbor_control = "not_evaluated";
  std::string uniform_h_checkerboard_status = "not_evaluated";
};

int pressure_periodic_index(int index, int n)
{
  return (index % n + n) % n;
}

int pressure_cell_row(int i, int j, int nx)
{
  return i + nx * j;
}

int pressure_x_face_count(int nx, int ny, bool periodic_x)
{
  return (periodic_x ? nx : nx + 1) * ny;
}

int pressure_y_face_count(int nx, int ny, bool periodic_y)
{
  return nx * (periodic_y ? ny : ny + 1);
}

int pressure_x_face_index(int face_i, int j, int nx, bool periodic_x)
{
  const int face_count_x = periodic_x ? nx : nx + 1;
  const int mapped_i = periodic_x ? pressure_periodic_index(face_i, nx) : face_i;
  return mapped_i + face_count_x * j;
}

int pressure_y_face_index(int i, int face_j, int nx, int ny, bool periodic_y)
{
  const int mapped_j = periodic_y ? pressure_periodic_index(face_j, ny) : face_j;
  return i + nx * mapped_j;
}

int pressure_x_face_east(int i, int j, int nx, bool periodic_x)
{
  return pressure_x_face_index(periodic_x ? i : i + 1, j, nx, periodic_x);
}

int pressure_x_face_west(int i, int j, int nx, bool periodic_x)
{
  return pressure_x_face_index(periodic_x ? i - 1 : i, j, nx, periodic_x);
}

int pressure_y_face_north(int i, int j, int nx, int ny, bool periodic_y)
{
  return pressure_y_face_index(i, periodic_y ? j : j + 1, nx, ny, periodic_y);
}

int pressure_y_face_south(int i, int j, int nx, int ny, bool periodic_y)
{
  return pressure_y_face_index(i, periodic_y ? j - 1 : j, nx, ny, periodic_y);
}

struct HostGmresResult {
  Eigen::VectorXd solution;
  int iterations = 0;
  bool converged = false;
  bool breakdown = false;
  bool nonfinite = false;
  amrex::Real residual_linf = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real relative_residual_linf = std::numeric_limits<amrex::Real>::infinity();
};

HostGmresResult solve_host_sparse_gmres(const HostSparseMatrix& matrix,
                                        const Eigen::VectorXd& rhs,
                                        const Eigen::VectorXd& initial_guess,
                                        int max_iterations,
                                        amrex::Real tolerance)
{
  HostGmresResult result;
  result.solution = initial_guess;
  if (!rhs.allFinite() || !initial_guess.allFinite()) {
    result.nonfinite = true;
    return result;
  }

  const amrex::Real rhs_norm_inf =
      std::max(rhs.lpNorm<Eigen::Infinity>(), amrex::Real(1.0));
  auto update_best = [&](const Eigen::VectorXd& trial, int iterations) {
    if (!trial.allFinite()) {
      result.nonfinite = true;
      return;
    }
    const Eigen::VectorXd residual = matrix * trial - rhs;
    if (!residual.allFinite()) {
      result.nonfinite = true;
      return;
    }
    const amrex::Real residual_linf = residual.lpNorm<Eigen::Infinity>();
    const amrex::Real relative_linf = residual_linf / rhs_norm_inf;
    if (residual_linf < result.residual_linf) {
      result.solution = trial;
      result.iterations = iterations;
      result.residual_linf = residual_linf;
      result.relative_residual_linf = relative_linf;
    }
    if (residual_linf <= tolerance || relative_linf <= tolerance) {
      result.solution = trial;
      result.iterations = iterations;
      result.residual_linf = residual_linf;
      result.relative_residual_linf = relative_linf;
      result.converged = true;
    }
  };

  update_best(initial_guess, 0);
  if (result.converged || result.nonfinite) {
    return result;
  }

  const int n = static_cast<int>(matrix.rows());
  const int iterations_limit = std::max(0, std::min(max_iterations, n));
  if (iterations_limit == 0) {
    return result;
  }

  const Eigen::VectorXd initial_residual = rhs - matrix * initial_guess;
  if (!initial_residual.allFinite()) {
    result.nonfinite = true;
    return result;
  }
  const amrex::Real beta = initial_residual.norm();
  if (!std::isfinite(beta)) {
    result.nonfinite = true;
    return result;
  }
  if (beta == amrex::Real(0.0)) {
    result.converged = true;
    result.residual_linf = amrex::Real(0.0);
    result.relative_residual_linf = amrex::Real(0.0);
    return result;
  }

  std::vector<Eigen::VectorXd> basis;
  basis.reserve(static_cast<std::size_t>(iterations_limit + 1));
  basis.push_back(initial_residual / beta);
  Eigen::MatrixXd hessenberg =
      Eigen::MatrixXd::Zero(iterations_limit + 1, iterations_limit);
  const amrex::Real breakdown_tol =
      std::numeric_limits<amrex::Real>::epsilon() * amrex::Real(100.0) *
      std::max(beta, amrex::Real(1.0));

  for (int j = 0; j < iterations_limit; ++j) {
    Eigen::VectorXd w = matrix * basis[static_cast<std::size_t>(j)];
    if (!w.allFinite()) {
      result.nonfinite = true;
      return result;
    }
    for (int i = 0; i <= j; ++i) {
      hessenberg(i, j) = basis[static_cast<std::size_t>(i)].dot(w);
      w -= hessenberg(i, j) * basis[static_cast<std::size_t>(i)];
    }
    const amrex::Real next_norm = w.norm();
    if (!std::isfinite(next_norm)) {
      result.nonfinite = true;
      return result;
    }
    hessenberg(j + 1, j) = next_norm;
    if (next_norm > breakdown_tol && j + 1 < iterations_limit + 1) {
      basis.push_back(w / next_norm);
    }

    const int rows = j + 2;
    const int cols = j + 1;
    Eigen::VectorXd target = Eigen::VectorXd::Zero(rows);
    target(0) = beta;
    const Eigen::VectorXd y =
        hessenberg.topLeftCorner(rows, cols).colPivHouseholderQr().solve(target);
    if (!y.allFinite()) {
      result.nonfinite = true;
      return result;
    }

    Eigen::VectorXd trial = initial_guess;
    for (int i = 0; i < cols; ++i) {
      trial += y(i) * basis[static_cast<std::size_t>(i)];
    }
    update_best(trial, j + 1);
    if (result.converged || result.nonfinite) {
      return result;
    }
    if (next_norm <= breakdown_tol) {
      result.breakdown = true;
      return result;
    }
  }

  return result;
}
