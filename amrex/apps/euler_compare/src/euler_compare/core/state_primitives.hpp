// Conserved/primitive state algebra and thermodynamic helpers.
// Initial conditions, exact states, and physical/ghost boundary fills.
amrex::Real cell_volume(const amrex::Geometry& geom)
{
  const auto dx = geom.CellSizeArray();
  return dx[0] * dx[1];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState load_cons(
    const amrex::Array4<const amrex::Real>& arr, int i, int j, int k) noexcept
{
  return ConservedState{arr(i, j, k, Rho), arr(i, j, k, Mx), arr(i, j, k, My), arr(i, j, k, E)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void store_cons(
    const amrex::Array4<amrex::Real>& arr, int i, int j, int k, const ConservedState& u) noexcept
{
  arr(i, j, k, Rho) = u.rho;
  arr(i, j, k, Mx) = u.mx;
  arr(i, j, k, My) = u.my;
  arr(i, j, k, E) = u.E;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState add_cons(const ConservedState& a, const ConservedState& b) noexcept
{
  return ConservedState{a.rho + b.rho, a.mx + b.mx, a.my + b.my, a.E + b.E};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState sub_cons(const ConservedState& a, const ConservedState& b) noexcept
{
  return ConservedState{a.rho - b.rho, a.mx - b.mx, a.my - b.my, a.E - b.E};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState scale_cons(amrex::Real s, const ConservedState& a) noexcept
{
  return ConservedState{s * a.rho, s * a.mx, s * a.my, s * a.E};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState add_primitive(
    const PrimitiveState& a, const PrimitiveState& b) noexcept
{
  return PrimitiveState{a.rho + b.rho, a.u + b.u, a.v + b.v, a.p + b.p};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState sub_primitive(
    const PrimitiveState& a, const PrimitiveState& b) noexcept
{
  return PrimitiveState{a.rho - b.rho, a.u - b.u, a.v - b.v, a.p - b.p};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState scale_primitive(
    amrex::Real s, const PrimitiveState& a) noexcept
{
  return PrimitiveState{s * a.rho, s * a.u, s * a.v, s * a.p};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real minmod(amrex::Real a, amrex::Real b) noexcept
{
  if (a * b <= amrex::Real(0.0)) {
    return amrex::Real(0.0);
  }
  return std::abs(a) < std::abs(b) ? a : b;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState minmod(
    const PrimitiveState& a, const PrimitiveState& b) noexcept
{
  return PrimitiveState{minmod(a.rho, b.rho), minmod(a.u, b.u), minmod(a.v, b.v), minmod(a.p, b.p)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real limited_slope(
    amrex::Real backward, amrex::Real forward, SlopeLimiterKind limiter) noexcept
{
  if (backward * forward <= amrex::Real(0.0)) {
    return amrex::Real(0.0);
  }
  switch (limiter) {
    case SlopeLimiterKind::Minmod:
      return minmod(backward, forward);
    case SlopeLimiterKind::MonotonizedCentral:
      return minmod(minmod(amrex::Real(2.0) * backward, amrex::Real(0.5) * (backward + forward)),
                    amrex::Real(2.0) * forward);
    case SlopeLimiterKind::VanLeer: {
      const amrex::Real denominator = backward + forward;
      if (std::abs(denominator) <= amrex::Real(1.0e-14)) {
        return amrex::Real(0.0);
      }
      return amrex::Real(2.0) * backward * forward / denominator;
    }
  }
  return minmod(backward, forward);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState limited_slope(
    const PrimitiveState& backward, const PrimitiveState& forward, SlopeLimiterKind limiter) noexcept
{
  return PrimitiveState{limited_slope(backward.rho, forward.rho, limiter),
                        limited_slope(backward.u, forward.u, limiter),
                        limited_slope(backward.v, forward.v, limiter),
                        limited_slope(backward.p, forward.p, limiter)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState to_primitive(const ConservedState& u, amrex::Real gamma) noexcept
{
  const amrex::Real inv_rho = amrex::Real(1.0) / u.rho;
  const amrex::Real ux = u.mx * inv_rho;
  const amrex::Real uy = u.my * inv_rho;
  const amrex::Real kinetic = amrex::Real(0.5) * u.rho * (ux * ux + uy * uy);
  const amrex::Real p = (gamma - amrex::Real(1.0)) * (u.E - kinetic);
  return PrimitiveState{u.rho, ux, uy, p};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState to_conserved(const PrimitiveState& q, amrex::Real gamma) noexcept
{
  const amrex::Real kinetic = amrex::Real(0.5) * q.rho * (q.u * q.u + q.v * q.v);
  return ConservedState{q.rho, q.rho * q.u, q.rho * q.v, q.p / (gamma - amrex::Real(1.0)) + kinetic};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState bdltv20_paper_to_primitive(
    const ConservedState& u, amrex::Real gamma, amrex::Real epsilon) noexcept
{
  const amrex::Real inv_rho = amrex::Real(1.0) / u.rho;
  const amrex::Real ux = u.mx * inv_rho;
  const amrex::Real uy = u.my * inv_rho;
  const amrex::Real kinetic = amrex::Real(0.5) * u.rho * (ux * ux + uy * uy);
  const amrex::Real p =
      (gamma - amrex::Real(1.0)) * (u.E - epsilon * kinetic);
  return PrimitiveState{u.rho, ux, uy, p};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState bdltv20_paper_to_conserved(
    const PrimitiveState& q, amrex::Real gamma, amrex::Real epsilon) noexcept
{
  const amrex::Real kinetic =
      amrex::Real(0.5) * q.rho * (q.u * q.u + q.v * q.v);
  return ConservedState{q.rho, q.rho * q.u, q.rho * q.v,
                        q.p / (gamma - amrex::Real(1.0)) +
                            epsilon * kinetic};
}

struct DualCloseoutCheckCell {
  amrex::Real solved_pressure = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real cons_total_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real cons_pressure = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real cons_internal_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real eos_total_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real eos_internal_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real eos_minus_cons_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real linearized_total_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real linearized_pressure = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real linearized_internal_energy = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real full_minus_linearized_kinetic = std::numeric_limits<amrex::Real>::quiet_NaN();
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE DualCloseoutCheckCell compute_dual_closeout_check_cell(
    const ConservedState& old_state, const ConservedState& trial, amrex::Real solved_pressure,
    amrex::Real gamma) noexcept
{
  DualCloseoutCheckCell out;
  out.solved_pressure = solved_pressure;
  out.cons_total_energy = trial.E;
  const bool trial_density_valid = std::isfinite(trial.rho) && trial.rho != amrex::Real(0.0);
  const bool old_density_valid = std::isfinite(old_state.rho) && old_state.rho != amrex::Real(0.0);

  amrex::Real full_kinetic = std::numeric_limits<amrex::Real>::quiet_NaN();
  if (trial_density_valid && std::isfinite(trial.mx) && std::isfinite(trial.my)) {
    full_kinetic =
        amrex::Real(0.5) * (trial.mx * trial.mx + trial.my * trial.my) / trial.rho;
  }
  if (std::isfinite(trial.E) && std::isfinite(full_kinetic)) {
    out.cons_internal_energy = trial.E - full_kinetic;
    out.cons_pressure = (gamma - amrex::Real(1.0)) * out.cons_internal_energy;
  }
  if (std::isfinite(solved_pressure)) {
    out.eos_internal_energy = solved_pressure / (gamma - amrex::Real(1.0));
    if (std::isfinite(full_kinetic)) {
      out.eos_total_energy = out.eos_internal_energy + full_kinetic;
      if (std::isfinite(trial.E)) {
        out.eos_minus_cons_energy = out.eos_total_energy - trial.E;
      }
    }
  }
  if (old_density_valid && std::isfinite(old_state.mx) && std::isfinite(old_state.my) &&
      std::isfinite(trial.mx) && std::isfinite(trial.my) && std::isfinite(solved_pressure) &&
      std::isfinite(full_kinetic)) {
    const amrex::Real linearized_kinetic =
        amrex::Real(0.5) * (old_state.mx * trial.mx + old_state.my * trial.my) / old_state.rho;
    out.linearized_total_energy =
        solved_pressure / (gamma - amrex::Real(1.0)) + linearized_kinetic;
    out.linearized_internal_energy = out.linearized_total_energy - full_kinetic;
    out.linearized_pressure = (gamma - amrex::Real(1.0)) * out.linearized_internal_energy;
    out.full_minus_linearized_kinetic = full_kinetic - linearized_kinetic;
  }
  return out;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real normal_velocity(const PrimitiveState& q, int axis) noexcept
{
  return axis == 0 ? q.u : q.v;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real sound_speed(const PrimitiveState& q, amrex::Real gamma) noexcept
{
  return std::sqrt(gamma * q.p / q.rho);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real reference_sound_speed(const RunConfig& cfg) noexcept
{
  return std::sqrt(cfg.gamma * cfg.pressure / cfg.density_outer);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real velocity_scale_from_mach(const RunConfig& cfg) noexcept
{
  return cfg.mach * reference_sound_speed(cfg);
}

