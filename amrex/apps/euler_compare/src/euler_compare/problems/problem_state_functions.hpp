// Analytic states and problem-specific scalar helper functions.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState gresho_state(
    amrex::Real x, amrex::Real y, const RunConfig& cfg) noexcept
{
  const amrex::Real dx = x - cfg.vortex_cx;
  const amrex::Real dy = y - cfg.vortex_cy;
  const amrex::Real r = std::sqrt(dx * dx + dy * dy);
  const amrex::Real rho = cfg.density_outer;
  const amrex::Real velocity_scale = velocity_scale_from_mach(cfg);

  amrex::Real tangential_factor = amrex::Real(0.0);
  amrex::Real pressure_factor = amrex::Real(-2.0) + amrex::Real(4.0) * std::log(amrex::Real(2.0));

  if (r < amrex::Real(0.2)) {
    tangential_factor = amrex::Real(5.0) * r;
    pressure_factor = amrex::Real(12.5) * r * r;
  } else if (r < amrex::Real(0.4)) {
    tangential_factor = amrex::Real(2.0) - amrex::Real(5.0) * r;
    pressure_factor = amrex::Real(12.5) * r * r +
                      amrex::Real(4.0) *
                          (amrex::Real(1.0) - amrex::Real(5.0) * r - std::log(amrex::Real(0.2)) + std::log(r));
  }

  const amrex::Real vtheta = velocity_scale * tangential_factor;
  const amrex::Real pressure = cfg.pressure + rho * velocity_scale * velocity_scale * pressure_factor;
  if (r < amrex::Real(1.0e-14)) {
    return PrimitiveState{rho, amrex::Real(0.0), amrex::Real(0.0), pressure};
  }
  return PrimitiveState{rho, -vtheta * dy / r, vtheta * dx / r, pressure};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
periodic_nearest_delta(amrex::Real x, amrex::Real center,
                       amrex::Real lo, amrex::Real hi) noexcept
{
  const amrex::Real length = hi - lo;
  amrex::Real delta = x - center;
  if (length > amrex::Real(0.0)) {
    delta -= length * std::floor(delta / length + amrex::Real(0.5));
  }
  return delta;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState isentropic_vortex_state(
    amrex::Real x, amrex::Real y, const RunConfig& cfg,
    amrex::Real prob_lo_x, amrex::Real prob_hi_x,
    amrex::Real prob_lo_y, amrex::Real prob_hi_y,
    amrex::Real time) noexcept
{
  const amrex::Real pi =
      amrex::Real(3.141592653589793238462643383279502884L);
  const amrex::Real gamma = cfg.gamma;
  const amrex::Real strength = cfg.isentropic_vortex_strength;
  const amrex::Real cx = cfg.vortex_cx + cfg.velocity_x * time;
  const amrex::Real cy = cfg.vortex_cy + cfg.velocity_y * time;
  const amrex::Real dx = periodic_nearest_delta(x, cx, prob_lo_x, prob_hi_x);
  const amrex::Real dy = periodic_nearest_delta(y, cy, prob_lo_y, prob_hi_y);
  const amrex::Real r2 = dx * dx + dy * dy;
  const amrex::Real exp_half = std::exp(amrex::Real(0.5) * (amrex::Real(1.0) - r2));
  const amrex::Real exp_full = std::exp(amrex::Real(1.0) - r2);
  const amrex::Real temp =
      amrex::Real(1.0) -
      (gamma - amrex::Real(1.0)) * strength * strength /
          (amrex::Real(8.0) * gamma * pi * pi) * exp_full;
  const amrex::Real rho =
      cfg.density_outer * std::pow(temp, amrex::Real(1.0) / (gamma - amrex::Real(1.0)));
  const amrex::Real pressure =
      cfg.pressure * std::pow(temp, gamma / (gamma - amrex::Real(1.0)));
  const amrex::Real velocity_factor = strength / (amrex::Real(2.0) * pi) * exp_half;
  return PrimitiveState{rho, cfg.velocity_x - velocity_factor * dy,
                        cfg.velocity_y + velocity_factor * dx, pressure};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState
riemann_quadrant_state(amrex::Real x, amrex::Real y, const RunConfig& cfg) noexcept
{
  const bool right = x >= cfg.riemann_interface_x;
  const bool top = y >= cfg.riemann_interface_y;

  // Liska-Wendroff/Clawpack quadrant Riemann state set. The ordering below is
  // upper-right, upper-left, lower-left, lower-right.
  if (right && top) {
    return PrimitiveState{amrex::Real(1.5), amrex::Real(0.0),
                          amrex::Real(0.0), amrex::Real(1.5)};
  }
  if (!right && top) {
    return PrimitiveState{amrex::Real(0.532258064516129),
                          amrex::Real(1.206045378311055),
                          amrex::Real(0.0), amrex::Real(0.3)};
  }
  if (!right && !top) {
    return PrimitiveState{amrex::Real(0.137992831541219),
                          amrex::Real(1.206045378311055),
                          amrex::Real(1.206045378311055),
                          amrex::Real(0.029032258064516)};
  }
  return PrimitiveState{amrex::Real(0.532258064516129), amrex::Real(0.0),
                        amrex::Real(1.206045378311055), amrex::Real(0.3)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real kinetic_energy_density(const PrimitiveState& q) noexcept
{
  return amrex::Real(0.5) * q.rho * (q.u * q.u + q.v * q.v);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bdltv20_specific_enthalpy(
    amrex::Real rho, amrex::Real pressure, amrex::Real gamma) noexcept
{
  return gamma * pressure / ((gamma - amrex::Real(1.0)) * rho);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bdltv20_eq46_face_enthalpy(
    amrex::Real rho_l, amrex::Real p_l, amrex::Real q_l,
    amrex::Real rho_r, amrex::Real p_r, amrex::Real q_r,
    amrex::Real gamma, int& guard_count) noexcept
{
  const amrex::Real h_l = bdltv20_specific_enthalpy(rho_l, p_l, gamma);
  const amrex::Real h_r = bdltv20_specific_enthalpy(rho_r, p_r, gamma);
  const amrex::Real denom = q_l + q_r;
  const amrex::Real scale =
      amrex::max(amrex::Real(1.0), std::abs(q_l) + std::abs(q_r));
  if (std::abs(denom) <= amrex::Real(1.0e-12) * scale) {
    ++guard_count;
    return amrex::Real(0.5) * (h_l + h_r);
  }
  return (h_l * q_l + h_r * q_r) / denom;
}

