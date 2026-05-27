// Spatial reconstruction and explicit numerical flux functions.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState primitive_slope(
    const amrex::Array4<const amrex::Real>& arr, int i, int j, int k, int axis, amrex::Real gamma,
    SlopeLimiterKind limiter) noexcept
{
  const PrimitiveState q_m =
      axis == 0 ? to_primitive(load_cons(arr, i - 1, j, k), gamma) : to_primitive(load_cons(arr, i, j - 1, k), gamma);
  const PrimitiveState q_c = to_primitive(load_cons(arr, i, j, k), gamma);
  const PrimitiveState q_p =
      axis == 0 ? to_primitive(load_cons(arr, i + 1, j, k), gamma) : to_primitive(load_cons(arr, i, j + 1, k), gamma);
  return limited_slope(sub_primitive(q_c, q_m), sub_primitive(q_p, q_c), limiter);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState reconstruct_state(
    const amrex::Array4<const amrex::Real>& arr, int i, int j, int k, int axis, int side, int spatial_order,
    amrex::Real gamma, SlopeLimiterKind limiter) noexcept
{
  const PrimitiveState q = to_primitive(load_cons(arr, i, j, k), gamma);
  if (spatial_order <= 1) {
    return q;
  }
  const PrimitiveState slope = primitive_slope(arr, i, j, k, axis, gamma, limiter);
  return add_primitive(q, scale_primitive(amrex::Real(0.5 * side), slope));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState conserved_minmod_slope(
    const amrex::Array4<const amrex::Real>& arr, int i, int j, int k, int axis) noexcept
{
  const ConservedState u_m =
      axis == 0 ? load_cons(arr, i - 1, j, k) : load_cons(arr, i, j - 1, k);
  const ConservedState u_c = load_cons(arr, i, j, k);
  const ConservedState u_p =
      axis == 0 ? load_cons(arr, i + 1, j, k) : load_cons(arr, i, j + 1, k);
  const ConservedState backward = sub_cons(u_c, u_m);
  const ConservedState forward = sub_cons(u_p, u_c);
  return ConservedState{minmod(backward.rho, forward.rho),
                        minmod(backward.mx, forward.mx),
                        minmod(backward.my, forward.my),
                        minmod(backward.E, forward.E)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState reconstruct_bdltv20_eq87_89_conserved_minmod_state(
    const amrex::Array4<const amrex::Real>& arr, int i, int j, int k, int axis, int side,
    int spatial_order, amrex::Real gamma) noexcept
{
  const ConservedState u = load_cons(arr, i, j, k);
  if (spatial_order <= 1) {
    return to_primitive(u, gamma);
  }
  const ConservedState slope = conserved_minmod_slope(arr, i, j, k, axis);
  return to_primitive(add_cons(u, scale_cons(amrex::Real(0.5 * side), slope)), gamma);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState reconstruct_bdltv20_star_state(
    const amrex::Array4<const amrex::Real>& arr, int i, int j, int k, int axis, int side,
    bool use_eq87_89_conserved_minmod, int spatial_order, amrex::Real gamma,
    SlopeLimiterKind limiter) noexcept
{
  if (use_eq87_89_conserved_minmod) {
    return reconstruct_bdltv20_eq87_89_conserved_minmod_state(
        arr, i, j, k, axis, side, spatial_order, gamma);
  }
  return reconstruct_state(arr, i, j, k, axis, side, spatial_order, gamma, limiter);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState physical_flux(
    const PrimitiveState& q, int axis, amrex::Real gamma) noexcept
{
  const ConservedState u = to_conserved(q, gamma);
  if (axis == 0) {
    return ConservedState{u.mx, u.mx * q.u + q.p, u.my * q.u, (u.E + q.p) * q.u};
  }
  return ConservedState{u.my, u.mx * q.v, u.my * q.v + q.p, (u.E + q.p) * q.v};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real max_characteristic_speed(
    const PrimitiveState& q, int axis, amrex::Real gamma) noexcept
{
  return std::abs(normal_velocity(q, axis)) + sound_speed(q, gamma);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState rusanov_flux(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma) noexcept
{
  const ConservedState ul = to_conserved(left, gamma);
  const ConservedState ur = to_conserved(right, gamma);
  const ConservedState fl = physical_flux(left, axis, gamma);
  const ConservedState fr = physical_flux(right, axis, gamma);
  const amrex::Real smax =
      amrex::max(max_characteristic_speed(left, axis, gamma), max_characteristic_speed(right, axis, gamma));
  return sub_cons(scale_cons(amrex::Real(0.5), add_cons(fl, fr)), scale_cons(amrex::Real(0.5) * smax, sub_cons(ur, ul)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState rusanov_flux_with_alpha(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma,
    amrex::Real alpha) noexcept
{
  const ConservedState ul = to_conserved(left, gamma);
  const ConservedState ur = to_conserved(right, gamma);
  const ConservedState fl = physical_flux(left, axis, gamma);
  const ConservedState fr = physical_flux(right, axis, gamma);
  return sub_cons(scale_cons(amrex::Real(0.5), add_cons(fl, fr)),
                  scale_cons(amrex::Real(0.5) * alpha, sub_cons(ur, ul)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState build_hllc_star_state(
    const PrimitiveState& q, int axis, amrex::Real contact_speed, amrex::Real star_density,
    amrex::Real star_energy) noexcept
{
  const amrex::Real tangential_velocity = axis == 0 ? q.v : q.u;
  if (axis == 0) {
    return ConservedState{star_density, star_density * contact_speed, star_density * tangential_velocity, star_energy};
  }
  return ConservedState{star_density, star_density * tangential_velocity, star_density * contact_speed, star_energy};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool hllc_star_state_is_degenerate(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma) noexcept
{
  const amrex::Real rho_l = left.rho;
  const amrex::Real rho_r = right.rho;
  const amrex::Real un_l = normal_velocity(left, axis);
  const amrex::Real un_r = normal_velocity(right, axis);
  const amrex::Real p_l = left.p;
  const amrex::Real p_r = right.p;
  const amrex::Real c_l = sound_speed(left, gamma);
  const amrex::Real c_r = sound_speed(right, gamma);

  const amrex::Real s_l = amrex::min(un_l - c_l, un_r - c_r);
  const amrex::Real s_r = amrex::max(un_l + c_l, un_r + c_r);
  if (!(s_l < s_r)) {
    return true;
  }
  if (s_l >= amrex::Real(0.0) || s_r <= amrex::Real(0.0)) {
    return false;
  }

  const amrex::Real denominator = rho_l * (s_l - un_l) - rho_r * (s_r - un_r);
  if (std::abs(denominator) <= amrex::Real(1.0e-12)) {
    return true;
  }

  const amrex::Real s_m =
      (p_r - p_l + rho_l * un_l * (s_l - un_l) - rho_r * un_r * (s_r - un_r)) / denominator;
  const amrex::Real left_star_denominator = s_l - s_m;
  const amrex::Real right_star_denominator = s_r - s_m;
  return std::abs(left_star_denominator) <= amrex::Real(1.0e-12) ||
         std::abs(right_star_denominator) <= amrex::Real(1.0e-12);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState hllc_flux(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma) noexcept
{
  const ConservedState ul = to_conserved(left, gamma);
  const ConservedState ur = to_conserved(right, gamma);
  const ConservedState fl = physical_flux(left, axis, gamma);
  const ConservedState fr = physical_flux(right, axis, gamma);

  const amrex::Real rho_l = left.rho;
  const amrex::Real rho_r = right.rho;
  const amrex::Real un_l = normal_velocity(left, axis);
  const amrex::Real un_r = normal_velocity(right, axis);
  const amrex::Real p_l = left.p;
  const amrex::Real p_r = right.p;
  const amrex::Real c_l = sound_speed(left, gamma);
  const amrex::Real c_r = sound_speed(right, gamma);

  const amrex::Real s_l = amrex::min(un_l - c_l, un_r - c_r);
  const amrex::Real s_r = amrex::max(un_l + c_l, un_r + c_r);
  if (!(s_l < s_r)) {
    return rusanov_flux(left, right, axis, gamma);
  }
  if (s_l >= amrex::Real(0.0)) {
    return fl;
  }
  if (s_r <= amrex::Real(0.0)) {
    return fr;
  }

  const amrex::Real denominator = rho_l * (s_l - un_l) - rho_r * (s_r - un_r);
  if (std::abs(denominator) <= amrex::Real(1.0e-12)) {
    return rusanov_flux(left, right, axis, gamma);
  }

  const amrex::Real s_m =
      (p_r - p_l + rho_l * un_l * (s_l - un_l) - rho_r * un_r * (s_r - un_r)) / denominator;
  const amrex::Real left_star_denominator = s_l - s_m;
  const amrex::Real right_star_denominator = s_r - s_m;
  if (std::abs(left_star_denominator) <= amrex::Real(1.0e-12) ||
      std::abs(right_star_denominator) <= amrex::Real(1.0e-12)) {
    return rusanov_flux(left, right, axis, gamma);
  }

  const amrex::Real p_star_l = p_l + rho_l * (s_l - un_l) * (s_m - un_l);
  const amrex::Real p_star_r = p_r + rho_r * (s_r - un_r) * (s_m - un_r);
  const amrex::Real p_star = amrex::Real(0.5) * (p_star_l + p_star_r);

  const amrex::Real rho_star_l = rho_l * (s_l - un_l) / left_star_denominator;
  const amrex::Real rho_star_r = rho_r * (s_r - un_r) / right_star_denominator;
  const amrex::Real e_star_l = ((s_l - un_l) * ul.E - p_l * un_l + p_star * s_m) / left_star_denominator;
  const amrex::Real e_star_r = ((s_r - un_r) * ur.E - p_r * un_r + p_star * s_m) / right_star_denominator;

  const ConservedState u_star_l = build_hllc_star_state(left, axis, s_m, rho_star_l, e_star_l);
  const ConservedState u_star_r = build_hllc_star_state(right, axis, s_m, rho_star_r, e_star_r);

  if (s_m >= amrex::Real(0.0)) {
    return add_cons(fl, scale_cons(s_l, sub_cons(u_star_l, ul)));
  }
  return add_cons(fr, scale_cons(s_r, sub_cons(u_star_r, ur)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void low_mach_correct_velocity_states(
    const PrimitiveState& left, const PrimitiveState& right, amrex::Real gamma, PrimitiveState& corrected_left,
    PrimitiveState& corrected_right) noexcept
{
  const amrex::Real c_l = sound_speed(left, gamma);
  const amrex::Real c_r = sound_speed(right, gamma);
  const amrex::Real speed_l = std::sqrt(left.u * left.u + left.v * left.v);
  const amrex::Real speed_r = std::sqrt(right.u * right.u + right.v * right.v);
  const amrex::Real m_l = speed_l / amrex::max(c_l, amrex::Real(1.0e-12));
  const amrex::Real m_r = speed_r / amrex::max(c_r, amrex::Real(1.0e-12));
  const amrex::Real z = amrex::min(amrex::max(amrex::max(m_l, m_r), amrex::Real(0.0)), amrex::Real(1.0));

  corrected_left = left;
  corrected_right = right;
  const amrex::Real u_avg = amrex::Real(0.5) * (left.u + right.u);
  const amrex::Real v_avg = amrex::Real(0.5) * (left.v + right.v);
  corrected_left.u = u_avg + amrex::Real(0.5) * z * (left.u - right.u);
  corrected_right.u = u_avg + amrex::Real(0.5) * z * (right.u - left.u);
  corrected_left.v = v_avg + amrex::Real(0.5) * z * (left.v - right.v);
  corrected_right.v = v_avg + amrex::Real(0.5) * z * (right.v - left.v);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState low_mach_hllc_flux(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma) noexcept
{
  PrimitiveState corrected_left;
  PrimitiveState corrected_right;
  low_mach_correct_velocity_states(left, right, gamma, corrected_left, corrected_right);
  return hllc_flux(corrected_left, corrected_right, axis, gamma);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real clamp01(amrex::Real value) noexcept
{
  return amrex::min(amrex::max(value, amrex::Real(0.0)), amrex::Real(1.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real xie_cubed_sensor(amrex::Real value) noexcept
{
  const amrex::Real f = clamp01(value);
  return f * f * f;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState xie_am_hllc_p_flux(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma,
    amrex::Real sensor_f) noexcept
{
  const ConservedState ul = to_conserved(left, gamma);
  const ConservedState ur = to_conserved(right, gamma);
  const ConservedState fl = physical_flux(left, axis, gamma);
  const ConservedState fr = physical_flux(right, axis, gamma);

  const amrex::Real rho_l = left.rho;
  const amrex::Real rho_r = right.rho;
  const amrex::Real q_l = normal_velocity(left, axis);
  const amrex::Real q_r = normal_velocity(right, axis);
  const amrex::Real p_l = left.p;
  const amrex::Real p_r = right.p;
  const amrex::Real c_l = sound_speed(left, gamma);
  const amrex::Real c_r = sound_speed(right, gamma);
  const amrex::Real s_l = amrex::min(q_l - c_l, q_r - c_r);
  const amrex::Real s_r = amrex::max(q_l + c_l, q_r + c_r);

  if (!(s_l < s_r)) {
    return hllc_flux(left, right, axis, gamma);
  }
  if (s_l >= amrex::Real(0.0)) {
    return fl;
  }
  if (s_r <= amrex::Real(0.0)) {
    return fr;
  }

  const amrex::Real alpha_l = rho_l * (s_l - q_l);
  const amrex::Real alpha_r = rho_r * (s_r - q_r);
  const amrex::Real alpha_denom = alpha_r - alpha_l;
  if (std::abs(alpha_denom) <= amrex::Real(1.0e-12)) {
    return hllc_flux(left, right, axis, gamma);
  }

  const amrex::Real s_star =
      (alpha_r * q_r - alpha_l * q_l + p_l - p_r) / alpha_denom;
  const amrex::Real rho_star_l_denom = s_l - s_star;
  const amrex::Real rho_star_r_denom = s_r - s_star;
  if (std::abs(rho_star_l_denom) <= amrex::Real(1.0e-12) ||
      std::abs(rho_star_r_denom) <= amrex::Real(1.0e-12) ||
      std::abs(alpha_l) <= amrex::Real(1.0e-12) ||
      std::abs(alpha_r) <= amrex::Real(1.0e-12)) {
    return hllc_flux(left, right, axis, gamma);
  }

  const amrex::Real rho_star_l = alpha_l / rho_star_l_denom;
  const amrex::Real rho_star_r = alpha_r / rho_star_r_denom;
  const amrex::Real e_l = ul.E / rho_l;
  const amrex::Real e_r = ur.E / rho_r;
  const amrex::Real e_star_l = e_l + (s_star - q_l) * (s_star + p_l / alpha_l);
  const amrex::Real e_star_r = e_r + (s_star - q_r) * (s_star + p_r / alpha_r);
  const amrex::Real p_star =
      (alpha_r * p_l - alpha_l * p_r - alpha_l * alpha_r * (q_l - q_r)) / alpha_denom;

  const amrex::Real sqrt_rho_l = std::sqrt(rho_l);
  const amrex::Real sqrt_rho_r = std::sqrt(rho_r);
  const amrex::Real roe_denom = sqrt_rho_l + sqrt_rho_r;
  if (!(roe_denom > amrex::Real(0.0))) {
    return hllc_flux(left, right, axis, gamma);
  }
  const amrex::Real u_hat = (sqrt_rho_l * left.u + sqrt_rho_r * right.u) / roe_denom;
  const amrex::Real v_hat = (sqrt_rho_l * left.v + sqrt_rho_r * right.v) / roe_denom;
  const amrex::Real q_hat = axis == 0 ? u_hat : v_hat;
  const amrex::Real h_l = (ul.E + p_l) / rho_l;
  const amrex::Real h_r = (ur.E + p_r) / rho_r;
  const amrex::Real h_hat = (sqrt_rho_l * h_l + sqrt_rho_r * h_r) / roe_denom;
  const amrex::Real c_hat_sq =
      (gamma - amrex::Real(1.0)) * (h_hat - amrex::Real(0.5) * (u_hat * u_hat + v_hat * v_hat));
  if (!(c_hat_sq > amrex::Real(0.0)) || !std::isfinite(c_hat_sq)) {
    return hllc_flux(left, right, axis, gamma);
  }

  const amrex::Real c_hat = std::sqrt(c_hat_sq);
  const amrex::Real m_hat = q_hat / amrex::max(c_hat, amrex::Real(1.0e-12));
  const amrex::Real speed_l = std::sqrt(left.u * left.u + left.v * left.v);
  const amrex::Real speed_r = std::sqrt(right.u * right.u + right.v * right.v);
  const amrex::Real theta =
      clamp01(amrex::max(speed_l / amrex::max(c_l, amrex::Real(1.0e-12)),
                         speed_r / amrex::max(c_r, amrex::Real(1.0e-12))));
  const amrex::Real f = clamp01(sensor_f);
  const amrex::Real p_star_star =
      theta * p_star + (amrex::Real(1.0) - theta) * amrex::Real(0.5) * (p_l + p_r);
  const amrex::Real p_star_star_star = f * p_star_star + (amrex::Real(1.0) - f) * p_star;
  const amrex::Real beta_p =
      (f - amrex::Real(1.0)) * s_l * s_r / (s_r - s_l) * (p_r - p_l) /
      (c_hat_sq * (amrex::Real(1.0) + std::abs(m_hat)));
  const ConservedState phi_p{
      beta_p,
      beta_p * u_hat,
      beta_p * v_hat,
      amrex::Real(0.5) * beta_p * (u_hat * u_hat + v_hat * v_hat)};

  const amrex::Real nx = axis == 0 ? amrex::Real(1.0) : amrex::Real(0.0);
  const amrex::Real ny = axis == 1 ? amrex::Real(1.0) : amrex::Real(0.0);
  const auto star_flux = [=] AMREX_GPU_HOST_DEVICE(const PrimitiveState& q,
                                                   amrex::Real q_normal,
                                                   amrex::Real rho_star,
                                                   amrex::Real e_star) noexcept {
    return ConservedState{
        rho_star * s_star,
        rho_star * s_star * (q.u + nx * (s_star - q_normal)) + nx * p_star_star_star,
        rho_star * s_star * (q.v + ny * (s_star - q_normal)) + ny * p_star_star_star,
        s_star * (rho_star * e_star + p_star)};
  };

  const ConservedState flux =
      s_star >= amrex::Real(0.0)
          ? add_cons(star_flux(left, q_l, rho_star_l, e_star_l), phi_p)
          : add_cons(star_flux(right, q_r, rho_star_r, e_star_r), phi_p);
  if (!std::isfinite(flux.rho) || !std::isfinite(flux.mx) || !std::isfinite(flux.my) ||
      !std::isfinite(flux.E)) {
    return hllc_flux(left, right, axis, gamma);
  }
  return flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState interface_flux(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma,
    RiemannKind riemann, amrex::Real xie_sensor_f = amrex::Real(1.0)) noexcept
{
  switch (riemann) {
    case RiemannKind::Rusanov:
      return rusanov_flux(left, right, axis, gamma);
    case RiemannKind::Hllc:
      return hllc_flux(left, right, axis, gamma);
    case RiemannKind::LowMachHllc:
      return low_mach_hllc_flux(left, right, axis, gamma);
    case RiemannKind::XieAmHllcP:
      return xie_am_hllc_p_flux(left, right, axis, gamma, xie_sensor_f);
  }
  return hllc_flux(left, right, axis, gamma);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real xie_pressure_ratio_from_states(
    const PrimitiveState& left, const PrimitiveState& right) noexcept
{
  const amrex::Real floor = amrex::Real(1.0e-14);
  const amrex::Real p_l = std::isfinite(left.p) ? amrex::max(left.p, floor) : floor;
  const amrex::Real p_r = std::isfinite(right.p) ? amrex::max(right.p, floor) : floor;
  return clamp01(amrex::min(p_l / p_r, p_r / p_l));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real xie_x_face_pressure_ratio(
    const amrex::Array4<const amrex::Real>& arr, int face_i, int j, int k, int spatial_order,
    amrex::Real gamma, SlopeLimiterKind limiter) noexcept
{
  const PrimitiveState left =
      reconstruct_state(arr, face_i - 1, j, k, 0, +1, spatial_order, gamma, limiter);
  const PrimitiveState right =
      reconstruct_state(arr, face_i, j, k, 0, -1, spatial_order, gamma, limiter);
  return xie_pressure_ratio_from_states(left, right);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real xie_y_face_pressure_ratio(
    const amrex::Array4<const amrex::Real>& arr, int i, int face_j, int k, int spatial_order,
    amrex::Real gamma, SlopeLimiterKind limiter) noexcept
{
  const PrimitiveState left =
      reconstruct_state(arr, i, face_j - 1, k, 1, +1, spatial_order, gamma, limiter);
  const PrimitiveState right =
      reconstruct_state(arr, i, face_j, k, 1, -1, spatial_order, gamma, limiter);
  return xie_pressure_ratio_from_states(left, right);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real xie_multiface_pressure_sensor(
    const amrex::Array4<const amrex::Real>& arr, int face_i, int face_j, int k, int axis,
    int spatial_order, amrex::Real gamma, SlopeLimiterKind limiter) noexcept
{
  amrex::Real f = amrex::Real(1.0);
  if (axis == 0) {
    f = xie_x_face_pressure_ratio(arr, face_i, face_j, k, spatial_order, gamma, limiter);
    f = amrex::min(f, xie_y_face_pressure_ratio(arr, face_i - 1, face_j, k, spatial_order, gamma, limiter));
    f = amrex::min(f, xie_y_face_pressure_ratio(arr, face_i - 1, face_j + 1, k, spatial_order, gamma, limiter));
    f = amrex::min(f, xie_y_face_pressure_ratio(arr, face_i, face_j, k, spatial_order, gamma, limiter));
    f = amrex::min(f, xie_y_face_pressure_ratio(arr, face_i, face_j + 1, k, spatial_order, gamma, limiter));
  } else {
    f = xie_y_face_pressure_ratio(arr, face_i, face_j, k, spatial_order, gamma, limiter);
    f = amrex::min(f, xie_x_face_pressure_ratio(arr, face_i, face_j - 1, k, spatial_order, gamma, limiter));
    f = amrex::min(f, xie_x_face_pressure_ratio(arr, face_i + 1, face_j - 1, k, spatial_order, gamma, limiter));
    f = amrex::min(f, xie_x_face_pressure_ratio(arr, face_i, face_j, k, spatial_order, gamma, limiter));
    f = amrex::min(f, xie_x_face_pressure_ratio(arr, face_i + 1, face_j, k, spatial_order, gamma, limiter));
  }
  return xie_cubed_sensor(f);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState reconstructed_interface_flux(
    const amrex::Array4<const amrex::Real>& arr, int face_i, int face_j, int k, int axis,
    int spatial_order, amrex::Real gamma, SlopeLimiterKind limiter, RiemannKind riemann) noexcept
{
  PrimitiveState left;
  PrimitiveState right;
  if (axis == 0) {
    left = reconstruct_state(arr, face_i - 1, face_j, k, 0, +1, spatial_order, gamma, limiter);
    right = reconstruct_state(arr, face_i, face_j, k, 0, -1, spatial_order, gamma, limiter);
  } else {
    left = reconstruct_state(arr, face_i, face_j - 1, k, 1, +1, spatial_order, gamma, limiter);
    right = reconstruct_state(arr, face_i, face_j, k, 1, -1, spatial_order, gamma, limiter);
  }
  const amrex::Real xie_sensor_f =
      riemann == RiemannKind::XieAmHllcP
          ? xie_multiface_pressure_sensor(arr, face_i, face_j, k, axis, spatial_order, gamma, limiter)
          : amrex::Real(1.0);
  return interface_flux(left, right, axis, gamma, riemann, xie_sensor_f);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState pressure_split_explicit_flux(
    const PrimitiveState& q, int axis) noexcept
{
  const amrex::Real mx = q.rho * q.u;
  const amrex::Real my = q.rho * q.v;
  const amrex::Real kinetic = amrex::Real(0.5) * q.rho * (q.u * q.u + q.v * q.v);
  if (axis == 0) {
    return ConservedState{mx, mx * q.u, my * q.u, kinetic * q.u};
  }
  return ConservedState{my, mx * q.v, my * q.v, kinetic * q.v};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real pressure_split_lf_speed(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma,
    ImexPredictorDissipationKind dissipation) noexcept
{
  const amrex::Real material_speed =
      amrex::max(std::abs(normal_velocity(left, axis)), std::abs(normal_velocity(right, axis)));
  if (dissipation == ImexPredictorDissipationKind::Acoustic) {
    return amrex::max(max_characteristic_speed(left, axis, gamma), max_characteristic_speed(right, axis, gamma));
  }
  return material_speed;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState imex_pressure_split_lf_flux(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma,
    ImexPredictorDissipationKind dissipation) noexcept
{
  const ConservedState ul = to_conserved(left, gamma);
  const ConservedState ur = to_conserved(right, gamma);
  const ConservedState fl = pressure_split_explicit_flux(left, axis);
  const ConservedState fr = pressure_split_explicit_flux(right, axis);
  const amrex::Real smax = pressure_split_lf_speed(left, right, axis, gamma, dissipation);
  ConservedState flux =
      sub_cons(scale_cons(amrex::Real(0.5), add_cons(fl, fr)), scale_cons(amrex::Real(0.5) * smax, sub_cons(ur, ul)));
  const amrex::Real momentum_coeff = amrex::Real(0.5) * (amrex::Real(3.0) - gamma);
  flux.mx *= momentum_coeff;
  flux.my *= momentum_coeff;
  flux.E = amrex::Real(0.0);
  return flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState bdltv20_pressure_split_lf_flux(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma,
    ImexPredictorDissipationKind dissipation) noexcept
{
  const ConservedState ul = to_conserved(left, gamma);
  const ConservedState ur = to_conserved(right, gamma);
  const ConservedState fl = pressure_split_explicit_flux(left, axis);
  const ConservedState fr = pressure_split_explicit_flux(right, axis);
  const amrex::Real smax = pressure_split_lf_speed(left, right, axis, gamma, dissipation);
  return sub_cons(scale_cons(amrex::Real(0.5), add_cons(fl, fr)),
                  scale_cons(amrex::Real(0.5) * smax, sub_cons(ur, ul)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ConservedState bdltv20_t1_s2_pressure_split_lf_flux(
    const PrimitiveState& left, const PrimitiveState& right, int axis, amrex::Real gamma,
    ImexPredictorDissipationKind dissipation) noexcept
{
  const ConservedState ul = to_conserved(left, gamma);
  const ConservedState ur = to_conserved(right, gamma);
  const ConservedState fl = pressure_split_explicit_flux(left, axis);
  const ConservedState fr = pressure_split_explicit_flux(right, axis);
  const amrex::Real smax = pressure_split_lf_speed(left, right, axis, gamma, dissipation);
  ConservedState flux =
      sub_cons(scale_cons(amrex::Real(0.5), add_cons(fl, fr)),
               scale_cons(amrex::Real(0.5) * smax, sub_cons(ur, ul)));

  // The clean T1/S2 path transports kinetic energy in the explicit
  // Toro-Vazquez energy flux; avoid diffusing the large pressure/internal-energy
  // part of total energy in this clean-path star update.
  flux.E = amrex::Real(0.5) * (fl.E + fr.E) -
           amrex::Real(0.5) * smax *
               (kinetic_energy_density(right) - kinetic_energy_density(left));
  return flux;
}

