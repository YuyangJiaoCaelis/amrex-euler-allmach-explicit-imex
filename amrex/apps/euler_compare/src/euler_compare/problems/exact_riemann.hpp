// Toro states and exact one-dimensional Riemann sampling.
ToroState toro_state(int test_id)
{
  switch (test_id) {
    case 1:
      return ToroState{1.0, 0.0, 0.0, 1.0, 0.125, 0.0, 0.0, 0.1};
    case 2:
      return ToroState{1.0, -2.0, 0.0, 0.4, 1.0, 2.0, 0.0, 0.4};
    case 3:
      return ToroState{1.0, 0.0, 0.0, 1000.0, 1.0, 0.0, 0.0, 0.01};
    case 4:
      return ToroState{1.0, 0.0, 0.0, 0.01, 1.0, 0.0, 0.0, 100.0};
    case 5:
      return ToroState{5.99924, 19.5975, 0.0, 460.894, 5.99242, -6.19633, 0.0, 46.0950};
    case 93:
      return ToroState{0.445, 1.698, 0.0, 3.528, 0.5, 0.0, 0.0, 0.571};
    case 94:
      return ToroState{1000.0, 1.0, 0.0, 1.0e5, 0.01, 1.0, 0.0, 1.0e5};
    case 95:
      return ToroState{1.0, 0.0, 0.0, 1.01, 1.0, 0.0, 0.0, 1.0};
    case 96:
      return ToroState{1.0, 0.0, 0.0, 1.0001, 1.0, 0.0, 0.0, 1.0};
    case 97:
      return ToroState{1.0, 2.0, 0.0, 0.1, 1.0, -2.0, 0.0, 0.1};
    case 90:
      return ToroState{1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0};
    case 91:
      return ToroState{2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0};
    case 92:
      return ToroState{2.0, 0.1, 0.0, 1.0, 1.0, 0.1, 0.0, 1.0};
    default:
      amrex::Abort("Unsupported Toro test id.");
  }
  return ToroState{};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
exact_riemann_sound_speed(const PrimitiveState& q, amrex::Real gamma) noexcept
{
  return std::sqrt(gamma * q.p / q.rho);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
exact_riemann_pressure_function(amrex::Real p, const PrimitiveState& q,
                                amrex::Real gamma, amrex::Real& value,
                                amrex::Real& derivative) noexcept
{
  const amrex::Real a = exact_riemann_sound_speed(q, gamma);
  if (p > q.p) {
    const amrex::Real acoef = amrex::Real(2.0) / ((gamma + amrex::Real(1.0)) * q.rho);
    const amrex::Real bcoef =
        (gamma - amrex::Real(1.0)) / (gamma + amrex::Real(1.0)) * q.p;
    const amrex::Real root = std::sqrt(acoef / (p + bcoef));
    value = (p - q.p) * root;
    derivative = root * (amrex::Real(1.0) -
                         amrex::Real(0.5) * (p - q.p) / (p + bcoef));
    return;
  }

  const amrex::Real ratio = p / q.p;
  value = amrex::Real(2.0) * a / (gamma - amrex::Real(1.0)) *
          (std::pow(ratio, (gamma - amrex::Real(1.0)) /
                               (amrex::Real(2.0) * gamma)) -
           amrex::Real(1.0));
  derivative = amrex::Real(1.0) / (q.rho * a) *
               std::pow(ratio, -(gamma + amrex::Real(1.0)) /
                                    (amrex::Real(2.0) * gamma));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
exact_riemann_star(const PrimitiveState& left, const PrimitiveState& right,
                   amrex::Real gamma, amrex::Real& p_star,
                   amrex::Real& u_star) noexcept
{
  const amrex::Real a_l = exact_riemann_sound_speed(left, gamma);
  const amrex::Real a_r = exact_riemann_sound_speed(right, gamma);
  amrex::Real p = amrex::max(
      amrex::Real(1.0e-12),
      amrex::Real(0.5) * (left.p + right.p) -
          amrex::Real(0.125) * (right.u - left.u) * (left.rho + right.rho) *
              (a_l + a_r));

  for (int iter = 0; iter < 80; ++iter) {
    amrex::Real f_l = amrex::Real(0.0);
    amrex::Real d_l = amrex::Real(0.0);
    amrex::Real f_r = amrex::Real(0.0);
    amrex::Real d_r = amrex::Real(0.0);
    exact_riemann_pressure_function(p, left, gamma, f_l, d_l);
    exact_riemann_pressure_function(p, right, gamma, f_r, d_r);
    const amrex::Real residual = f_l + f_r + right.u - left.u;
    amrex::Real p_next = p - residual / (d_l + d_r);
    if (!(p_next > amrex::Real(0.0)) || !std::isfinite(p_next)) {
      p_next = amrex::Real(0.5) * p;
    }
    if (std::abs(p_next - p) <=
        amrex::Real(1.0e-12) * amrex::max(amrex::Real(1.0), p_next)) {
      p = p_next;
      break;
    }
    p = p_next;
  }

  amrex::Real f_l = amrex::Real(0.0);
  amrex::Real d_l = amrex::Real(0.0);
  amrex::Real f_r = amrex::Real(0.0);
  amrex::Real d_r = amrex::Real(0.0);
  exact_riemann_pressure_function(p, left, gamma, f_l, d_l);
  exact_riemann_pressure_function(p, right, gamma, f_r, d_r);
  p_star = p;
  u_star = amrex::Real(0.5) * (left.u + right.u + f_r - f_l);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PrimitiveState
exact_riemann_sample(amrex::Real x, amrex::Real time, const PrimitiveState& left,
                     const PrimitiveState& right, amrex::Real gamma,
                     amrex::Real interface_x) noexcept
{
  if (time <= amrex::Real(0.0)) {
    return x < interface_x ? left : right;
  }

  amrex::Real p_star = amrex::Real(0.0);
  amrex::Real u_star = amrex::Real(0.0);
  exact_riemann_star(left, right, gamma, p_star, u_star);

  const amrex::Real xi = (x - interface_x) / time;
  const amrex::Real a_l = exact_riemann_sound_speed(left, gamma);
  const amrex::Real a_r = exact_riemann_sound_speed(right, gamma);
  const amrex::Real gm1 = gamma - amrex::Real(1.0);
  const amrex::Real gp1 = gamma + amrex::Real(1.0);

  const auto shock_density = [=] AMREX_GPU_HOST_DEVICE(const PrimitiveState& q,
                                                       amrex::Real p) noexcept {
    const amrex::Real ratio = p / q.p;
    return q.rho * ((ratio + gm1 / gp1) / ((gm1 / gp1) * ratio + amrex::Real(1.0)));
  };
  const auto rarefaction_density = [=] AMREX_GPU_HOST_DEVICE(const PrimitiveState& q,
                                                            amrex::Real p) noexcept {
    return q.rho * std::pow(p / q.p, amrex::Real(1.0) / gamma);
  };

  if (xi <= u_star) {
    if (p_star > left.p) {
      const amrex::Real speed =
          left.u - a_l * std::sqrt(gp1 / (amrex::Real(2.0) * gamma) * p_star / left.p +
                                   gm1 / (amrex::Real(2.0) * gamma));
      if (xi <= speed) {
        return left;
      }
      return PrimitiveState{shock_density(left, p_star), u_star, left.v, p_star};
    }

    const amrex::Real speed_head = left.u - a_l;
    const amrex::Real a_star =
        a_l * std::pow(p_star / left.p, gm1 / (amrex::Real(2.0) * gamma));
    const amrex::Real speed_tail = u_star - a_star;
    if (xi <= speed_head) {
      return left;
    }
    if (xi >= speed_tail) {
      return PrimitiveState{rarefaction_density(left, p_star), u_star, left.v, p_star};
    }
    const amrex::Real u =
        amrex::Real(2.0) / gp1 *
        (a_l + amrex::Real(0.5) * gm1 * left.u + xi);
    const amrex::Real a =
        amrex::Real(2.0) / gp1 *
        (a_l + amrex::Real(0.5) * gm1 * (left.u - xi));
    const amrex::Real rho = left.rho * std::pow(a / a_l, amrex::Real(2.0) / gm1);
    const amrex::Real p =
        left.p * std::pow(a / a_l, amrex::Real(2.0) * gamma / gm1);
    return PrimitiveState{rho, u, left.v, p};
  }

  if (p_star > right.p) {
    const amrex::Real speed =
        right.u + a_r * std::sqrt(gp1 / (amrex::Real(2.0) * gamma) * p_star / right.p +
                                  gm1 / (amrex::Real(2.0) * gamma));
    if (xi >= speed) {
      return right;
    }
    return PrimitiveState{shock_density(right, p_star), u_star, right.v, p_star};
  }

  const amrex::Real speed_head = right.u + a_r;
  const amrex::Real a_star =
      a_r * std::pow(p_star / right.p, gm1 / (amrex::Real(2.0) * gamma));
  const amrex::Real speed_tail = u_star + a_star;
  if (xi >= speed_head) {
    return right;
  }
  if (xi <= speed_tail) {
    return PrimitiveState{rarefaction_density(right, p_star), u_star, right.v, p_star};
  }
  const amrex::Real u =
      amrex::Real(2.0) / gp1 *
      (-a_r + amrex::Real(0.5) * gm1 * right.u + xi);
  const amrex::Real a =
      amrex::Real(2.0) / gp1 *
      (a_r - amrex::Real(0.5) * gm1 * (right.u - xi));
  const amrex::Real rho = right.rho * std::pow(a / a_r, amrex::Real(2.0) / gm1);
  const amrex::Real p =
      right.p * std::pow(a / a_r, amrex::Real(2.0) * gamma / gm1);
  return PrimitiveState{rho, u, right.v, p};
}
