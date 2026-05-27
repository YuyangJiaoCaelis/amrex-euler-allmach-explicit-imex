#!/usr/bin/env python3

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class RiemannState:
    rho: float
    u: float
    p: float


@dataclass(frozen=True)
class StarRegion:
    p: float
    u: float
    iterations: int
    residual: float


@dataclass(frozen=True)
class WaveInfo:
    left_wave: str
    right_wave: str
    p_star: float
    u_star: float
    left_shock_speed: float | None
    right_shock_speed: float | None
    left_rarefaction_head_speed: float | None
    left_rarefaction_tail_speed: float | None
    right_rarefaction_head_speed: float | None
    right_rarefaction_tail_speed: float | None
    contact_speed: float


def sound_speed(state: RiemannState, gamma: float) -> float:
    if state.rho <= 0.0 or state.p <= 0.0:
        raise ValueError("Riemann states require positive density and pressure.")
    return math.sqrt(gamma * state.p / state.rho)


def pressure_function(p: float, state: RiemannState, gamma: float) -> tuple[float, float]:
    a = sound_speed(state, gamma)
    if p > state.p:
        acoef = 2.0 / ((gamma + 1.0) * state.rho)
        bcoef = (gamma - 1.0) / (gamma + 1.0) * state.p
        root = math.sqrt(acoef / (p + bcoef))
        value = (p - state.p) * root
        derivative = root * (1.0 - 0.5 * (p - state.p) / (p + bcoef))
        return value, derivative

    ratio = p / state.p
    value = 2.0 * a / (gamma - 1.0) * (ratio ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)
    derivative = 1.0 / (state.rho * a) * ratio ** (-(gamma + 1.0) / (2.0 * gamma))
    return value, derivative


def solve_star_region(left: RiemannState, right: RiemannState, gamma: float) -> StarRegion:
    a_l = sound_speed(left, gamma)
    a_r = sound_speed(right, gamma)
    vacuum_limit = 2.0 * (a_l + a_r) / (gamma - 1.0)
    if right.u - left.u >= vacuum_limit:
        raise ValueError("The selected Riemann problem would generate vacuum.")

    p_pvrs = 0.5 * (left.p + right.p) - 0.125 * (right.u - left.u) * (left.rho + right.rho) * (a_l + a_r)
    p = max(1.0e-12, p_pvrs)
    residual = math.inf

    for iteration in range(1, 81):
        f_l, d_l = pressure_function(p, left, gamma)
        f_r, d_r = pressure_function(p, right, gamma)
        residual = f_l + f_r + right.u - left.u
        p_next = p - residual / (d_l + d_r)
        if p_next <= 0.0 or not math.isfinite(p_next):
            p_next = 0.5 * p
        if abs(p_next - p) <= 1.0e-12 * max(1.0, p_next):
            p = p_next
            break
        p = p_next
    else:
        raise RuntimeError("Exact Riemann p_star solve did not converge.")

    f_l, _ = pressure_function(p, left, gamma)
    f_r, _ = pressure_function(p, right, gamma)
    u = 0.5 * (left.u + right.u + f_r - f_l)
    residual = f_l + f_r + right.u - left.u
    return StarRegion(p=p, u=u, iterations=iteration, residual=residual)


def shock_density(state: RiemannState, gamma: float, p_star: float) -> float:
    ratio = p_star / state.p
    return state.rho * ((ratio + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * ratio + 1.0))


def rarefaction_density(state: RiemannState, gamma: float, p_star: float) -> float:
    return state.rho * (p_star / state.p) ** (1.0 / gamma)


def wave_info(left: RiemannState, right: RiemannState, gamma: float, star: StarRegion) -> WaveInfo:
    a_l = sound_speed(left, gamma)
    a_r = sound_speed(right, gamma)
    left_shock_speed = None
    right_shock_speed = None
    left_head = None
    left_tail = None
    right_head = None
    right_tail = None

    if star.p > left.p:
        left_wave = "shock"
        left_shock_speed = left.u - a_l * math.sqrt(
            (gamma + 1.0) / (2.0 * gamma) * star.p / left.p + (gamma - 1.0) / (2.0 * gamma)
        )
    else:
        left_wave = "rarefaction"
        a_star_l = a_l * (star.p / left.p) ** ((gamma - 1.0) / (2.0 * gamma))
        left_head = left.u - a_l
        left_tail = star.u - a_star_l

    if star.p > right.p:
        right_wave = "shock"
        right_shock_speed = right.u + a_r * math.sqrt(
            (gamma + 1.0) / (2.0 * gamma) * star.p / right.p + (gamma - 1.0) / (2.0 * gamma)
        )
    else:
        right_wave = "rarefaction"
        a_star_r = a_r * (star.p / right.p) ** ((gamma - 1.0) / (2.0 * gamma))
        right_head = right.u + a_r
        right_tail = star.u + a_star_r

    return WaveInfo(
        left_wave=left_wave,
        right_wave=right_wave,
        p_star=star.p,
        u_star=star.u,
        left_shock_speed=left_shock_speed,
        right_shock_speed=right_shock_speed,
        left_rarefaction_head_speed=left_head,
        left_rarefaction_tail_speed=left_tail,
        right_rarefaction_head_speed=right_head,
        right_rarefaction_tail_speed=right_tail,
        contact_speed=star.u,
    )


def sample(
    x: float,
    time: float,
    left: RiemannState,
    right: RiemannState,
    gamma: float,
    interface_x: float,
    star: StarRegion,
) -> RiemannState:
    if time <= 0.0:
        return left if x < interface_x else right

    xi = (x - interface_x) / time
    a_l = sound_speed(left, gamma)
    a_r = sound_speed(right, gamma)

    if xi <= star.u:
        if star.p > left.p:
            speed = left.u - a_l * math.sqrt(
                (gamma + 1.0) / (2.0 * gamma) * star.p / left.p + (gamma - 1.0) / (2.0 * gamma)
            )
            if xi <= speed:
                return left
            return RiemannState(shock_density(left, gamma, star.p), star.u, star.p)

        speed_head = left.u - a_l
        a_star = a_l * (star.p / left.p) ** ((gamma - 1.0) / (2.0 * gamma))
        speed_tail = star.u - a_star
        if xi <= speed_head:
            return left
        if xi >= speed_tail:
            return RiemannState(rarefaction_density(left, gamma, star.p), star.u, star.p)
        u = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * left.u + xi)
        a = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * (left.u - xi))
        rho = left.rho * (a / a_l) ** (2.0 / (gamma - 1.0))
        pressure = left.p * (a / a_l) ** (2.0 * gamma / (gamma - 1.0))
        return RiemannState(rho, u, pressure)

    if star.p > right.p:
        speed = right.u + a_r * math.sqrt(
            (gamma + 1.0) / (2.0 * gamma) * star.p / right.p + (gamma - 1.0) / (2.0 * gamma)
        )
        if xi >= speed:
            return right
        return RiemannState(shock_density(right, gamma, star.p), star.u, star.p)

    speed_head = right.u + a_r
    a_star = a_r * (star.p / right.p) ** ((gamma - 1.0) / (2.0 * gamma))
    speed_tail = star.u + a_star
    if xi >= speed_head:
        return right
    if xi <= speed_tail:
        return RiemannState(rarefaction_density(right, gamma, star.p), star.u, star.p)
    u = 2.0 / (gamma + 1.0) * (-a_r + 0.5 * (gamma - 1.0) * right.u + xi)
    a = 2.0 / (gamma + 1.0) * (a_r - 0.5 * (gamma - 1.0) * (right.u - xi))
    rho = right.rho * (a / a_r) ** (2.0 / (gamma - 1.0))
    pressure = right.p * (a / a_r) ** (2.0 * gamma / (gamma - 1.0))
    return RiemannState(rho, u, pressure)


def wave_locations(info: WaveInfo, interface_x: float, time: float) -> dict[str, float | str]:
    def loc(speed: float | None) -> float | str:
        if speed is None:
            return ""
        return interface_x + speed * time

    return {
        "exact_left_wave": info.left_wave,
        "exact_right_wave": info.right_wave,
        "exact_contact_x": interface_x + info.contact_speed * time,
        "exact_left_shock_x": loc(info.left_shock_speed),
        "exact_right_shock_x": loc(info.right_shock_speed),
        "exact_left_rarefaction_head_x": loc(info.left_rarefaction_head_speed),
        "exact_left_rarefaction_tail_x": loc(info.left_rarefaction_tail_speed),
        "exact_right_rarefaction_head_x": loc(info.right_rarefaction_head_speed),
        "exact_right_rarefaction_tail_x": loc(info.right_rarefaction_tail_speed),
    }


def major_wave_positions(info: WaveInfo, interface_x: float, time: float) -> list[float]:
    locations = wave_locations(info, interface_x, time)
    positions: list[float] = []
    for key, value in locations.items():
        if key.startswith("exact_") and key.endswith("_x") and isinstance(value, float):
            positions.append(value)
    return positions

