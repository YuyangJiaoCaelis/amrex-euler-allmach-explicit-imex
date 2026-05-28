// String parsers and small configuration predicates.
// ParmParse readers, option parsers, and consistency checks.
std::string lower_string(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

ProblemKind parse_problem(const std::string& value)
{
  const std::string lower = lower_string(value);
  if (lower == "advection_blob" || lower == "blob") {
    return ProblemKind::AdvectionBlob;
  }
  if (lower == "toro1" || lower == "toro_test_1" || lower == "sod" || lower == "sod_x") {
    return ProblemKind::Toro1;
  }
  if (lower == "gresho" || lower == "gresho_vortex") {
    return ProblemKind::GreshoVortex;
  }
  if (lower == "riemann_quadrant" || lower == "quadrant_riemann" ||
      lower == "liska_wendroff_quadrant") {
    return ProblemKind::RiemannQuadrant;
  }
  if (lower == "isentropic_vortex" || lower == "bdltv20_isentropic_vortex") {
    return ProblemKind::IsentropicVortex;
  }
  if (lower == "shock_density_bubble_2d" ||
      lower == "same_gamma_shock_density_bubble" ||
      lower == "shock_density_bubble") {
    return ProblemKind::ShockDensityBubble2D;
  }
  if (lower == "shock_density_bubble_cylindrical_clawpack" ||
      lower == "shock_density_bubble_cylindrical" ||
      lower == "clawpack_shock_density_bubble") {
    return ProblemKind::ShockDensityBubbleCylindrical;
  }
  amrex::Abort("Unknown euler.problem=" + value);
  return ProblemKind::AdvectionBlob;
}

RiemannKind parse_riemann(const std::string& value)
{
  const std::string lower = lower_string(value);
  if (lower == "rusanov" || lower == "local_lax_friedrichs" || lower == "llf") {
    return RiemannKind::Rusanov;
  }
  if (lower == "hllc") {
    return RiemannKind::Hllc;
  }
  if (lower == "lowmach_hllc" || lower == "lm_hllc") {
    return RiemannKind::LowMachHllc;
  }
  if (lower == "xie_am_hllc_p" || lower == "am_hllc_p" || lower == "xie_hllc_p") {
    return RiemannKind::XieAmHllcP;
  }
  amrex::Abort("Unknown euler.riemann=" + value);
  return RiemannKind::Hllc;
}

MethodKind parse_method(const std::string& value)
{
  const std::string lower = lower_string(value);
  if (lower == "explicit") {
    return MethodKind::Explicit;
  }
  if (lower == "imex") {
    return MethodKind::Imex;
  }
  amrex::Abort("Unknown euler.method=" + value);
  return MethodKind::Explicit;
}

ImexFormKind parse_imex_form(const std::string& value)
{
  const std::string lower = lower_string(value);
  if (lower == "bdltv20_t1_s2_source_map_picard" ||
      lower == "bdltv20-t1-s2-source-map-picard") {
    return ImexFormKind::Bdltv20T1S2SourceMapPicard;
  }
  amrex::Abort("Unknown euler.imex_form=" + value +
               "; this code pack supports bdltv20_t1_s2_source_map_picard.");
  return ImexFormKind::Bdltv20T1S2SourceMapPicard;
}

SlopeLimiterKind parse_slope_limiter(const std::string& value)
{
  const std::string lower = lower_string(value);
  if (lower == "minmod") {
    return SlopeLimiterKind::Minmod;
  }
  if (lower == "mc" || lower == "monotonized_central" || lower == "monotonized-central") {
    return SlopeLimiterKind::MonotonizedCentral;
  }
  if (lower == "vanleer" || lower == "van_leer" || lower == "van-leer") {
    return SlopeLimiterKind::VanLeer;
  }
  amrex::Abort("Unknown euler.slope_limiter=" + value);
  return SlopeLimiterKind::Minmod;
}

ImexPredictorDissipationKind parse_imex_predictor_dissipation(const std::string& value)
{
  const std::string lower = lower_string(value);
  if (lower == "material" || lower == "material_speed") {
    return ImexPredictorDissipationKind::Material;
  }
  if (lower == "acoustic" || lower == "acoustic_speed") {
    return ImexPredictorDissipationKind::Acoustic;
  }
  amrex::Abort("Unknown euler.imex_predictor_dissipation=" + value);
  return ImexPredictorDissipationKind::Material;
}

const char* to_string(ProblemKind problem)
{
  switch (problem) {
    case ProblemKind::AdvectionBlob:
      return "advection_blob";
    case ProblemKind::Toro1:
      return "toro1";
    case ProblemKind::GreshoVortex:
      return "gresho_vortex";
    case ProblemKind::RiemannQuadrant:
      return "riemann_quadrant";
    case ProblemKind::IsentropicVortex:
      return "isentropic_vortex";
    case ProblemKind::ShockDensityBubble2D:
      return "shock_density_bubble_2d";
    case ProblemKind::ShockDensityBubbleCylindrical:
      return "shock_density_bubble_cylindrical_clawpack";
  }
  return "advection_blob";
}

const char* to_string(RiemannKind riemann)
{
  switch (riemann) {
    case RiemannKind::Rusanov:
      return "rusanov";
    case RiemannKind::Hllc:
      return "hllc";
    case RiemannKind::LowMachHllc:
      return "lowmach_hllc";
    case RiemannKind::XieAmHllcP:
      return "xie_am_hllc_p";
  }
  return "hllc";
}

const char* to_string(MethodKind method)
{
  switch (method) {
    case MethodKind::Explicit:
      return "explicit";
    case MethodKind::Imex:
      return "imex";
  }
  return "explicit";
}

const char* to_string(ImexFormKind)
{
  return "bdltv20_t1_s2_source_map_picard";
}

const char* to_string(SlopeLimiterKind limiter)
{
  switch (limiter) {
    case SlopeLimiterKind::Minmod:
      return "minmod";
    case SlopeLimiterKind::MonotonizedCentral:
      return "mc";
    case SlopeLimiterKind::VanLeer:
      return "vanleer";
  }
  return "minmod";
}

const char* to_string(ImexPredictorDissipationKind kind)
{
  switch (kind) {
    case ImexPredictorDissipationKind::Material:
      return "material";
    case ImexPredictorDissipationKind::Acoustic:
      return "acoustic";
  }
  return "material";
}

const char* imex_predictor_flux_label(ImexFormKind)
{
  return "bdltv20_pressure_split_lf";
}

const char* imex_route_tag(const RunConfig& cfg)
{
  if (cfg.method != MethodKind::Imex) {
    return "not_applicable";
  }
  if (cfg.bdltv20_paper_t1_s2 != "off") {
    return "bdltv20_direct_paper_driver";
  }
  return "bdltv20_source_map";
}

bool is_bdltv20_t1_s2_source_map_form(ImexFormKind)
{
  return true;
}

bool is_bdltv20_base_source_map_form(ImexFormKind)
{
  return true;
}

bool is_bdltv20_source_map_form(ImexFormKind)
{
  return true;
}

bool is_bdltv20_pressure_unknown_form(ImexFormKind)
{
  return true;
}

bool has_base_bdltv20_paper_controls(const RunConfig& cfg)
{
  return cfg.imex_predictor_dissipation == ImexPredictorDissipationKind::Material &&
         cfg.imex_pressure_stabilization == "off" &&
         cfg.imex_acoustic_startup == 0 &&
         cfg.imex_acoustic_cfl_cap == amrex::Real(0.0);
}

bool has_bdltv20_o1_eq79_selector_controls(const RunConfig&)
{
  return false;
}

std::string pressure_boundary_policy_label(const amrex::Geometry& geom)
{
  const auto dim_label = [&](int dim) {
    return geom.isPeriodic(dim) ? std::string("periodic_pressure_row_wrap")
                                : std::string("outflow_neumann_pressure_ghost_rhs");
  };
  return "x=" + dim_label(0) + ";y=" + dim_label(1);
}
