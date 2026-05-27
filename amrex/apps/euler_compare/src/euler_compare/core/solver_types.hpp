// Core enums and primitive/conserved state records used by the Euler solvers.
enum class ProblemKind {
  AdvectionBlob,
  Toro1,
  GreshoVortex,
  RiemannQuadrant,
  IsentropicVortex,
  ShockDensityBubble2D,
  ShockDensityBubbleCylindrical
};
enum class RiemannKind { Rusanov, Hllc, LowMachHllc, XieAmHllcP };
enum class MethodKind { Explicit, Imex };
enum class ImexFormKind { Bdltv20T1S2SourceMapPicard };
enum class SlopeLimiterKind { Minmod, MonotonizedCentral, VanLeer };
enum class ImexPredictorDissipationKind { Material, Acoustic };

struct PrimitiveState {
  amrex::Real rho = 0.0;
  amrex::Real u = 0.0;
  amrex::Real v = 0.0;
  amrex::Real p = 0.0;
};

struct ConservedState {
  amrex::Real rho = 0.0;
  amrex::Real mx = 0.0;
  amrex::Real my = 0.0;
  amrex::Real E = 0.0;
};

