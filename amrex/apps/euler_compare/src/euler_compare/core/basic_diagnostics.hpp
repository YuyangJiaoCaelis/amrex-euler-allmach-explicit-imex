// Basic admissibility records shared by explicit and IMEX rows.
// Diagnostic records written by the solver and plotting workflows.
struct Diagnostics {
  amrex::Real rho_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real pressure_min = std::numeric_limits<amrex::Real>::infinity();
  int nonfinite_count = 0;
};

struct AdmissibilityDiagnostics {
  amrex::Real rho_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real pressure_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real internal_energy_min = std::numeric_limits<amrex::Real>::infinity();
  int nonfinite_count = 0;
  int density_failure_count = 0;
  int pressure_failure_count = 0;
  int internal_energy_failure_count = 0;
};

struct HighTrialRunDiagnostics {
  amrex::Real rho_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real pressure_min = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real internal_energy_min = std::numeric_limits<amrex::Real>::infinity();
  int nonfinite_step_count = 0;
  int density_failure_step_count = 0;
  int pressure_failure_step_count = 0;
  int internal_energy_failure_step_count = 0;
};
