// Problem-specific metric records and exact Riemann state parameters.
struct GreshoMetrics {
  amrex::Real density_l1_error = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real velocity_l1_error = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real pressure_l1_error = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real pressure_perturbation_l1_error = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real pressure_perturbation_l1_relative_error = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real kinetic_energy_initial = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real kinetic_energy_final = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real kinetic_energy_ratio = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real reference_sound_speed = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real velocity_scale = std::numeric_limits<amrex::Real>::quiet_NaN();
};

struct ImexTimeStepDiagnostics {
  amrex::Real material_rate_max = amrex::Real(0.0);
  amrex::Real acoustic_rate_max = amrex::Real(0.0);
  int acoustic_startup_used = 0;
  int acoustic_cap_used = 0;
};


struct ToroState {
  amrex::Real rho_l = 1.0;
  amrex::Real u_l = 0.0;
  amrex::Real v_l = 0.0;
  amrex::Real p_l = 1.0;
  amrex::Real rho_r = 0.125;
  amrex::Real u_r = 0.0;
  amrex::Real v_r = 0.0;
  amrex::Real p_r = 0.1;
};
