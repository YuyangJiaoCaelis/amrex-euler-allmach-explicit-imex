// Problem setup and finite-volume helper include hub.
// The implementation is split by numerical role to keep the AMReX driver readable.
#include "euler_compare/core/state_primitives.hpp"
#include "problem_state_functions.hpp"
#include "euler_compare/numerics/reconstruction_fluxes.hpp"
#include "exact_riemann.hpp"
#include "problem_initialization.hpp"
#include "euler_compare/numerics/diagnostic_metrics_timestep.hpp"
#include "boundary_conditions.hpp"
#include "advection_blob_transport.hpp"
