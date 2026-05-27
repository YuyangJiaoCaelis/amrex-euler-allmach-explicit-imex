#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>

#include "euler_compare/application.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

#if (AMREX_SPACEDIM != 2)
#error "This application requires 2D (AMREX_SPACEDIM == 2)."
#endif

enum ConservedComponent { Rho = 0, Mx = 1, My = 2, E = 3, NCons = 4 };

// Internal numerical modules are included here so the public entry point remains
// thin while the AMReX app keeps a simple GNUmake build.
#include "euler_compare/core/solver_types.hpp"
#include "euler_compare/core/diagnostics.hpp"
#include "euler_compare/core/run_config.hpp"
#include "euler_compare/core/config_reader.hpp"
#include "euler_compare/problems/problem_setup_boundary.hpp"
#include "euler_compare/numerics/explicit_update.hpp"
#include "euler_compare/imex/imex_source_map_helpers.hpp"
#include "euler_compare/imex/imex_pressure_helpers.hpp"
#include "euler_compare/imex/imex_bdltv20_updates.hpp"
#include "euler_compare/io/output_cases.hpp"

}  // namespace

#include "euler_compare/app/driver_main.hpp"
