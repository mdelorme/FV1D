#pragma once

#include <cmath>
#include <functional>
#include "INIReader.h"

namespace fv1d {

using real_t = double;

enum IVar : uint8_t {
  IR = 0,
  IM = 1,
  IU = 1,
  IP = 2,
  IE = 2
};

enum RiemannSolver {
  HLL,
  HLLC
};

enum BoundaryType {
  BC_ABSORBING,
  BC_REFLECTING,
  BC_REFLECTING_GRAVITY,
  BC_HSE
};

enum ReconstructionType {
  PCM,
  PCM_WB,
  PLM
};

// Run
real_t save_freq;
real_t tend;
std::string filename_out = "run.h5";
BoundaryType boundary = BC_REFLECTING;
ReconstructionType reconstruction = PCM; 
RiemannSolver riemann_solver = HLL;
real_t CFL = 0.1;

// Mesh
int Nx;      // Number of domain cells
int Ng;      // Number of ghosts
int N;       // Total number of cells
int ibeg ;   // First cell of the domain
int iend;    // First cell outside of the domain
real_t xmin; // Minimum boundary of the domain
real_t xmax; // Maximum boundary of the domain
real_t dx;   // Space step


// Run and physics
real_t epsilon = 1.0e-6;
real_t gamma0 = 5.0/3.0;
bool gravity = false;
bool no_mechanical_flux_at_bc = true;
bool well_balanced = false;
std::string problem;


// Polytropes and such
real_t m1 = 1.0;
real_t theta1 = 10.0;
real_t g = theta1*(m1+1.0);

// Helper to get the position in the mesh
real_t get_x(int i) {
  return (i-ibeg+0.5) * dx;
}

using State = std::array<real_t, 3>;
using Array = std::vector<State>;

void read_inifile(std::string filename) {
  INIReader reader(filename);

  // Mesh
  Nx = reader.GetInteger("mesh", "Nx", 32);
  Ng = reader.GetInteger("mesh", "Nghosts", 2);
  xmin = reader.GetFloat("mesh", "xmin", 0.0);
  xmax = reader.GetFloat("mesh", "xmax", 1.0);

  N    = Nx + 2*Ng;
  ibeg = Ng;
  iend = Ng+Nx;
  dx = (xmax-xmin) / Nx;
  
  // Run
  tend = reader.GetFloat("run", "tend", 1.0);
  save_freq = reader.GetFloat("run", "save_freq", 1.0e-1);
  filename_out = reader.Get("run", "outpupt_filename", "run.h5");

  std::string tmp;
  tmp = reader.Get("run", "boundaries", "reflecting");
  std::map<std::string, BoundaryType> bc_map{
    {"reflecting_gravity", BC_REFLECTING_GRAVITY},
    {"reflecting",         BC_REFLECTING},
    {"absorbing",          BC_ABSORBING},
    {"hse",                BC_HSE}
  };
  boundary = bc_map[tmp];

  tmp = reader.Get("solvers", "reconstruction", "pcm");
  std::map<std::string, ReconstructionType> recons_map{
    {"pcm",    PCM},
    {"pcm_wb", PCM_WB},
    {"plm",    PLM}
  };
  reconstruction = recons_map[tmp];

  tmp = reader.Get("solvers", "riemann_solver", "hllc");
  std::map<std::string, RiemannSolver> riemann_map{
    {"hll", HLL},
    {"hllc", HLLC}
  };
  riemann_solver = riemann_map[tmp];
  CFL = reader.GetFloat("solvers", "CFL", 0.8);

  // Physics
  epsilon = reader.GetFloat("misc", "epsilon", 1.0e-6);
  gamma0  = reader.GetFloat("physics", "gamma0", 5.0/3.0);
  gravity = reader.GetBoolean("physics", "gravity", false);
  g       = reader.GetFloat("physics", "g", 0.0);
  m1      = reader.GetFloat("polytrope", "m1", 1.0);
  theta1  = reader.GetFloat("polytrope", "theta1", 10.0);
  problem = reader.Get("physics", "problem", "blast");

  no_mechanical_flux_at_bc = reader.GetBoolean("physics", "no_mechanical_flux_at_bc", false);
} 

// Function to apply for the boundaries
std::function<State(Array&, int, int, real_t)> bc_function;
}

// All states operations
#include "States.h"
