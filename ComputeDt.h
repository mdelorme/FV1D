#pragma once

#include "SimInfo.h"
#include "ThermalConduction.h"
#include "Viscosity.h"

namespace fv1d {

real_t compute_dt(Array &Q, real_t max_dt, bool diag) {

  real_t inv_dt_hyp  = 0.0;
  real_t inv_dt_tc   = 0.0;
  real_t inv_dt_visc = 0.0;

  #pragma omp parallel reduction(max:inv_dt_hyp) reduction(max:inv_dt_tc) reduction(max:inv_dt_visc)
  for (int i=ibeg; i<iend; ++i) {
    real_t x = get_x(i);

    // Hyperbolic CFL
    real_t cs = speed_of_sound(Q[i]);
    inv_dt_hyp = std::max(inv_dt_hyp, (cs + std::fabs(Q[i][IU]))/dx);
  
    // Parabolic
    if (thermal_conductivity_active)
      inv_dt_tc = std::max(inv_dt_tc, 2.0 * compute_kappa(x) / (dx*dx));
    
    if (viscosity_active)
      inv_dt_visc = std::max(inv_dt_visc, 2.0 * compute_mu(x) / (dx*dx));
  }

  if (diag) {
    std::cout << "Computing timesteps : dt_hyp=" << 1.0/inv_dt_hyp 
              << "; dt_TC="   << 1.0/inv_dt_tc 
              << "; dt_visc=" << 1.0/inv_dt_visc << std::endl; 
  }

  real_t dt = CFL / (std::max({inv_dt_hyp, inv_dt_tc, inv_dt_visc}));

  return std::min(dt, max_dt);
}

}