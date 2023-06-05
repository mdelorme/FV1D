#pragma once

#include "SimInfo.h"

namespace fv1d {

real_t compute_kappa(real_t x) {
  switch (thermal_conductivity_mode) {
    default: return kappa1;
  }
}

void apply_thermal_conduction(Array &Q, Array &Unew, real_t dt) {
  for (int i=ibeg; i < iend; ++i) {
    real_t x = get_x(i);
    real_t kappaL = compute_kappa(x-0.5*dx);
    real_t kappaR = compute_kappa(x+0.5*dx);

    // Ideal EOS with R = 1 assumed. T = P/rho
    real_t TC = Q[i][IP]   / Q[i][IR];
    real_t TL = Q[i-1][IP] / Q[i-1][IR];
    real_t TR = Q[i+1][IP] / Q[i+1][IR];

    // Computing thermal flux
    real_t FL = kappaL * (TC - TL) / dx;
    real_t FR = kappaR * (TR - TC) / dx;

    /** 
     * Boundaries treatment
     * IMPORTANT NOTE :
     * To be accurate, in the case of fixed temperature, since the temperature is taken at the interface
     * the value of kappa should either be averaged between the cell-centered value and the interface
     * or be evaluated at x=0.25dx / x=xmax-0.25dx
     */
    if (i==ibeg && bctc_xmin != BCTC_NONE) {
      switch (bctc_xmin) {
        case BCTC_FIXED_TEMPERATURE: FL = kappaL * 2.0 * (TC-bctc_xmin_value) / dx; break;
        case BCTC_FIXED_GRADIENT:    FL = kappaL * bctc_xmin_value; break;
        default: break;
      }
    }

    if (i==iend-1 && bctc_xmax != BCTC_NONE) {
      switch (bctc_xmax) {
        case BCTC_FIXED_TEMPERATURE: FR = kappaR * 2.0 * (bctc_xmax_value-TC) / dx; break;
        case BCTC_FIXED_GRADIENT:    FR = kappaR * bctc_xmax_value; break;       
        default: break;
      }
    }

    // And updating using a Godunov-like scheme
    Unew[i][IE] += dt/dx * (FL - FR);
  }
}

}