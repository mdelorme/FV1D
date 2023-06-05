#pragma once

#include "SimInfo.h"

namespace fv1d {

real_t compute_mu(real_t x) {
  switch (viscosity_mode) {
    default: return mu; break;
  }
}

void apply_viscosity(Array &Q, Array &Unew, real_t dt) {
  for (int i=ibeg; i < iend; ++i) {
    real_t x = get_x(i);
    real_t muL = compute_mu(x-0.5*dx);
    real_t muR = compute_mu(x+0.5*dx);

    // Computing thermal flux
    constexpr real_t four_thirds = 4.0/3.0;
    real_t FL = muL * four_thirds * (Q[i][IU]-Q[i-1][IU])/dx;
    real_t FR = muR * four_thirds * (Q[i+1][IU]-Q[i][IU])/dx;

    // And updating using a Godunov-like scheme
    Unew[i][IE] += dt/dx * (FL - FR);
  }
}

}