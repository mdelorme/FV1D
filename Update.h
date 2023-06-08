#pragma once 

#include "SimInfo.h"
#include "RiemannSolvers.h"
#include "BoundaryConditions.h"
#include "ThermalConduction.h"
#include "Viscosity.h"

namespace fv1d {

void compute_slopes(const Array &Q, Array &slopes) {
  for (int i=1; i < N-1; ++i) {
    for (int ivar=0; ivar < 3; ++ivar) {
      real_t dL = Q[i][ivar] - Q[i-1][ivar];
      real_t dR = Q[i+1][ivar] - Q[i][ivar];

      // Min mod
      if (dL * dR < 0.0)
        slopes[i][ivar] = 0.0;
      else if (fabs(dL) < fabs(dR))
        slopes[i][ivar] = dL;
      else
        slopes[i][ivar] = dR;
    }
  }
}

State reconstruct(State &q, State &slope, real_t sign) {
  State res;
  switch (reconstruction) {
    case PLM: res = q + slope * sign * 0.5; break; // Piecewise Linear
    case PCM_WB: // Piecewise constant + Well-balancing
      res[IR] = q[IR];
      res[IU] = q[IU];
      res[IP] = q[IP] + sign * q[IR] * g * dx * 0.5;
      break;
    default:  res = q; // Piecewise Constant
  }

  return res;
}

void compute_fluxes_and_update(Array &Q, Array &slopes, Array &Unew, real_t dt) {
  #pragma omp parallel for
  for (int i=ibeg; i <= iend; ++i) {
    State qCL = reconstruct(Q[i],   slopes[i], -1.0);
    State qCR = reconstruct(Q[i],   slopes[i],  1.0);
    State qL  = reconstruct(Q[i-1], slopes[i-1],  1.0);
    State qR  = reconstruct(Q[i+1], slopes[i+1], -1.0);

    auto riemann = [&](State qL, State qR, State &flux, real_t &pout) {
      switch (riemann_solver) {
        case HLL: hll(qL, qR, flux, pout); break;
        default: hllc(qL, qR, flux, pout); break;
      }
    };

    // Calculating flux left and right of the cell
    State fluxL, fluxR;
    real_t poutL, poutR;

    riemann(qL, qCL, fluxL, poutL);
    riemann(qCR, qR, fluxR, poutR);

    // Remove mechanical flux in a well-balanced fashion
    if (well_balanced_flux_at_bc && (i==ibeg || i==iend-1)) {
      if (i==ibeg)
        fluxL = State{0.0, poutR - Q[i][IR]*g*dx, 0.0};
      else 
        fluxR = State{0.0, poutL + Q[i][IR]*g*dx, 0.0};
    }

    // Godunov update
    Unew[i] += dt/dx * (fluxL - fluxR);

    if (gravity) {
      // Update momentum
      Unew[i][IM] += dt * Q[i][IR] * g;

      // Update energy
      Unew[i][IE]   += dt * 0.5 * (fluxL[IR] + fluxR[IR]) * g;
    }

    Unew[i][IR] = std::max(1.0e-6, Unew[i][IR]);
  }

}


void update(Array &Q, Array &Unew, real_t dt) {
  // First filling up boundaries for ghosts terms
  fill_boundaries(Q, dt);

  // Hyperbolic update
  Array slopes{(size_t)N};
  compute_slopes(Q, slopes);
  compute_fluxes_and_update(Q, slopes, Unew, dt);

  // Splitted terms
  if (thermal_conductivity_active)
    apply_thermal_conduction(Q, Unew, dt);
  if (viscosity_active)
    apply_viscosity(Q, Unew, dt);
}

}