#pragma once

#include <map>
#include <cassert>

#include "SimInfo.h"

namespace fv1d {

State fill_absorbing(Array &Q, int i, int iref, real_t dt) {
  return Q[iref];
}

State fill_reflecting(Array &Q, int i, int iref, real_t dt) {
  int ipiv = (i < iref ? ibeg : iend);
  int isym = 2*ipiv - i - 1;
  State q = Q[isym];
  q[IU] *= -1.0;

  return q;
}

State fill_reflecting_gravity(Array &Q, int i, int iref, real_t dt) {
  int ipiv = (i < iref ? ibeg : iend);
  int isym = 2*ipiv - i - 1;
  
  State q = Q[iref];

  q[IU] = -Q[isym][IU];

  if (gravity)
    q[IU] += dt*g;

  return q;  
}


State fill_hse(Array &Q, int i, int iref, real_t dt) {
  real_t x = get_x(i);
  
  int isym = (i < iref ? iref+1 : iref-1);
  State qsym = Q[isym];
  real_t xsym = get_x(isym);

  real_t rho = std::pow(1.0 + x*theta1, m1);
  real_t sign = (i < iref ? -1.0 : 1.0);
  real_t prs = qsym[IP] + 0.5 * dx * g * (rho + qsym[IR]);

  return State{rho, 0.0, prs};
}


/**
 * One time initialization for the boundary conditions
 **/
void init_boundaries() {
  assert(Ng == 2); // Most functions have been coded with Ng = 2 !!!!

  switch (boundary) {
    case BC_HSE: bc_function = fill_hse; break;
    case BC_ABSORBING: bc_function = fill_absorbing; break;
    case BC_REFLECTING_GRAVITY: bc_function = fill_reflecting_gravity; break;
    default: bc_function = fill_reflecting; break;
  }
}

/**
 * Fills the boundary with a value. The value is returned from a user-defed method.
 * Note that this function is overly complicated for this type of code.
 * There is a real reason to that: to replicate the way it is called in dyablo
 */
void fill_boundaries(Array &Q, real_t dt) {
  for (int i=0; i < Ng; ++i) {
    int itop = i;
    int ibot = iend+i;
    int iref_top = ibeg;
    int iref_bot = iend-1;

    // And apply
    Q[itop] = bc_function(Q, itop, iref_top, dt);
    Q[ibot] = bc_function(Q, ibot, iref_bot, dt);
  }
}


}