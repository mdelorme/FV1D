#pragma once

#include "SimInfo.h"

namespace fv1d {

real_t compute_dt(Array &Q, real_t max_dt) {

  real_t inv_dt = 0.0;

  for (int i=ibeg; i<iend; ++i) {
    real_t cs = speed_of_sound(Q[i]);
    inv_dt = std::max(inv_dt, (cs + std::fabs(Q[i][IU]))/dx);
  }

  return std::min(CFL / inv_dt, max_dt);
}

}