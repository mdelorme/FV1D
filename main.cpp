#include <iostream>
#include <array>

#include "SimInfo.h"
#include "Init.h"
#include "ComputeDt.h"
#include "SaveSolution.h"
#include "Update.h"

using namespace fv1d;

int main(int argc, char **argv) {
  read_inifile(argv[1]);

  real_t t = 0.0;
  int ite = 0;

  Array U, Q;
  U.resize(N);
  Q.resize(N);
  
  init(Q);
  primToCons(Q, U);

  real_t dt;
  t = 0.0;
  real_t next_save = 0.0;
  while (t + epsilon < tend) {
    Array Unew{U};

    bool save_needed = (t + epsilon > next_save);

    consToPrim(U, Q);
    dt = compute_dt(Q, (ite == 0 ? save_freq : next_save-t), save_needed);

    if (save_needed) {
      std::cout << " - Saving at time " << t << std::endl;
      save_solution(Q, ite++, t, dt);
      next_save += save_freq;
    }

    update(Q, Unew, dt);
    U = std::move(Unew);

    t += dt;
  }

  save_solution(Q, ite++, t, dt);

  return 0;
}