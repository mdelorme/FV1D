#pragma once

#include "SimInfo.h"
#include "BoundaryConditions.h"


namespace fv1d {

void init_sod(Array &Q) {
  std::cout << "Initializing SOD" << std::endl;
  
  for (int i=ibeg; i < iend; ++i) {
    if (get_x(i) <= 0.5) {
      Q[i][IR] = 1.0;
      Q[i][IP] = 1.0;
      Q[i][IU] = 0.0;
    }
    else {
      Q[i][IR] = 0.125;
      Q[i][IP] = 0.1;
      Q[i][IU] = 0.0;
    }
  }
}

void init_blast(Array &Q) {
  std::cout << "Initializing blast" << std::endl;

  real_t xmid = 0.5 * (xmin+xmax);

  for (int i=ibeg; i < iend; ++i) {
    real_t x = get_x(i);
    real_t r = fabs(xmid - x);

    std::cout << "xmid = " << xmid << "; x = " << x << "; r = " << r << std::endl;

    if (r < 0.1) {
      Q[i][IR] = 1.0;
      Q[i][IU] = 0.0;
      Q[i][IP] = 10.0;
    }
    else {
      Q[i][IR] = 1.2;
      Q[i][IU] = 0.0;
      Q[i][IP] = 0.1;
    }
  }
}

void init_C91(Array &Q) {
  std::cout << "Initializing C91" << std::endl;

  for (int i=ibeg; i < iend; ++i) {
    real_t x = get_x(i);
    real_t T = (1.0 + theta1*x);
    real_t rho = std::pow(T, m1);
    real_t prs = std::pow(T, m1+1.0);

    Q[i][IR] = rho;
    Q[i][IU] = 0.0;
    Q[i][IP] = prs;
  }
}

void init_C91_hse(Array &Q) {
  std::cout << "Initializing C91 - HSE" << std::endl;

  for (int i=ibeg; i < iend; ++i) {
    real_t x = get_x(i);
    real_t T = (1.0 + theta1*x);
    real_t rho = (1.0/(theta1*(m1+1))*(  std::pow(1.0+theta1*(x+0.5*dx), m1+1) 
                                       - std::pow(1.0+theta1*(x-0.5*dx), m1+1)));
    
    real_t prs = (1.0/(theta1*(m1+2))*(  std::pow(1.0+theta1*(x+0.5*dx), m1+2) 
                                       - std::pow(1.0+theta1*(x-0.5*dx), m1+2)));

    Q[i][IR] = rho;
    Q[i][IU] = 0.0;
    Q[i][IP] = prs;
  }
}

void init(Array &Q) {
  if (problem == "SOD")
    init_sod(Q);
  else if (problem == "blast")
    init_blast(Q);
  else if (problem == "C91")
    init_C91(Q);
  else if (problem == "C91_hse")
    init_C91_hse(Q);

  // One time init stuff
  init_boundaries();
}


}