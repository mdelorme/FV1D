#pragma once

#include <fstream>

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
    real_t rho = (1.0/(dx*theta1*(m1+1))*(  std::pow(1.0+theta1*(x+0.5*dx), m1+1) 
                                          - std::pow(1.0+theta1*(x-0.5*dx), m1+1)));
    
    real_t prs = (1.0/(dx*theta1*(m1+2))*(  std::pow(1.0+theta1*(x+0.5*dx), m1+2) 
                                          - std::pow(1.0+theta1*(x-0.5*dx), m1+2)));

    Q[i][IR] = rho;
    Q[i][IU] = 0.0;
    Q[i][IP] = prs;
  }
}

void init_B02(Array &Q) {
  std::cout << "Initializing B02" << std::endl;

  for (int i=ibeg; i < iend; ++i) {
    real_t x = get_x(i);

    real_t T0   = 1.0;
    real_t rho0 = std::pow(T0, m1);
    real_t p0   = std::pow(T0, m1+1);

    real_t T1   = T0 + theta1*b02_xmid;
    real_t rho1 = std::pow(T1, m1);
    real_t p1   = std::pow(T1, m1+1);

    if (x < b02_xmid) {
      real_t T = T0 + x*theta1;
      Q[i][IR] = rho0 * std::pow(T, m1);
      Q[i][IP] = p0 * std::pow(T, m1+1);
    }
    else {
      real_t T = T1 + (x-b02_xmid)*theta2;
      Q[i][IR] = rho1 * std::pow(T/T1, m2);
      Q[i][IP] = p1 * std::pow(T/T1, m2+1);
    }

    Q[i][IU] = 0.0;
  }
}

void init_B02_restart(Array &Q) {
  std::cout << "Restarting B02 from file ics.dat" << std::endl;
  
  std::string filename_in = "ics.dat";
  std::ifstream f_in(filename_in);

  // Skipping first line
  std::string tmp;
  std::getline(f_in, tmp);
 
  for (int i=ibeg; i < iend; ++i) {
    f_in >> tmp >> Q[i][IR] >> Q[i][IU] >> Q[i][IP];
    constexpr real_t pert_amp = 1.0e-3;
    real_t pert = ((float)rand() / (float)RAND_MAX - 0.5)*2.0;
    real_t x = get_x(i);
    const real_t pert_margin = dx;
    if (x < pert_margin || x > (b02_xmid-pert_margin))
      pert = 0.0;

    Q[i][IP] *= (1.0 + pert);
  }
  
  f_in.close();
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
  else if (problem == "B02")
    init_B02(Q);
  else if (problem == "B02_restart")
    init_B02_restart(Q);

  // One time init stuff
  init_boundaries();
}


}