#pragma once

#include <highfive/H5Easy.hpp>
#include <ostream>
#include <iomanip>

#include "SimInfo.h"

using namespace H5Easy;

namespace fv1d {

void save_solution(const Array &Q, int iteration, real_t t, real_t dt) {
  std::ostringstream oss;
  
  std::setw(4);
  std::setfill('0');
  oss << "ite_" << iteration;
  std::string path = oss.str();
  auto flag = (iteration == 0 ? File::Truncate : File::ReadWrite);

  File file(filename_out, flag);

  if (iteration == 0) {
    file.createAttribute("N", N);
    file.createAttribute("Nx", Nx);
    file.createAttribute("ibeg", ibeg);
    file.createAttribute("iend", iend);
    file.createAttribute("problem", problem);

    std::vector<real_t> x;
    for (int i=ibeg; i < iend; ++i)
      x.push_back(get_x(i));
    file.createDataSet("x", x);
  }

  

  using Vector = std::vector<real_t>;

  Vector vrho, vvel, vP;

  for (int i=ibeg; i<iend; ++i) {
    real_t rho = Q[i][IR];
    real_t vel = Q[i][IU];
    real_t Ek = 0.5 * rho * vel*vel;
    real_t p = Q[i][IP];

    vrho.push_back(rho);
    vvel.push_back(vel);
    vP.push_back(p);
  }

  auto group = file.createGroup(path);
  group.createDataSet("rho", vrho);
  group.createDataSet("vel", vvel);
  group.createDataSet("prs", vP);
  group.createAttribute("time", t);
}

}