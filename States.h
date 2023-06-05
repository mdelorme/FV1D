#pragma once

namespace fv1d {

State primToCons(State &q) {
  State res;
  res[IR] = q[IR];
  res[IM] = q[IR]*q[IU];

  real_t Ek = 0.5 * res[IM] * q[IU];
  res[IE] = (Ek + q[IP] / (gamma0-1.0));
  return res;
}

State consToPrim(State &u) {
  State res;
  res[IR] = u[IR];
  res[IU] = u[IM] / u[IR];

  real_t Ek = 0.5 * res[IU] * u[IM];
  res[IP] = (u[IE] - Ek) * (gamma0-1.0);
  return res; 
}

void consToPrim(Array &U, Array &Q) {
  for (int i=0; i < N; ++i)
    Q[i] = consToPrim(U[i]);
}

void primToCons(Array &Q, Array &U) {
  for (int i=0; i < N; ++i)
    U[i] = primToCons(Q[i]);
}

real_t speed_of_sound(State &q) {
  return std::sqrt(q[IP] * gamma0 / q[IR]);
}


State& operator+=(State &a, State b) {
  for (int i=0; i < 3; ++i)
    a[i] += b[i];
  return a;
}

State& operator-=(State &a, State b) {
  for (int i=0; i < 3; ++i)
    a[i] -= b[i];
  return a;
}

State operator*(const State &a, real_t q) {
  State res;
  for (int i=0; i < 3; ++i)
    res[i] = a[i]*q;
  return res;
}

State operator/(const State &a, real_t q) {
  State res;
  for (int i=0; i < 3; ++i)
    res[i] = a[i]/q;
  return res;
}

State operator*(real_t q, const State &a) {
  return a*q;
}

State operator+(const State &a, const State &b) {
  State res;
  for (int i=0; i < 3; ++i)
    res[i] = a[i]+b[i];
  return res;
}

State operator-(const State &a, const State &b) {
  State res;
  for (int i=0; i < 3; ++i)
    res[i] = a[i]-b[i];
  return res;
}

State compute_flux(State &q) {
  const real_t Ek = 0.5 * q[IR] * q[IU] * q[IU];
  const real_t E = (q[IP] / (gamma0-1.0) + Ek);

  State fout {
    q[IR]*q[IU],
    q[IR]*q[IU]*q[IU] + q[IP],
    (q[IP] + E) * q[IU]
  };

  return fout;
}
}