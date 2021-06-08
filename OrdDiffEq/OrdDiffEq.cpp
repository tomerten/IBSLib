#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "Models.hpp"
#include "NumericFunctions.hpp"
#include "RadiationDamping.hpp"
#include "Twiss.hpp"
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <vector>
using namespace std;

void ODE(double twiss[], int nrows, double (*twissdata)[12], double harmon[],
         double voltages[], double stepsize, int maxsteps, vector<double> &ex,
         vector<double> &ey, vector<double> &sigs, vector<double> sige) {
  // Radiation integrals
  double gamma = twiss[0];
  double pc = twiss[1];
  double gammatr = twiss[2];
  double mass = twiss[3];
  double charge = twiss[4];
  double q1 = twiss[5];
  double len = twiss[6];
  double twiss_rad[6];
  /*
  for (int i = 0; i < 6; i++) {
    twiss_rad[i] = twiss[i];
  }
  twiss_rad[2] = twiss[6];

  double *radint;
  radint = RadiationDampingLattice(nrows, twissdata);
  printradint(radint);

  double *equi = RadiationDampingGrowthRatesAndEquilibriumEmittances(
      twiss_rad, radint, emass / pmass);

  double tauradx, taurady, taurads;
  tauradx = equi[0];
  taurady = equi[1];
  taurads = equi[2];
  double sigeoe2 = equi[5];

  double U0 = RadiationLossesPerTurn(twiss_rad, radint[0], emass / pmass);

  // Longitudinal
  double phis = SynchronuousPhase(0.0, 173, U0, -1, 1, harmon, voltages, 1e-6);
  vector<double> sige2;
  double betar = sqrt(1 - 1 / (gamma * gamma));
  double trev = len / (betar * clight);
  double frev = 1 / trev;
  double omega = 2 * pi * frev;
  double qs = SynchrotronTune(omega, U0, -1, 1, harmon, voltages, phis,
                              eta(gamma, gammatr), twiss_rad[1]);
  double omegas = omega * qs;
  printf("omega %12.6e\n", omega);
  printf("qs: %12.6e\n", qs);
  printf("omegas: %12.6e\n", omegas);

  sige.push_back(SigeFromRFAndSigs(voltages[0], harmon[0], sigs[0], U0, gamma,
                                   gammatr, pc, len, phis, false));
  printf("Sige : %12.6e\n", sige[0]);
  printf("Sigs : %12.6e\n", sigsfromsige(sige[0], gamma, gammatr, omegas));
  sige2.push_back(sige[0] * sige[0]);

  int steps = 0;
  */
  /*

do {
steps++;
printf("step : %10i", steps);
} while (steps < maxsteps ||
  fabs(ex[steps] - ex[steps - 1] + ey[steps] - ey[steps - 1]) < 1e-3);
  */
}