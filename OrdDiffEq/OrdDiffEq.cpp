#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "Models.hpp"
#include "NumericFunctions.hpp"
#include "RadiationDamping.hpp"
#include "Twiss.hpp"
#include <algorithm>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

void red() { printf("\033[1;31m"); }
void yellow() { printf("\033[1;33m"); }
void green() { printf("\033[1;32m"); }
void blue() { printf("\033[1;34m"); }
void cyan() { printf("\033[1;36m"); }
void reset() { printf("\033[0m"); }

void ODE(map<string, double> &twiss, map<string, vector<double>> &twissdata,
         int nrf, double harmon[], double voltages[], double stepsize,
         int maxsteps, vector<double> &ex, vector<double> &ey,
         vector<double> &sigs, vector<double> sige) {
  // Radiation integrals
  double gamma = twiss["GAMMA"];
  double pc = twiss["PC"];
  double gammatr = twiss["GAMMATR"];
  double mass = twiss["MASS"];
  double charge = twiss["CHARGE"];
  double q1 = twiss["Q1"];
  double len = twiss["LENGTH"];
  double twiss_rad[6];

  double aatom = emass / pmass;
  double betar = BetaRelativisticFromGamma(gamma);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double neta = eta(gamma, gammatr);
  double epsilon = 1.0e-6;

  double *radint;
  radint = RadiationDampingLattice(twissdata);
  printradint(radint);

  double U0 = RadiationLossesPerTurn(twiss, radint[1], aatom);
  double phis =
      SynchronuousPhase(0.0, 173, U0, charge, nrf, harmon, voltages, epsilon);
  double qs =
      SynchrotronTune(omega, U0, charge, nrf, harmon, voltages, phis, neta, pc);
  double omegas = qs * omega;

  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twiss, radint, aatom, omegas);

  double tauradx, taurady, taurads;
  tauradx = equi[0];
  taurady = equi[1];
  taurads = equi[2];
  double sigeoe2 = equi[5];

  cyan();
  printf("%-30s %20.6e (%s)\n", "Tx :", tauradx, "");
  printf("%-30s %20.6e (%s)\n", "Ty :", taurady, "");
  printf("%-30s %20.6e (%s)\n", "Ts :", taurads, "");
  blue();
  printf("%-30s %20.6e (%s)\n", "qs      :", qs, "");
  printf("%-30s %20.6e (%s)\n", "synch freq :", omegas, "");
  printf("%-30s %20.6e (%s)\n", "SigEOE2 :", sigeoe2, "");
  printf("%-30s %20.6e (%s)\n", "SigEOE  :", sqrt(sigeoe2), "");
  printf("%-30s %20.6e (%s)\n", "eta  :", eta(gamma, gammatr), "");

  printf("%-30s %20.6e (%s)\n", "Sigs  :", sigs[0], "");
  printf("%-30s %20.6e (%s)\n", "Sigsinf  :", equi[6], "");
  double sige0 = sigefromsigs(omega, equi[6], qs, gamma, gammatr);
  printf("%-30s %20.6e (%s)\n", "Sige0  :", sige0, "");
  sige0 = SigeFromRFAndSigs(equi[6], U0, charge, nrf, harmon, voltages, gamma,
                            gammatr, pc, len, phis, true);
  printf("%-30s %20.6e (%s)\n", "Sige0  :", sige0, "");

  reset();

  /*

    double U0 = RadiationLossesPerTurn(twiss_rad, radint[0], emass / pmass);

    // Longitudinal
    double phis = SynchronuousPhase(0.0, 173, U0, -1, 1, harmon, voltages,
    1e-6); vector<double> sige2; double betar = sqrt(1 - 1 / (gamma * gamma));
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