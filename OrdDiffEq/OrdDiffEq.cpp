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
         int nrf, double harmon[], double voltages[], double dt, int maxsteps,
         vector<double> &ex, vector<double> &ey, vector<double> &sigs,
         vector<double> sige) {
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
          twiss, radint, aatom, qs);

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
                            gammatr, pc, len, phis, false);
  printf("%-30s %20.6e (%s)\n", "Sige0  :", sige0, "");

  reset();
  sige0 = SigeFromRFAndSigs(sigs[0], U0, charge, nrf, harmon, voltages, gamma,
                            gammatr, pc, len, phis, false);
  sige.push_back(sige0);

  vector<double> sige2;
  sige2.push_back(sige[0] * sige[0]);

  int i = 0;

  vector<double> exx(maxsteps);
  exx[0] = ex[0];
  do {
    i++;
    printf("step : %10i\n", i);
    ex.push_back(ex[0] * exp(-2 * i * dt * (1 / tauradx)) +
                 equi[3] * (1 - exp(-2 * i * dt / tauradx)));
    ey.push_back(ey[0] * exp(-2 * i * dt * (1 / taurady)) +
                 equi[4] * (1 - exp(-2 * i * dt / taurady)));

    sige2.push_back(sige2[0] * exp(-2 * i * dt * (1 / taurads)) +
                    equi[5] * (1 - exp(-2 * i * dt / taurads)));
    sige.push_back(sqrt(sige2[i]));
    sigs.push_back(sigsfromsige(sige[i], gamma, gammatr, omegas));

    printf("--------------------\n");
    printf("step : %10i\n", i);
    cyan();
    printf("ex   : %12.6e\n", ex[i]);
    printf("ey   : %12.6e\n", ey[i]);
    printf("sigs : %12.6e\n", sigs[i]);
    yellow();
    printf("exdiff   %12.6e\n", fabs((ex[i] - ex[i - 1]) / ex[i - 1]));
    printf("eydiff   %12.6e\n", fabs((ey[i] - ey[i - 1]) / ey[i - 1]));
    printf("sigsdiff %12.6e\n", fabs((sigs[i] - sigs[i - 1]) / sigs[i - 1]));
    reset();
  } while (i < maxsteps &&
           (fabs((ex[i] - ex[i - 1]) / ex[i - 1]) > 1e-3 ||
            fabs((ey[i] - ey[i - 1]) / ey[i - 1]) > 1e-3 ||
            fabs((sigs[i] - sigs[i - 1]) / sigs[i - 1]) > 1e-3));
}