#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "Models.hpp"
#include "NumericFunctions.hpp"
#include "RadiationDamping.hpp"
#include "Twiss.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
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

void WriteToFile(string filename, vector<double> &ex, vector<double> &ey,
                 vector<double> &sigs) {
  ofstream csvfile(filename);
  if (csvfile.is_open()) {
    int num_of_rows = min({ex.size(), ey.size(), sigs.size()});
    csvfile << "ex,ey,sigs" << endl;
    for (int i = 0; i < num_of_rows; i++) {
      csvfile << ex[i] << "," << ey[i] << "," << sigs[i] << endl;
    }
  } else {
    cout << "File could not be opened";
  }
  csvfile.close();
}

void ODE(map<string, double> &twiss, map<string, vector<double>> &twissdata,
         int nrf, double harmon[], double voltages[], double dt, int maxsteps,
         vector<double> &ex, vector<double> &ey, vector<double> &sigs,
         vector<double> sige, int model, double pnumber) {
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

  double *ibs;
  double aes, aex, aey;
  vector<double> extemp, eytemp, sige2temp;
  extemp.push_back(ex[0]);
  eytemp.push_back(ey[0]);
  sige2temp.push_back(sige2[0]);

  switch (model) {
  case 1:
    ibs = PiwinskiSmooth(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, r0);
    break;
  case 2:
    ibs = PiwinskiLattice(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                          twissdata, r0);
    break;
  case 3:
    ibs = PiwinskiLatticeModified(pnumber, ex[0], ey[0], sigs[0], sige[0],
                                  twiss, twissdata, r0);
    break;

  case 4:
    ibs = Nagaitsev(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata,
                    r0);
    break;
  }

  printouts(ibs);
  do {
    switch (model) {
    case 1:
      ibs = PiwinskiSmooth(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 2:
      ibs = PiwinskiLattice(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                            twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 3:
      ibs = PiwinskiLatticeModified(pnumber, ex[0], ey[0], sigs[0], sige[0],
                                    twiss, twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 4:
      ibs = Nagaitsev(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata,
                      r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    }

    // double dex, dey, dsige2;
    cyan();
    printf("Previous Values \n");
    printf("================\n");

    printf("ex   : %12.6e\n", ex[i]);
    printf("ey   : %12.6e\n", ey[i]);
    printf("sigs : %12.6e\n\n", sigs[i]);
    reset();

    i++;
    printf("IBS amplitude growth rates\n");
    printf("==========================\n");
    printf("step : %10i\n", i);
    printf("aes : %10.6f\n", aes);
    printf("aex : %10.6f\n", aex);
    printf("aey : %10.6f\n\n", aey);

    ex.push_back(extemp[i - 1] * exp(2 * dt * (-1 / tauradx + aex)) +
                 equi[3] * (1 - exp(-2 * dt * i / tauradx)));
    extemp.push_back(extemp[i - 1] * exp(2 * dt * (-1 / tauradx + aex)));

    ey.push_back(eytemp[i - 1] * exp(2 * dt * (-1 / taurady + aey)) +
                 equi[4] * (1 - exp(-2 * dt * i / taurady)));
    eytemp.push_back(eytemp[i - 1] * exp(2 * dt * (-1 / taurady + aey)));

    sige2.push_back(sige2temp[i - 1] * exp(2 * dt * (-1 / taurads + aes)) +
                    equi[5] * (1 - exp(-2 * i * dt / taurads)));
    sige2temp.push_back(sige2temp[i - 1] * exp(2 * dt * (-1 / taurads + aes)));

    sige.push_back(sqrt(sige2[i]));
    sigs.push_back(sigsfromsige(sige[i], gamma, gammatr, omegas));

    printf("--------------------\n");
    printf("step : %10i\n", i);
    cyan();
    printf("ex   : %12.6e\n", ex[i]);
    printf("ey   : %12.6e\n", ey[i]);
    printf("sigs : %12.6e\n", sigs[i]);
    green();
    printf("exdiff abs  : %12.6e\n", (ex[i] - ex[i - 1]));
    printf("eydiff abs  : %12.6e\n", (ey[i] - ey[i - 1]));
    printf("sigsdiff abs: %12.6e\n", (sigs[i] - sigs[i - 1]));
    yellow();
    printf("exdiff   : %12.6e\n", fabs((ex[i] - ex[i - 1]) / ex[i - 1]));
    printf("eydiff   : %12.6e\n", fabs((ey[i] - ey[i - 1]) / ey[i - 1]));
    printf("sigsdiff : %12.6e\n", fabs((sigs[i] - sigs[i - 1]) / sigs[i - 1]));
    reset();
  } while (i < maxsteps &&
           (fabs((ex[i] - ex[i - 1]) / ex[i - 1]) > 1e-3 ||
            fabs((ey[i] - ey[i - 1]) / ey[i - 1]) > 1e-3 ||
            fabs((sigs[i] - sigs[i - 1]) / sigs[i - 1]) > 1e-3));
  blue();
  printf("%20s %12.6e\n", "Final ex :", ex[ex.size() - 1]);
  printf("%20s %12.6e\n", "Final ey :", ey[ey.size() - 1]);
  printf("%20s %12.6e\n", "Final sigs :", sigs[sigs.size() - 1]);
  reset();
}