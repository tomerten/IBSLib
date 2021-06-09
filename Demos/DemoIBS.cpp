#include "CoulombLogFunctions.hpp"
#include "Models.hpp"
#include "NumericFunctions.hpp"
#include "RadiationDamping.hpp"
#include "Twiss.hpp"
#include <map>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

void red() { printf("\033[1;31m"); }
void yellow() { printf("\033[1;33m"); }
void green() { printf("\033[1;32m"); }
void blue() { printf("\033[1;34m"); }
void cyan() { printf("\033[1;36m"); }
void reset() { printf("\033[0m"); }

int main() {
  string twissfilename = "b2_design_lattice_1996.twiss";
  map<string, double> twissheadermap;
  map<string, vector<double>> twisstablemap;

  twissheadermap = GetTwissHeader(twissfilename);
  twisstablemap = GetTwissTableAsMap(twissfilename);

  double aatom = emass / pmass;
  // VeffRFeV
  int nrf = 1;
  double harmon[1];
  double voltages[1];
  harmon[0] = 400.;
  voltages[0] = -4. * 375e3;

  /*
  ================================================================================
  BASIC NUMERIC FUNCTIONS
  ================================================================================
  */
  // parameters needed in this section
  double gammar = twissheadermap["GAMMA"];
  double gammatr = twissheadermap["GAMMATR"];
  double pc = twissheadermap["PC"];
  double len = twissheadermap["LENGTH"];
  // fmohl
  double a = 5.709563671168914e-04;
  double b = 2.329156389696222e-01;
  double q = 2.272866910079534e00;
  int npp = 1000;
  // synchronuous phase eps
  double epsilon = 1.0e-6;
  double rho = 4.35;

  // updateTwiss
  updateTwiss(twisstablemap, rho);
  // printTwissMap("I2", twisstablemap);
  double betar = BetaRelativisticFromGamma(gammar);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double Lpwd = 1.0e6;
  double sigs = 0.005;
  double neta = eta(gammar, gammatr);
  double VrfEffeV =
      EffectiveRFVoltageInElectronVolt(173, -1, 1, harmon, voltages);
  double VrfEffeVp =
      EffectiveRFVoltageInElectronVoltPrime(173, -1, 1, harmon, voltages);
  double *radint;
  radint = RadiationDampingLattice(twisstablemap);
  double U0 = RadiationLossesPerTurn(twissheadermap, radint[1], emass / pmass);
  double VrfEffeVU0 = VeffRFeVRadlosses(173, U0, -1, 1, harmon, voltages);
  double phis =
      SynchronuousPhase(0.0, 173, U0, -1, 1, harmon, voltages, epsilon);
  double qs =
      SynchrotronTune(omega, U0, -1, 1, harmon, voltages, phis, neta, pc);
  double sige0 = sigefromsigs(omega, sigs, qs, gammar, gammatr);
  double omegas = qs * omega;
  double VrfEffeVPWD = VeffRFeVPotentialWellDistortion(
      173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc);
  double VrfEffeVPWDp = VeffRFeVPotentialWellDistortionPrime(
      173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc);
  double phisPWD = SynchronuousPhaseWithPWD(
      0.0, 173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc, epsilon);
  double qsPWD = SynchrotronTunePWD(omega, U0, -1, 1, harmon, voltages, Lpwd, 1,
                                    0.005, phisPWD, neta, pc);
  double bxavg = len / (2.0 * pi * twissheadermap["Q1"]);
  double byavg = len / (2.0 * pi * twissheadermap["Q2"]);

  // va
  double pnumber = 1e10;
  double ex = 5e-9;
  double ey = 1e-10;
  double dxavg = twissheadermap["DXRMS"];
  double dyavg = twissheadermap["DYRMS"];
  double charge = twissheadermap["CHARGE"];
  double en0 = twissheadermap["ENERGY"];
  double mass = twissheadermap["MASS"];
  double sigt = sigs / clight;
  bool printout = true;
  radint = RadiationDampingLattice(twisstablemap);
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twissheadermap, radint, aatom, omegas);

  printf("\n\n");
  blue();
  printf("IBS Models\n");
  printf("==========\n");
  green();
  printf("\nPiwinski Smooth ...\n");
  reset();

  double *res;
  sige0 = 1e-4;
  res = PiwinskiSmooth(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                       r0);
  printouts(res);

  green();
  printf("\nPiwinski Lattice ...\n");
  reset();

  res = PiwinskiLattice(1e10, equi[3], equi[4], sigs, sige0, twissheadermap,
                        twisstablemap, r0);
  printouts(res);

  green();
  printf("\nPiwinski Lattice modified...\n");
  reset();

  res = PiwinskiLatticeModified(pnumber, equi[3], equi[4], sigs, sige0,
                                twissheadermap, twisstablemap, r0);
  printouts(res);

  green();
  printf("\nNagaitsev...\n");
  reset();

  res = Nagaitsev(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                  twisstablemap, r0);
  printouts(res);

  cyan();
  printf("\nNagaitsev... with tailcut\n");
  reset();

  res = Nagaitsevtailcut(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                         twisstablemap, r0, aatom);
  printouts(res);
  /*

  ================================================================================

  IBS METHODS WITH DEDICATED INTEGRATORS

  ================================================================================
  */
  green();
  printf("\nIBS Madx...\n");
  red();
  reset();

  res = ibsmadx(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                twisstablemap, r0, true);

  cyan();
  printf("\nIBS Madx... with tailcut\n");
  red();
  reset();

  res = ibsmadxtailcut(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                       twisstablemap, r0, aatom);
  printouts(res);

  red();
  printf("\nIBS Bjorken-Mtingwa... Failing\n");
  red();
  printf("Vertical growth rates failed using standards Simpson integration "
         "!!!.\n");
  reset();
  res = BjorkenMtingwa2(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                        twisstablemap, r0);
  printouts(res);

  green();
  printf(
      "\nIBS Bjorken-Mtingwa... now using simpson per decade integration \n");
  reset();
  res = BjorkenMtingwa(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                       twisstablemap, r0);
  printouts(res);

  cyan();
  printf(
      "\nIBS Bjorken-Mtingwa... now using simpson per decade integration with "
      "tailcut \n");
  reset();
  res = BjorkenMtingwatailcut(pnumber, equi[3], equi[4], sigs, sige0,
                              twissheadermap, twisstablemap, r0, aatom);
  printouts(res);

  green();
  printf("\nIBS Conte-Martini...  using simpson per decade integration \n");
  reset();
  res = ConteMartini(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                     twisstablemap, r0);
  printouts(res);

  cyan();
  printf("\nIBS Conte-Martini...  using simpson per decade integration with "
         "tailcut\n");
  reset();
  res = ConteMartinitailcut(pnumber, equi[3], equi[4], sigs, sige0,
                            twissheadermap, twisstablemap, r0, aatom);
  printouts(res);

  green();
  printf("\nIBS MADX...  using simpson per decade integration \n");
  printf("IBS MADX...  uses a slightly different way of calculating the "
         "integrand than ibsmadx.\n");
  reset();
  res = MadxIBS(pnumber, equi[3], equi[4], sigs, sige0, twissheadermap,
                twisstablemap, r0);
  printouts(res);
  return 0;
}