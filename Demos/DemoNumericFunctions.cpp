#include "NumericFunctions.hpp"
#include "RadiationDamping.hpp"
#include "Twiss.hpp"
#include <map>
#include <math.h>
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
  CHECK CONSTANTS
  ================================================================================
  */
  blue();
  printf("CONSTANTS\n");
  printf("=========\n");
  reset();
  printf("%-30s %20.6e (%s)\n", "Speed of light :", clight, "m/s");
  printf("%-30s %20.6e (%s)\n", "Reduced Planck constant :", hbar, "GeV");
  printf("%-30s %20.6e (%s)\n", "Electron mass :", emass, "GeV");
  printf("%-30s %20.6e (%s)\n", "Proton mass :", pmass, "GeV");
  printf("%-30s %20.6e (%s)\n", "Neutron mass :", nmass, "GeV");
  printf("%-30s %20.6e (%s)\n", "Atomic Mass Unit :", mumass, "GeV");
  printf("%-30s %20.6e (%s)\n", "Pi :", pi, "rad");
  printf("%-30s %20.6e (%s)\n", "Electric Charge :", ec, "C");
  printf("%-30s %20.6e (%s)\n", "Classical electron radius :", erad, "m");
  printf("%-30s %20.6e (%s)\n\n", "Converted proton radius :", prad, "m");

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
  double neta = eta(gammar, gammatr);
  double *radint;
  radint = RadiationDampingLattice(twisstablemap);
  double U0 = RadiationLossesPerTurn(twissheadermap, radint[1], emass / pmass);
  double phis =
      SynchronuousPhase(0.0, 172, U0, -1, 1, harmon, voltages, epsilon);
  double VrfEffeV =
      EffectiveRFVoltageInElectronVolt(phis, -1, 1, harmon, voltages);
  double VrfEffeVp =
      EffectiveRFVoltageInElectronVoltPrime(phis, -1, 1, harmon, voltages);
  double VrfEffeVU0 = VeffRFeVRadlosses(phis, U0, -1, 1, harmon, voltages);
  double qs =
      SynchrotronTune(omega, U0, -1, 1, harmon, voltages, phis, neta, pc);
  double sige0 = sigefromsigs(omega, 0.005, qs, gammar, gammatr);
  double omegas = qs * omega;
  double VrfEffeVPWD = VeffRFeVPotentialWellDistortion(
      173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc);
  double VrfEffeVPWDp = VeffRFeVPotentialWellDistortionPrime(
      173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc);
  double phisPWD = SynchronuousPhaseWithPWD(
      0.0, 173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc, epsilon);
  double qsPWD = SynchrotronTunePWD(omega, U0, -1, 1, harmon, voltages, Lpwd, 1,
                                    0.005, phisPWD, neta, pc);

  blue();
  printf("Basic Functions\n");
  printf("===============\n");
  reset();

  printf("%-30s %20.6e (%s)\n", "SigEfromSigs :", sige0, "");
  printf("%-30s %20.6e (%s)\n", "Eta :", neta, "");
  printf("%-30s %20.6e (%s)\n", "Fmohl          :", fmohl(a, b, q, npp), "");
  printf("%-30s %20.6e (%s)\n", "Particle Radius :", r0, "m");
  printf("%-30s %20.6e (%s)\n", "Relativistic beta from gamma :", betar, "");
  printf("%-30s %20.6e (%s)\n", "Rds (Nagaitsev) :", rds(1, 2, 3), "");
  printf("%-30s %20.6e (%s)\n", "phis :", phis, "deg");
  printf("%-30s %20.6e (%s)\n", "VeffRFeV :", VrfEffeV, "eV");
  printf("%-30s %20.6e (%s)\n", "VeffRFeVPrime:", VrfEffeVp, "eV");
  printf("%-30s %20.6e (%s)\n", "U0:", U0, "eV");
  printf("%-30s %20.6e (%s)\n", "VeffRFeVRadlosses :", VrfEffeVU0, "eV");
  printf("%-30s %20.6e (%s)\n", "synchronuousphase :", phis, "rad");
  printf("%-30s %20.6e (%s)\n", "synchrotronTune :", qs, "");
  printf("%-30s %20.6e (%s)\n", "synchrotron angularfrequency :", qs * omega,
         "Hz");
  printf("%-30s %20.6e (%s)\n", "VeffRFeVPWD :", VrfEffeVPWD, "eV");
  printf("%-30s %20.6e (%s)\n", "VeffRFeVPWDprime:", VrfEffeVPWDp, "");
  printf("%-30s %20.6e (%s)\n", "synchronuousphasewithPWD :", phisPWD, "");
  printf("%-30s %20.6e (%s)\n", "synchrotronTunePWD :", qsPWD, "");

  return 0;
}