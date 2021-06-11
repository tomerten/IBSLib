#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "Models.hpp"
#include "NumericFunctions.hpp"
#include "OrdDiffEq.hpp"
#include "RadiationDamping.hpp"
#include "Twiss.hpp"
#include <map>
#include <stdio.h>
#include <string>
#include <vector>

// using namespace std;

void red() { printf("\033[1;31m"); }
void yellow() { printf("\033[1;33m"); }
void green() { printf("\033[1;32m"); }
void blue() { printf("\033[1;34m"); }
void cyan() { printf("\033[1;36m"); }
void reset() { printf("\033[0m"); }

int main() {
  /*
  ================================================================================
  INPUT
  ================================================================================
  */
  string twissfilename = "b2_design_lattice_1996.twiss";
  double aatom = emass / pmass;
  // VeffRFeV
  int nrf = 1;
  double harmon[1];
  double voltages[1];
  harmon[0] = 400.;
  voltages[0] = -4. * 375e3;
  double rho = 4.354479242994256;

  /*
  ================================================================================
  READ FROM INPUT FILES
  ================================================================================
  */
  // read twiss header as map and the table as map of vectors
  map<string, double> twissheadermap;
  map<string, vector<double>> twisstablemap;

  twissheadermap = GetTwissHeader(twissfilename);
  twisstablemap = GetTwissTableAsMap(twissfilename);

  double twissheaderpiwismooth[6] = {
      3.32681701e03, // gamma
      2.40000000e02, // length
      36.96878449,   // gammatr
      1.78499656e01, // q1
      6.74303147e00, // q2
  };

  double twissheadernagaitsev[5] = {
      3.32681701e03, // gamma
      -1.0,          // charge
      2.40000000e02, // length
      1.7,           // energy
      5.10998950e-04 // mass
  };

  double twissheaderrad[6] = {
      3.32681701e03,  // gamma
      1.69999992e00,  // pc
      2.40000000e02,  // length
      5.10998950e-04, // mass
      -1.00000000e00, // charge
      1.78499656e01,  // q1
  };

  double twissode[7] = {
      3.32681701e03,  // gamma
      1.69999992e00,  // pc
      36.96878449,    // gammatr
      5.10998950e-04, // mass
      -1.00000000e00, // charge
      1.78499656e01,  // q1
      2.40000000e02,  // length
  };

  vector<vector<double>> twisstable_piwilattice;
  vector<string> cols /* */ {"L", "BETX", "BETY", "DX"};
  twisstable_piwilattice = GetTable(twissfilename, cols);

  cols.clear();
  vector<vector<double>> twisstable_piwimodified;
  cols.push_back("L");
  cols.push_back("BETX");
  cols.push_back("BETY");
  cols.push_back("DX");
  cols.push_back("DPX");
  cols.push_back("ALFX");
  twisstable_piwimodified = GetTable("b2_design_lattice_1996.twiss", cols);

  cols.clear();
  vector<vector<double>> twisstable_bm;
  cols.push_back("L");
  cols.push_back("BETX");
  cols.push_back("BETY");
  cols.push_back("ALFX");
  cols.push_back("ALFY");
  cols.push_back("DX");
  cols.push_back("DPX");
  cols.push_back("DY");
  cols.push_back("DPY");

  twisstable_bm = GetTable("b2_design_lattice_1996.twiss", cols);

  cols.clear();
  vector<vector<double>> twisstable_rad;
  cols.push_back("L");
  cols.push_back("BETX");
  cols.push_back("BETY");
  cols.push_back("DX");
  cols.push_back("DPX");
  cols.push_back("DY");
  cols.push_back("DPY");
  cols.push_back("ALFX");
  cols.push_back("ALFY");
  cols.push_back("ANGLE");
  cols.push_back("K1L");
  cols.push_back("K1SL");

  twisstable_rad = GetTable("b2_design_lattice_1996.twiss", cols);

  int nrows = static_cast<int>(twisstable_piwilattice.size());
  double twiss_piwilattice[nrows][4];

  // copy to arr
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < twisstable_piwilattice[0].size(); ++j) {
      twiss_piwilattice[i][j] = twisstable_piwilattice[i][j];
      // cout << twissarr[i][j] << " ";
    }
    // cout << endl;
  }

  nrows = static_cast<int>(twisstable_piwimodified.size());
  double twiss_piwimodified[nrows][6];

  // copy to arr
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < twisstable_piwimodified[0].size(); ++j) {
      twiss_piwimodified[i][j] = twisstable_piwimodified[i][j];
      // cout << twissarr[i][j] << " ";
    }
    // cout << endl;
  }

  nrows = static_cast<int>(twisstable_bm.size());
  double twiss_bm[nrows][9];

  // copy to arr
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < twisstable_bm[0].size(); ++j) {
      twiss_bm[i][j] = twisstable_bm[i][j];
      // cout << twiss_bm[i][j] << " ";
    }
    // cout << endl;
  }

  nrows = static_cast<int>(twisstable_rad.size());
  double twiss_rad[nrows][12];

  // copy to arr
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < twisstable_rad[0].size(); ++j) {
      twiss_rad[i][j] = twisstable_rad[i][j];
      // cout << twiss_bm[i][j] << " ";
    }
    // cout << endl;
  }

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

  // updateTwiss
  updateTwiss(twisstablemap);
  // printTwissMap("I2", twisstablemap);
  double betar = BetaRelativisticFromGamma(gammar);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double Lpwd = 1.0e6;
  double neta = eta(gammar, gammatr);
  double sige0 = sigefromsigs(2.0 * pi * 1.2e6, 0.001, 5e-3, gammar, gammatr);
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
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twissheadermap, radint, aatom, omegas);

  /*
   ================================================================================
   COULOMB LOG FUNCTIONS
   ================================================================================
   */
  double pnumber = 1e10;
  double ex = 5e-9;
  double ey = 1e-10;
  double dxavg = twissheadermap["DXRMS"];
  double dyavg = twissheadermap["DYRMS"];
  double charge = twissheadermap["CHARGE"];
  double sige = 1e-4;
  double sigt = 0.005;
  double en0 = twissheadermap["ENERGY"];
  double mass = twissheadermap["MASS"];

  // r0 = erad;
  bool printout = true;

  /*
  ================================================================================
  IBS Models
  ==============================================================================
  */
  printf("\n\n");
  blue();
  printf("ODE Model\n");
  printf("==========\n");
  reset();

  /*

  ================================================================================

  IBS METHODS WITH DEDICATED INTEGRATORS

  ================================================================================
  */

  /*
  ================================================================================
  ODE
  ================================================================================
  */
  vector<double> exa;
  vector<double> eya;
  vector<double> sigsa;
  vector<double> sigea;
  int maxsteps = 20;
  double timestep = 5.0e-3;

  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, timestep, maxsteps,
      exa, eya, sigsa, sigea);
  return 0;
}