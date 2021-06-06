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
  // read twiss header as map and the table as map of vectors
  string twissfilename = "b2_design_lattice_1996.twiss";
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
  CHECK CONSTANTS
  ================================================================================
  */
  blue();
  printf("CONSTANTS\n");
  printf("=========\n");
  reset();
  printf("%-30s %20.6f (%s)\n", "Speed of light :", clight, "m/s");
  printf("%-30s %20.6e (%s)\n", "Reduced Planck constant :", hbar, "GeV");
  printf("%-30s %20.6e (%s)\n", "Electron mass :", emass, "GeV");
  printf("%-30s %20.6e (%s)\n", "Proton mass :", pmass, "GeV");
  printf("%-30s %20.6e (%s)\n", "Neutron mass :", nmass, "GeV");
  printf("%-30s %20.6e (%s)\n", "Atomic Mass Unit :", mumass, "GeV");
  printf("%-30s %20.6e (%s)\n", "Pi :", pi, "");
  printf("%-30s %20.6e (%s)\n", "Electric Charge :", ec, "C");
  printf("%-30s %20.6e (%s)\n", "Classical electron radius :", erad, "m");
  printf("%-30s %20.6e (%s)\n\n", "Converted proton radius :", prad, "m");

  /*
  ================================================================================
  BASIC NUMERIC FUNCTIONS
  ================================================================================
  */

  blue();
  printf("Basic Functions\n");
  printf("===============\n");
  reset();

  // sigefromsigs
  printf("%-30s %20.6e (%s)\n",
         "SigEfromSigs :", sigefromsigs(2.0 * pi * 5e5, 0.001, 5e3, 10.0), "");

  // eta
  printf("%-30s %20.6e (%s)\n", "Eta :", eta(3600.0, 37.0), "");

  // fmohl
  double a = 5.709563671168914e-04;
  double b = 2.329156389696222e-01;
  double q = 2.272866910079534e00;
  int npp = 1000;
  printf("%-30s %20.6e (%s)\n", "Fmohl          :", fmohl(a, b, q, npp), "");

  // particle radius
  printf("%-30s %20.6e (%s)\n", "Paricle Radius :", particle_radius(1, 1), "m");

  // beta relativistic from gamma
  printf("%-30s %20.6e (%s)\n",
         "Relativistic beta from gamma :", BetaRelativisticFromGamma(1), "m/s");

  // rds from nagaitsev paper
  printf("%-30s %20.6e (%s)\n", "Rds (Nagaitsev) :", rds(1, 2, 3), "");

  // VeffRFeV
  double harmon[1];
  double voltages[1];
  harmon[0] = 400.;
  voltages[0] = -4. * 350e3;
  printf("%-30s %20.6e (%s)\n",
         "VeffRFeV :", VeffRFeV(173, -1, 1, harmon, voltages), "eV");

  // VeffRFeVprime
  printf("%-30s %20.6e (%s)\n",
         "VeffRFeVPrime:", VeffRFeVPrime(173, -1, 1, harmon, voltages), "eV");

  // updateTwiss
  updateTwiss(twisstablemap);
  // printTwissMap("I2", twisstablemap);

  // Energy loss per turn
  double *radint;
  radint = RadiationDampingLattice(nrows, twiss_rad);
  double U0 = RadiationLossesPerTurn(twissheaderrad, radint[0], emass / pmass);
  printf("%-30s %20.6e (%s)\n", "VeffRFeVRadlosses :",
         VeffRFeVRadlosses(173, U0, -1, 1, harmon, voltages), "eV");
  printf("%-30s %20.6e (%s)\n", "synchronuousphase :",
         synchronuousphase(0.0, 173, U0, -1, 1, harmon, voltages, 1e-3), "");

  double betar =
      sqrt(1 - 1 / (twissheaderpiwismooth[0] * twissheaderpiwismooth[0]));
  double trev = twissheaderpiwismooth[1] / (betar * clight);
  double frev = 1 / trev;
  double omega = 2 * pi * frev;
  double etaa = eta(twissheaderpiwismooth[0], twissheaderpiwismooth[2]);
  double Lpwd = 1e6;
  printf("%-30s %20.6e (%s)\n", "synchrotronTune :",
         synchrotronTune(
             omega, U0, -1, 1, harmon, voltages,
             synchronuousphase(0.0, 173, U0, -1, 1, harmon, voltages, 1e-3),
             etaa, twissheaderrad[1]),
         "");

  printf("%-30s %20.6e (%s)\n", "VeffRFeVPotentialWellDistortion :",
         VeffRFeVPotentialWellDistortion(173, U0, -1, 1, harmon, voltages, Lpwd,
                                         1, 0.005, twissheaderrad[1]),
         "eV");
  printf("%-30s %20.6e (%s)\n", "VeffRFeVPotentialWellDistortionPWD :",
         VeffRFeVPotentialWellDistortionPrime(173, U0, -1, 1, harmon, voltages,
                                              Lpwd, 1, 0.005,
                                              twissheaderrad[1]),
         "");
  printf("%-30s %20.6e (%s)\n", "synchronuousphasewithPWD :",
         synchronuousphasewithPWD(0.0, 173, U0, -1, 1, harmon, voltages, Lpwd,
                                  1, 0.005, twissheaderrad[1], 1e-3),
         "");
  printf("%-30s %20.6e (%s)\n", "synchrotronTunePWD :",
         synchrotronTunePWD(
             omega, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005,
             synchronuousphasewithPWD(0.0, 173, U0, -1, 1, harmon, voltages,
                                      Lpwd, 1, 0.005, twissheaderrad[1], 1e-3),
             etaa, twissheaderrad[1]),
         "");
  /*
  ================================================================================
  RADIATION DAMPING METHODS
  ================================================================================
  */
  yellow();
  printf("Radiation Damping\n");
  printf("=================\n");
  green();
  printf("Radiation Smooth Ring Approximation\n");
  reset();
  radint = RadiationDampingApprox(twissheadernagaitsev[2],
                                  twissheaderpiwismooth[2], 4.35, 2.13, 5.66);
  printradint(radint);

  green();
  printf("Radiation Damping element by element\n");
  reset();
  radint = RadiationDampingLattice(nrows, twiss_rad);
  printradint(radint);

  green();
  printf("Radiation Damping Equib \n");
  reset();
  double *equi = RadiationDampingGrowthRatesAndEquilibriumEmittances(
      twissheaderrad, radint, emass / pmass);
  printf("%-30s %10.6e (%s)\n", "Taux :", equi[0], "s");
  printf("%-30s %10.6e (%s)\n", "Tauy :", equi[1], "s");
  printf("%-30s %10.6e (%s)\n", "Taus :", equi[2], "s");
  printf("%-30s %10.6e (%s)\n", "exinf :", equi[3], "");
  printf("%-30s %10.6e (%s)\n", "eyinf :", equi[4], "");
  printf("%-30s %10.6e (%s)\n", "sgeoe2 :", equi[5], "");
  printf("%-30s %10.6e (%s)\n", "jx :", equi[6], "");
  printf("%-30s %10.6e (%s)\n", "jy :", equi[7], "");

  green();
  printf("Radiation Losses per turn. \n");
  reset();
  printf("%-30s %10.6e (%s)\n", "U0 :",
         RadiationLossesPerTurn(twissheaderrad, radint[0], emass / pmass),
         "eV/Turn");

  green();
  printf("Critical Energy Calculations. \n");
  reset();
  double *critical;
  critical = RadiationCriticalEnergy(4.35, twissheaderrad[0], 2 * pi * 1.25e6);
  printf("%-30s %10.6e (%s)\n", "omega_crit :", critical[0], "");
  printf("%-30s %10.6e (%s)\n", "Theta_crit :", critical[1], "");
  printf("%-30s %10.6e (%s)\n", "E_crit :", critical[2], "eV");
  printf("%-30s %10.6e (%s)\n", "E_per_photon_avg :", critical[3], "eV/photon");
  printf("%-30s %10.6e (%s)\n", "N_photons_avg_per_turn :", critical[4],
         "1/turn");
  /*
   ================================================================================
   COULOMB LOG FUNCTIONS
   ================================================================================
   */
  blue();
  printf("CoulombLog Functions\n");
  printf("====================\n");
  reset();
  double clog[2];
  twclog(1e10, 7, 12, 1, 0, 5e-9, 1e-10, 1e-18, 3600, 1, 1700, 0.5, 0.005, 2e-9,
         clog);
  printf("%-30s %20.6e (%s)\n", "CoulombLog value :", clog[0], "");
  printf("%-30s %20.6e (%s)\n\n", "CoulombLog constant :", clog[1], "");

  double pnumber = 1e10;
  double ex = 5e-9;
  double ey = 1e-10;
  double twissheader[9] = {
      3.32681701e03,  // gamma
      -1.00000000e00, // charge
      2.40000000e02,  // length
      1.7,            // energy
      5.10998950e-04, // mass
      1.78499656e01,  // q1
      6.74303147e00,  // q2
      2.24396462e-01, // dxrms
      0.00000000e00,  // dyrms
  };
  double sige = 1e-4;
  double sigt = 0.005;
  double r0 = erad;
  bool printout = true;
  CoulombLog(1e10, equi[3], equi[4], twissheader, sige, sigt, r0, printout,
             clog);

  yellow();
  printf("\nWith Tailcut...\n");
  reset();

  TailCutCoulombLog(1e10, equi[3], equi[4], twissheader, sige, sigt, equi[0],
                    equi[1], equi[2], r0, printout, clog);

  cyan();
  printf("\nRing averaged Coulomblog...\n");
  reset();
  double avgclog = 0.0;
  for (int i; i < nrows; ++i) {
    twclog(pnumber, twiss_rad[i][1], twiss_rad[i][2], twiss_rad[i][3],
           twiss_rad[i][4], equi[3], equi[4], r0, twissheader[0],
           twissheader[1], twissheader[3], twissheader[4], sige, sigt, clog);
    avgclog += twiss_rad[i][0] * clog[0];
  }
  avgclog /= twissheader[2];
  printf("CoulombLog avg : %12.6e", avgclog);

  cyan();
  printf("\nRing averaged Coulomblog with Tailcut...\n");
  reset();
  for (int i; i < nrows; ++i) {
    twclogtail(pnumber, twiss_rad[i][0], twiss_rad[i][1], twiss_rad[i][2],
               twiss_rad[i][3], twiss_rad[i][4], twiss_rad[i][5],
               twiss_rad[i][5], twiss_rad[i][7], twiss_rad[i][8],
               twiss_rad[i][9], twiss_rad[i][10], twiss_rad[i][11], equi[3],
               equi[4], r0, emass / pmass, twissheader[0], twissheader[3],
               twissheader[2], twissheader[4], twissheader[1], sige, sigt,
               clog);
    avgclog += twiss_rad[i][0] * clog[0];
  }
  avgclog /= twissheader[2];
  printf("CoulombLog avg : %12.6e", avgclog);
  /*
   ================================================================================
   INTEGRATOR FUNCTIONS
   ================================================================================
   */
  printf("\n\n");
  red();
  printf("Integrator Functions\n");
  printf("====================\n");
  reset();
  double tau[3];
  a = 1, b = 2;
  double c = 3;
  double c1 = 4, cx = 5, cy = 6, cprime = 7;
  double cyy = 8, tl1 = 9, tl2 = 10, tx1 = 11, tx2 = 12, ty1 = 13, ty2 = 14;
  SimpsonDecade(a, b, c, c1, cx, cy, cprime, cyy, tl1, tl2, tx1, tx2, ty1, ty2,
                tau);
  printf("SimpsonDecade madx...\n");
  printf("%-30s %20.6e (%s)\n", "al :", tau[0], "");
  printf("%-30s %20.6e (%s)\n", "ax :", tau[1], "");
  printf("%-30s %20.6e (%s)\n", "ay :", tau[2], "");
  printf("\n\n");

  printf("%-30s %20.6e (%s)\n",
         "Variable integrand :", IBSIntegralIntegrand(1, 2, 3, 4, 5, 6), "");
  printf("%-30s %20.6e (%s)\n", "Standard Simpson integral :",
         simpson(IBSIntegralIntegrand, 1, 2, 3, 4, 5, 6, 7, 10), "");

  printf("\n\n");
  printf("Standard SimpsonDecade ...\n");
  intSimpson(IBSIntegralIntegrand, 1, 2, 3, 4, 6, 7, 8, 9, 10, tau);
  printf("%-30s %20.6e (%s)\n", "al :", tau[0], "");
  printf("%-30s %20.6e (%s)\n", "ax :", tau[1], "");
  printf("%-30s %20.6e (%s)\n", "ay :", tau[2], "");
  printf("\n\n");

  /*
  ================================================================================
  IBS Models
  ==============================================================================
  */
  printf("\n\n");
  green();
  printf("IBS Models\n");
  printf("==========\n");
  cyan();
  printf("Piwinski Smooth ...\n");
  reset();
  double *res;
  res = PiwinskiSmooth(1e10, equi[3], equi[4], sigt, sige,
                       twissheaderpiwismooth, r0);
  printouts(res);

  cyan();
  printf("Piwinski Lattice ...\n");
  reset();

  res = PiwinskiLattice(1e10, equi[3], equi[4], sigt, sige,
                        twissheaderpiwismooth, 610, twiss_piwilattice, r0);
  printouts(res);

  cyan();
  printf("Piwinski Lattice modified...\n");
  reset();

  res = PiwinskiLatticeModified(1e10, equi[3], equi[4], sigt, sige,
                                twissheaderpiwismooth, 610, twiss_piwimodified,
                                r0);
  printouts(res);
  cyan();
  printf("Nagaitsev...\n");
  reset();

  res = Nagaitsev(1e10, equi[3], equi[4], sigt, sige, twissheadernagaitsev,
                  nrows, twiss_piwimodified, r0);
  printouts(res);

  cyan();
  printf("Nagaitsev... with tailcut\n");
  reset();

  res =
      Nagaitsevtailcut(1e10, equi[3], equi[4], sigt, sige, twissheadernagaitsev,
                       nrows, twiss_rad, r0, emass / pmass);
  printouts(res);
  /*

  ================================================================================

  IBS METHODS WITH DEDICATED INTEGRATORS

  ================================================================================
  */
  green();
  printf("IBS Madx...\n");
  red();
  reset();

  res = ibsmadx(1e10, equi[3], equi[4], sigt, sige, twissheadernagaitsev, nrows,
                twiss_bm, r0, true);
  green();
  printf("IBS Madx... with tailcut\n");
  red();
  reset();

  res = ibsmadxtailcut(1e10, equi[3], equi[4], sigt, sige, twissheadernagaitsev,
                       nrows, twiss_rad, r0, emass / pmass);
  printouts(res);

  yellow();
  printf("IBS Bjorken-Mtingwa... Failing\n");
  red();
  printf("Growth rates have no meaning - using standards Simpson integration "
         "does not work !!!.\n");
  reset();
  res = BjorkenMtingwa2(1e10, equi[3], equi[4], sigt, sige,
                        twissheadernagaitsev, nrows, twiss_piwimodified, r0);
  printouts(res);

  green();
  printf("IBS Bjorken-Mtingwa... now using simpson per decade integration \n");
  reset();
  res = BjorkenMtingwa(1e10, equi[3], equi[4], sigt, sige, twissheadernagaitsev,
                       nrows, twiss_bm, r0);
  printouts(res);

  green();
  printf("IBS Bjorken-Mtingwa... now using simpson per decade integration with "
         "tailcut \n");
  reset();
  res = BjorkenMtingwatailcut(1e10, equi[3], equi[4], sigt, sige,
                              twissheadernagaitsev, nrows, twiss_rad, r0,
                              emass / pmass);
  printouts(res);

  green();
  printf("IBS Conte-Martini...  using simpson per decade integration \n");
  reset();
  res = ConteMartini(1e10, equi[3], equi[4], sigt, sige, twissheadernagaitsev,
                     nrows, twiss_bm, r0);
  printouts(res);

  green();
  printf("IBS Conte-Martini...  using simpson per decade integration with "
         "tailcut\n");
  reset();
  res = ConteMartinitailcut(1e10, equi[3], equi[4], sigt, sige,
                            twissheadernagaitsev, nrows, twiss_rad, r0,
                            emass / pmass);
  printouts(res);

  green();
  printf("IBS MADX...  using simpson per decade integration \n");
  printf("IBS MADX...  uses a slightly different way of calculating the "
         "integrand than ibsmadx.\n");
  reset();
  res = MadxIBS(1e10, equi[3], equi[4], sigt, sige, twissheadernagaitsev, nrows,
                twiss_bm, r0);
  printouts(res);

  /*
  ================================================================================
  ODE
  ================================================================================
  */
  vector<double> exa;
  vector<double> eya;
  vector<double> sigsa;
  vector<double> sigea;
  int maxsteps = 10;
  double stepsize = 1.0;

  exa.push_back(equi[3]);
  eya.push_back(equi[4]);
  sigsa.push_back(sigt);

  ODE(twissode, nrows, twiss_rad, harmon, voltages, stepsize, maxsteps, exa,
      eya, sigsa, sigea);
  return 0;
}