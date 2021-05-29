#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "NumericFunctions.hpp"
#include <stdio.h>

int main() {
  /*
  ================================================================================
  CHECK CONSTANTS
  ================================================================================
  */
  printf("CONSTANTS\n");
  printf("=========\n");
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
  printf("Basic Functions\n");
  printf("===============\n");

  printf("%-30s %20.6e (%s)\n",
         "SigEfromSigs :", sigefromsigs(2.0 * pi * 5e5, 0.001, 5e3, 10.0), "");
  printf("%-30s %20.6e (%s)\n", "Eta :", eta(3600.0, 37.0), "");

  double a = 5.709563671168914e-04;
  double b = 2.329156389696222e-01;
  double q = 2.272866910079534e00;
  int npp = 1000;
  printf("%-30s %20.6e (%s)\n", "Fmohl :", fmohl(a, b, q, npp), "");
  printf("%-30s %20.6e (%s)\n", "Paricle Radius :", particle_radius(1, 1), "m");
  printf("%-30s %20.6e (%s)\n",
         "Relativistic beta from gamma :", BetaRelativisticFromGamma(1), "m/s");
  printf("%-30s %20.6e (%s)\n\n", "Rds (Nagaitsev) :", rds(1, 2, 3), "");

  /*
   ================================================================================
   COULOMB LOG FUNCTIONS
   ================================================================================
   */
  printf("CoulombLog Functions\n");
  printf("====================\n");
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
      1.69999992e00,  // pc
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
  CoulombLog(1e10, ex, ey, twissheader, sige, sigt, r0, printout, clog);

  printf("\nWith Tailcut...\n");
  TailCutCoulombLog(1e10, ex, ey, twissheader, sige, sigt, 0.005, 0.008, 0.01,
                    r0, printout, clog);
  /*
   ================================================================================
   COULOMB LOG FUNCTIONS
   ================================================================================
   */
  printf("\n\n");
  printf("Integrator Functions\n");
  printf("====================\n");
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
  printf("Standar SimpsonDecade ...\n");
  intSimpson(IBSIntegralIntegrand, 1, 2, 3, 4, 6, 7, 8, 9, 10, tau);
  printf("%-30s %20.6e (%s)\n", "al :", tau[0], "");
  printf("%-30s %20.6e (%s)\n", "ax :", tau[1], "");
  printf("%-30s %20.6e (%s)\n", "ay :", tau[2], "");
  printf("\n\n");
  return 0;
}