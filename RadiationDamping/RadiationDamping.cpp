#include "NumericFunctions.hpp"
#include <iostream>
#include <math.h>
#include <stdio.h>

void printradint(double out[6]) {
  printf("\n");
  printf("    Radiation Integrals \n");
  printf("    I2  = %15.6f\n", out[0]);
  printf("    I3  = %15.6f\n", out[1]);
  printf("    I4x = %15.6f\n", out[2]);
  printf("    I4y = %15.6f\n", out[3]);
  printf("    I5x = %15.6f\n", out[5]);
  printf("    I5y = %15.6f\n", out[6]);
}
/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE BUILD C LIBRARY FOR IBS
  AUTHOR    : TOM MERTENS
  DATE      : 04/02/2021
  COPYRIGHT : CERN / HZB

DESCRIPTION :
    CALCULATE APPROXIMATE RADIATION INTEGRALS USING LATTICE AVERAGES

DETAILS:
    average betax = 2 * pi * Qx
    average betay = 2 * pi * Qy

---------------------------------------------------------------------------------------------------------------
latticeLength (double)          : length of the accelerator lattice
gammaTransition (double)        : relativistic transition gamma of the lattice
dipoleBendingRadius (double)    : bending radius of the dipoles (average)
betax (double)                  : average betax
betay (double)                  : average betay
---------------------------------------------------------------------------------------------------------------
radiationIntegrals (double[6])  :
    0 -> I2
    1 -> I3
    2 -> I4x
    3 -> I4y
    4 -> I5x
    5 -> I5y
---------------------------------------------------------------------------------------------------------------

*/
double *RadiationDampingApprox(double latticeLength, double gammaTransition,
                               double dipoleBendingRadius, double betax,
                               double betay) {
  static double radiationIntegrals[6];

  // Courant-Snyder optical functions
  double alphax = 0.0f;
  double alphay = 0.0f;
  double gammax = (1.0f + alphax * alphax) / betax;
  double gammay = (1.0f + alphay * alphay) / betay;

  // average dispersion
  double dx = latticeLength / (2.0f * pi * gammaTransition * gammaTransition);
  double dy = 0.0f;
  double dpx = 0.1f;
  double dpy = 0.0f;

  // curly H
  double hx = betax * dpx + 2.0f * alphax * dx * dpx + gammax * dx;
  double hy = betay * dpy + 2.0f * alphay * dy * dpy + gammay * dy;

  // calculate radiation integrals
  radiationIntegrals[0] = 2.0f * pi / dipoleBendingRadius;
  radiationIntegrals[1] =
      2.0f * pi / (dipoleBendingRadius * dipoleBendingRadius);
  radiationIntegrals[2] = 0.0f;
  radiationIntegrals[3] = 0.0f;
  radiationIntegrals[4] =
      hx * 2.0f * pi / (dipoleBendingRadius * dipoleBendingRadius);
  radiationIntegrals[5] =
      hy * 2.0f * pi / (dipoleBendingRadius * dipoleBendingRadius);

  return radiationIntegrals;
}
/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE BUILD C LIBRARY FOR IBS
  AUTHOR    : TOM MERTENS
  DATE      : 04/02/2021
  COPYRIGHT : CERN / HZB

DESCRIPTION :
    CALCULATE RADIATION INTEGRALS AT EACH LATTICE ELEMENT AND AVERAGE
    OVER THE WHOLE LATTICE

DETAILS:
    TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
        0 -> l
        1 -> betx
        2 -> bety
        3 -> dx
        4 -> dpx
        5 -> dy
        6 -> dpy
        7 -> alfx
        8 -> alfy
        9 -> angle
        10-> k1l
        11-> k1ls
---------------------------------------------------------------------------------------------------------------
twissdata  (double[])       : twiss table data for the required columns (see
DETAILS above) rows (int)                  : number of rows in twissdata cols
(int)                  : number of columns in twissdata (6)
---------------------------------------------------------------------------------------------------------------
radiationIntegrals (double[6])  :
    0 -> I2
    1 -> I3
    2 -> I4x
    3 -> I4y
    4 -> I5x
    5 -> I5y
---------------------------------------------------------------------------------------------------------------
*/

double *RadiationDampingLattice(int rows, double (*twissdata)[12]) {
  static double radiationIntegrals[6];
  double I2 = 0.0;
  double I3 = 0.0;
  double I4x = 0.0;
  double I4y = 0.0;
  double I5x = 0.0;
  double I5y = 0.0;

  double rhoi, ki, gammax, gammay, hx, hy;
  double rhoi2, rhoi3;

#pragma omp parallel for shared(twissdata) reduction(+ : I2, I3, I4x, I4y, I5x,I5y)
  for (int i = 0; i < rows; i++) {
    // local copies
    double *l = &(twissdata[i][0]);
    double *bx = &(twissdata[i][1]);
    double *by = &(twissdata[i][2]);
    double *dx = &(twissdata[i][3]);
    double *dpx = &(twissdata[i][4]);
    double *dy = &(twissdata[i][5]);
    double *dpy = &(twissdata[i][6]);
    double *ax = &(twissdata[i][7]);
    double *ay = &(twissdata[i][8]);
    double *angle = &(twissdata[i][9]);
    double *k1l = &(twissdata[i][10]);
    double *k1sl = &(twissdata[i][11]);

    // calculate local bending radius
    rhoi = (*angle == 0.0) ? 0.0 : *l / *angle;
    rhoi2 = rhoi * rhoi;
    rhoi3 = rhoi2 * rhoi;

    // strength per length unit
    ki = (*l == 0.0) ? 0.0 : *k1l / *l;

    // first for integrals
    I2 += (rhoi == 0.0) ? 0.0 : *l / rhoi2;
    I3 += (rhoi == 0.0) ? 0.0 : *l / rhoi3;
    I4x += (rhoi == 0.0) ? 0.0
                         : (*dx / rhoi3) * *l +
                               (2.0f / rhoi) * (ki * *dx + *k1sl * *dy * *l);
    I4y += 0.0;

    // Courant-Snyder gamma
    gammax = (1.0f + *ax * *ax) / *bx;
    gammay = (1.0f + *ay * *ay) / *by;

    // curly H
    hx = *bx * *dpx * *dpx + 2.0f * *ax * *dx * *dpx + gammax * *dx * *dx;
    hy = *by * *dpy * *dpy + 2.0f * *ay * *dy * *dpy + gammay * *dy * *dy;

    I5x += (rhoi == 0) ? 0.0 : hx * *l / rhoi3;
    I5y += (rhoi == 0) ? 0.0 : hy * *l / rhoi3;
  }

  radiationIntegrals[0] = I2;
  radiationIntegrals[1] = I3;
  radiationIntegrals[2] = I4x;
  radiationIntegrals[3] = I4y;
  radiationIntegrals[4] = I5x;
  radiationIntegrals[5] = I5y;

  return radiationIntegrals;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE BUILD C LIBRARY FOR IBS
  AUTHOR    : TOM MERTENS
  DATE      : 05/02/2021
  COPYRIGHT : CERN / HZB

DESCRIPTION :
    CALCULATE RADIATION GROWTH RATES AND EQUILIBRIUM VALUES

DETAILS:
    TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
        0 -> gamma
        1 -> pc
        2 -> length
        3 -> mass
        4 -> charge
        5 -> qx
    USES RADIATION INTEGRALS
        0 -> I2
        1 -> I3
        2 -> I4x
        3 -> I4y
        4 -> I5x
        5 -> I5y
---------------------------------------------------------------------------------------------------------------
twiss  (double[])               : twiss table data for the required columns (see
DETAILS above) radiationIntegrals (double[6])  : radiation integrals aatom
(double)                  : atomic A - for electrons this is
electron_mass_energy_MeV / proton_mass_energy_MeV
---------------------------------------------------------------------------------------------------------------
output (double[8])
    0 -> alphax
    1 -> alphay
    2 -> alphas
    3 -> exinf
    4 -> eyinf
    5 -> sigEoE2
    6 -> jx
    7 -> jy
---------------------------------------------------------------------------------------------------------------
*/

double *RadiationDampingGrowthRatesAndEquilibriumEmittances(
    double twiss[5], double radiationIntegrals[6], double aatom) {
  const double c = 299792458.0f;
  const double hbar = 1.0545718176461565e-34;
  const double pi = 3.141592653589793f;
  const double electron_volt_joule_relationship = 1.602176634e-19;

  static double output[8];

  double gamma = twiss[0];
  double p0 = twiss[1] * 1.0e9;
  double len = twiss[2] * 1.0;
  double restE = twiss[3] * 1.0e9;
  double charge = twiss[4] * 1.0;

  double particle_radius = charge * charge / aatom * 1.54e-18;

  double CalphaEC = particle_radius * c / (3.0f * restE * restE * restE) *
                    (p0 * p0 * p0 / len);

  // transverse partition numbers
  double jx = 1.0f - radiationIntegrals[3] / radiationIntegrals[0];
  double jy = 1.0f - radiationIntegrals[4] / radiationIntegrals[0];
  double alphax = 2.0f * CalphaEC * radiationIntegrals[0] * jx;
  double alphay = 2.0f * CalphaEC * radiationIntegrals[0] * jy;
  double alphas = 2.0f * CalphaEC * radiationIntegrals[0] * (jx + jy);

  // mc**2 expressed in Joule to match units of cq
  double mass = restE * electron_volt_joule_relationship;
  double cq = 55.0f / (32.0f * sqrt(3)) * (hbar * c) / mass;

  double sigE0E2 = cq * gamma * gamma * radiationIntegrals[1] /
                   (2.0f * radiationIntegrals[0] + radiationIntegrals[3] +
                    radiationIntegrals[4]);
  double exinf =
      cq * gamma * gamma * radiationIntegrals[4] / (jx * radiationIntegrals[0]);
  double eyinf =
      cq * gamma * gamma * radiationIntegrals[5] / (jx * radiationIntegrals[0]);

  double betaAvg = len / (twiss[5] * 2 * pi);

  eyinf = (eyinf == 0.0f) ? cq * betaAvg * radiationIntegrals[1] /
                                (2.0f * jy * radiationIntegrals[0])
                          : eyinf;

  output[0] = 1.0f / alphax;
  output[1] = 1.0f / alphay;
  output[2] = 1.0f / alphas;
  output[3] = exinf;
  output[4] = eyinf;
  output[5] = sigE0E2;
  output[6] = jx;
  output[7] = jy;

  return output;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE BUILD C LIBRARY FOR IBS
  AUTHOR    : TOM MERTENS
  DATE      : 05/02/2021
  COPYRIGHT : CERN / HZB

DESCRIPTION :
    CALCULATE RADIATION LOSSES PER TURN

DETAILS:
    TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
        0 -> gamma
        1 -> pc
        2 -> length
        3 -> mass
        4 -> charge
    USES RADIATION INTEGRALS I2
---------------------------------------------------------------------------------------------------------------
twiss  (double[])               : twiss table data for the required columns (see
DETAILS above) I2 (double)                     : radiation integrals aatom
(double)                  : atomic A - for electrons this is
electron_mass_energy_MeV / proton_mass_energy_MeV
---------------------------------------------------------------------------------------------------------------
*/

double RadiationLossesPerTurn(double twiss[5], double I2, double aatom) {
  const double c = clight;

  double gamma = twiss[0];
  double p0 = twiss[1];
  double len = twiss[2];
  double mass = twiss[3];
  double charge = twiss[4];
  double particle_radius = charge * charge / aatom * 1.54e-18;
  double cgamma = (4.0 * pi / 3.0) * (particle_radius / (mass * mass * mass));
  double betar = sqrt(1.0 - 1.0 / (gamma * gamma));
  double vrev = c * betar;
  double trev = len / vrev;

  return (c * cgamma) / (2.0 * pi * len) * p0 * p0 * p0 * p0 * I2 * 1.0e9 *
         trev;
}

/*
*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : TOM MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 1.0 : ADD RADIATION SPECTRUM CALCULATIONS
    AUTHOR    : TOM MERTENS
    DATE      : 08/02/2021
    COPYRIGHT : HZB

    DESCRIPTION :
        CALCULATE THE CRITICAL FREQUENCY AND ANGLE OF EMITTED RADIATION IN
DIPOLES

    DETAILS:

---------------------------------------------------------------------------------------------------------------
output:
    0 -> omega critical
    1 -> theta critical
    2 -> critical photon energy for given omega
    3 -> average energy per photon
    4 -> average number of photons emitted per turn
---------------------------------------------------------------------------------------------------------------
*/

double *RadiationCriticalEnergy(double rho, double gamma, double omega) {
  const double c = clight;
  const double alphafine = 7.297352569300000e-03;
  const double h = 6.626070150000000e-34;
  const double twoOthree = 2.0f / 3.0f;
  const double gamma3 = gamma * gamma * gamma;

  static double output[5];

  output[0] = twoOthree * c / rho * gamma3;
  output[1] = 1.0f / gamma * pow(output[0] / omega, 1.0f / 3.0f);
  output[2] = h * output[0];
  output[3] = 1.0f / 3.0f * output[2];
  output[4] = 2 * pi * alphafine * gamma;

  return output;
}