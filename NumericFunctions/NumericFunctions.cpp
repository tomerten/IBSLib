#include "NumericFunctions.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;
/*
================================================================================
================================================================================
METHOD TO CALCULATE dE/E FROM sigma_s.

REMARK:
-------
  dE/E (do not confuse with dp/p, there is a factor beta**2) De/E = beta**2 dp/p

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 16/02/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - double : qs
        longitunal tune
    - double : sigs
        bunch length [m]
    - double : omega0
        angular RF frequency
    - double : eta
        phase slip fact


  Returns:
  --------
    double
      dE/E

================================================================================
================================================================================
*/
double sigefromsigs(double omega0, double sigs, double qs, double eta) {
  return qs * omega0 * (sigs / (fabs(eta) * clight));
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE sigma_s FROM dE/E.

REMARK:
-------
  dE/E (do not confuse with dp/p, there is a factor beta**2) De/E = beta**2 dp/p

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - double : sige
        dE/E
    - double : gamma
        relativistic gamma
    - double : gammatr
        relativistic gamma for transition
    - double : omegas
        synchrotron angular frequency (omega0 * qs)


  Returns:
  --------
    double
      sigma_s

================================================================================
================================================================================
*/

double sigsfromsige(double sige, double gamma, double gammatr, double omegas) {
  return clight * fabs(eta(gamma, gammatr)) / omegas * sige;
}
/*
================================================================================
================================================================================
METHOD TO CALCULATE PHASE SLIP FACTOR ETA.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - double : gamma
        relativistic gamma
    - double : gammatr
        relativistic gamma for transition

  Returns:
  --------
    double
      eta - phase slip factor (1/ gammatr**2  - 1/ gamma**2)

================================================================================
================================================================================
*/

double eta(double gamma, double gammatr) {
  return 1.0 / (gammatr * gammatr - 1.0 / (gamma * gamma));
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE FMOHL.

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    -  MICHAELA SCHAUMANN

  HISTORY:
    - 05/02/2021 COPYRIGHT : CERN / HZB
    - 08/06/2021 : initial cpp version (Tom)

  REFS:
    - HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS
    - BANE: A SIMPLIFIED MODEL OF INTRABEAM SCATTERING (EPAC 2002)
    - PIWINSKI  Tech. Rep. HEAC 74, Stanford, 1974.

================================================================================
  Arguments:
  ----------
    - double : a
        parameter 1
    - double : b
        parameter 2
    - double : q
        parameter 3
    - int : n

  Returns:
  --------
    double
      Fmohl

================================================================================
================================================================================
*/
double fmohl(double a, double b, double q, int n) {

  double u, cp, cq;
  double sum = 0.0;
  double const du = 1.0 / n; // careful needs to be 1.0 not 1

#pragma omp parallel for reduction(+ : sum)
  for (int i = 0; i <= n; i++) {
    double dsum = 0.0;
    u = i * du;
    cp = sqrt(a * a + (1.0 - a * a) * u * u);
    cq = sqrt(b * b + (1.0 - b * b) * u * u);
    dsum = 2.0 * log(q * (1.0 / cp + 1.0 / cq) / 2.0) - euler;
    dsum *= (1.0 - 3.0 * u * u) / (cp * cq);
    dsum = (i == 0) ? dsum /= 2.0 : dsum;
    dsum = (i == n) ? dsum /= 2.0 : dsum;
    sum += dsum;
  }

  sum *= 8 * pi * du;

  return sum;
}

/*
-------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
--------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES
(P-PB)
AUTHOR    : TOM MERTENS
DATE      : 05/02/2021
COPYRIGHT : CERN / HZB

    DESCRIPTION :
        CALCULATE CLASSICAL PARTICLE RADIUS

    REF:
        HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS (p. )

---------------------------------------------------------------------------------
charge (double)		: particle charge
aatom				: atomic number - for electrons this is
electron_mass_energy_MeV / proton_mass_energy_MeV
---------------------------------------------------------------------------------
*/
double particle_radius(double charge, double aatom) {
  return charge * charge / aatom * prad;
}

/*
---------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
---------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES
(P-PB)
AUTHOR    : TOM MERTENS
DATE      : 05/02/2021
COPYRIGHT : CERN / HZB

    DESCRIPTION :
        CALCULATE BETA RELATIVISTIC FROM GAMMA RELATIVISTC

---------------------------------------------------------------------------------
gamma (double)  : relativistic gamma
---------------------------------------------------------------------------------
*/
double BetaRelativisticFromGamma(double gamma) {
  return sqrt(1 - (1 / (gamma * gamma)));
}
/*
REF : https://stackoverflow.com/questions/3437404/min-and-max-in-c

#define GENERIC_MAX(x, y) ((x) > (y) ? (x) : (y))

#define ENSURE_int(i)   _Generic((i), int:   (i))
#define ENSURE_float(f) _Generic((f), float: (f))


#define MAX(type, x, y) \
  (type)GENERIC_MAX(ENSURE_##type(x), ENSURE_##type(y))
*/

#define max(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a > _b ? _a : _b;                                                         \
  })

#define min(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a < _b ? _a : _b;                                                         \
  })

/*
---------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
---------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES
(P-PB)
AUTHOR    : TOM MERTENS
DATE      : 06/02/2021
COPYRIGHT : CERN

    DESCRIPTION :
        FUNCTIONS TO CALCULATE NAGAITSEV FUNCTION  RDS

    REF:
        PRSTAB 8, 064403 (2005)
---------------------------------------------------------------------------------
x (double)  : rds parameter
y (double)  : rds parameter
z (double)  : rds parameter
---------------------------------------------------------------------------------
*/
double rds(double x, double y, double z) {
  // init
  double errtol = 0.05;
  // double tiny   = 1.0e-25;
  // double big    = 4.5e21;
  double c1 = 3.0 / 14.0;
  double c2 = 1.0 / 6.0;
  double c3 = 9.0 / 22.0;
  double c4 = 3.0 / 26.0;
  double c5 = 0.25 * c3;
  double c6 = 1.50 * c4;

  double xt = x;
  double yt = y;
  double zt = z;
  double sum = 0.0;
  double fac = 1.0;
  int iter = 0;

  double sqrtx, sqrty, sqrtz, alamb;
  double ave, delx, dely, delz;

  do {
    iter = iter + 1;
    sqrtx = sqrt(xt);
    sqrty = sqrt(yt);
    sqrtz = sqrt(zt);

    alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
    sum = sum + fac / (sqrtz * (zt + alamb));

    fac = 0.25 * fac;
    xt = 0.25 * (xt + alamb);
    yt = 0.25 * (yt + alamb);
    zt = 0.25 * (zt + alamb);
    ave = 0.20 * (xt + yt + 3.0 * zt);
    delx = (ave - xt) / ave;
    dely = (ave - yt) / ave;
    delz = (ave - zt) / ave;
  } while (max(max(fabs(delx), fabs(dely)), fabs(delz)) >= errtol);

  double ea = delx * dely;
  double eb = delz * delz;
  double ec = ea - eb;
  double ed = ea - 6.0 * eb;
  double ee = ed + ec + ec;

  return 3.0 * sum +
         fac *
             (1.0 + ed * (-c1 + c5 * ed - c6 * delz * ee) +
              delz * (c2 * ee + delz * (-c3 * ec + delz * c4 * ea))) /
             (ave * sqrt(ave));
}

double VeffRFeV(double phi, double charge, int nrf, double harmon[],
                double voltages[]) {
  double vrf = voltages[0] * sin(phi);
  for (int i = 1; i < nrf; i++) {
    vrf += voltages[i] * sin((harmon[i] / harmon[0]) * phi);
  }
  vrf *= charge;

  return vrf;
}

double VeffRFeVPrime(double phi, double charge, int nrf, double harmon[],
                     double voltages[]) {
  double vrf = voltages[0] * cos(phi);
  for (int i = 1; i < nrf; i++) {
    vrf += voltages[i] * (harmon[i] / harmon[0]) *
           cos((harmon[i] / harmon[0]) * phi);
  }
  vrf *= charge;

  return vrf;
}

double VeffRFeVRadlosses(double phi, double U0, double charge, int nrf,
                         double harmon[], double voltages[]) {
  double vrf = VeffRFeV(phi, charge, nrf, harmon, voltages) - U0;
  return vrf;
}

/*
================================================================================

REF :
https://www.quantstart.com/articles/Implied-Volatility-in-C-using-Template-Functions-and-Newton-Raphson/

================================================================================
*/

double synchronuousphase(double target, double init_phi, double U0,
                         double charge, int nrf, double harmon[],
                         double voltages[], double epsilon) {
  // Set the initial option prices and volatility
  double y = VeffRFeVRadlosses(init_phi, U0, charge, nrf, harmon, voltages);
  double x = init_phi;

  while (fabs(y - target) > epsilon) {
    double d_x = VeffRFeVPrime(x, charge, nrf, harmon, voltages);
    x += (target - y) / d_x;
    y = VeffRFeVRadlosses(x, U0, charge, nrf, harmon, voltages);
  }
  return x;
}

double VeffRFeVPotentialWellDistortion(double phi, double U0, double charge,
                                       int nrf, double harmon[],
                                       double voltages[], double L, double N,
                                       double sigs, double pc) {

  return VeffRFeVRadlosses(phi, U0, charge, nrf, harmon, voltages) +
         charge * L * N * ec * clight * clight /
             (sqrt(2 * pi) * sigs * sigs * sigs * pc * 1e9);
}

double VeffRFeVPotentialWellDistortionPrime(double phi, double U0,
                                            double charge, int nrf,
                                            double harmon[], double voltages[],
                                            double L, double N, double sigs,
                                            double pc) {
  return VeffRFeVPrime(phi, charge, nrf, harmon, voltages) +
         charge * ec * L * N * clight * clight /
             (sqrt(2.0 * pi) * sigs * sigs * sigs * pc * 1e9);
}

double synchronuousphasewithPWD(double target, double init_phi, double U0,
                                double charge, int nrf, double harmon[],
                                double voltages[], double L, double N,
                                double sigs, double pc, double epsilon) {
  // Set the initial option prices and volatility
  double y = VeffRFeVPotentialWellDistortion(init_phi, U0, charge, nrf, harmon,
                                             voltages, L, N, sigs, pc);
  double x = init_phi;

  while (fabs(y - target) > epsilon) {
    double d_x = VeffRFeVPotentialWellDistortionPrime(
        x, U0, charge, nrf, harmon, voltages, L, N, sigs, pc);
    x += (target - y) / d_x;
    y = VeffRFeVPotentialWellDistortion(x, U0, charge, nrf, harmon, voltages, L,
                                        N, sigs, pc);
  }
  return x;
}

double synchrotronTune(double omega0, double U0, double charge, int nrf,
                       double harmon[], double voltages[], double phis,
                       double eta, double pc) {
  return sqrt(harmon[0] * eta *
              fabs(VeffRFeVPrime(phis, charge, nrf, harmon, voltages)) /
              (2 * pi * pc * 1e9));
}

double synchrotronTunePWD(double omega0, double U0, double charge, int nrf,
                          double harmon[], double voltages[], double L,
                          double N, double sigs, double phis, double eta,
                          double pc) {
  return sqrt(harmon[0] * eta *
              fabs(VeffRFeVPotentialWellDistortionPrime(
                  phis, U0, charge, nrf, harmon, voltages, L, N, sigs, pc)) /
              (2 * pi * pc * 1e9));
}

double csige(double v0, double h0, double sigs, double U0, double gamma,
             double gammatr, double pc, double circ, double phis,
             bool printout) {

  double betar = sqrt(1 - 1 / (gamma * gamma));
  double trev = circ / (betar * clight);
  double frev = 1 / trev;
  double eta = (1.0 / (gammatr * gammatr)) - (1.0 / (gamma * gamma));
  double omega0 = 2 * pi / trev;

  double nus =
      sqrt(fabs(h0 * eta) / (2 * pi * betar * pc) * fabs(v0 * cos(phis)));
  double sige = nus * omega0 * sigs / (clight * eta);
  sige = 2 * pi * nus * sigs / (eta * clight);

  if (printout) {
    printf("Synchrotron Tune : %12.6e\n ", nus);
    printf("Synchrotron freq : %12.6e\n ", nus * omega0);
    printf("Sige             : %12.6e\n ", sige);
  }
  return sige;
}

/*
================================================================================
================================================================================
METHOD TO UPDATE THE MADX TWISS TABLE WITH COLUMNS (KEYS IN THE MAP) WITH
RECURRING VALUES IN THE IBS AND RADIATION DAMPING ROUTINES.

REMARK:
-------
  The code below is not optimized for speed but for readability.
  In principle this function is run once per simulation so the
  time gain of optimization would be minimal.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 06/08/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - map<string, <double>>& : table
      Twiss table

  Returns:
  --------
    void
      updates the Twiss table (as map) with new keys

================================================================================
================================================================================
 */
void updateTwiss(map<string, vector<double>> &table) {
  // get length of table to reserve the vector sizes
  int size = table["L"].size();

  vector<double> rho(size), k(size), I2(size), I3(size), I4x(size), I4y(size),
      I5x(size), I5y(size), gammax(size), gammay(size), hx(size), hy(size);

  // calculate the new columns
  for (int i = 0; i < size; i++) {
    double angle = table["ANGLE"][i];
    double l = table["L"][i];
    double bx = table["BETX"][i];
    double by = table["BETY"][i];
    double ax = table["ALFX"][i];
    double ay = table["ALFY"][i];
    double dx = table["DX"][i];
    double dpx = table["DPX"][i];
    double dy = table["DY"][i];
    double dpy = table["DPY"][i];
    double k1l = table["K1L"][i];
    double k1sl = table["K1SL"][i];
    double rhoi2, rhoi3;

    // local bending radius
    rho[i] = (angle == 0.0) ? 0.0 : l / angle;
    rhoi2 = rho[i] * rho[i];
    rhoi3 = rhoi2 * rho[i];

    k[i] = (l == 0.0) ? 0.0 : k1l / l;

    // first for integrals
    I2[i] = (rho[i] == 0.0) ? 0.0 : l / rhoi2;
    I3[i] = (rho[i] == 0.0) ? 0.0 : l / rhoi3;
    I4x[i] = (rho[i] == 0.0) ? 0.0
                             : (dx / rhoi3) * l +
                                   (2.0 / rho[i]) * (k[i] * dx + k1sl * dy * l);
    I4y[i] = 0.0;

    // Courant-Snyder gamma
    gammax[i] = (1.0 + ax * ax) / bx;
    gammay[i] = (1.0 + ay * ay) / by;

    // curly H
    hx[i] = bx * dpx * dpx + 2.0 * ax * dx * dpx + gammax[i] * dx * dx;
    hy[i] = by * dpy * dpy + 2.0 * ay * dy * dpy + gammay[i] * dy * dy;

    I5x[i] = (rho[i] == 0) ? 0.0 : hx[i] * l / rhoi3;
    I5y[i] = (rho[i] == 0) ? 0.0 : hy[i] * l / rhoi3;
  }

  // extend the map table
  table["rho"] = rho;
  table["k"] = k;

  table["gammax"] = gammax;
  table["gammay"] = gammay;
  table["hx"] = hx;
  table["hy"] = hy;

  table["I2"] = I2;
  table["I3"] = I3;
  table["I4x"] = I4x;
  table["I4y"] = I4y;
  table["I5x"] = I5x;
  table["I5y"] = I5y;
}

/*
================================================================================
================================================================================
METHOD TO PRINT A COLUMN FROM THE UPDATED TWISS TABLE (DEBUGGING)
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 06/08/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - string : key
      key in the map - column name  of the Twiss table
    - map<string, <double>> : table
      Twiss or extended Twiss table

  Returns:
  --------
    void

================================================================================
================================================================================
*/
void printTwissMap(string key, map<string, vector<double>> &table) {
  vector<double> toprint = table[key];

  for (auto elem : toprint) {
    cout << elem << endl;
  }
}