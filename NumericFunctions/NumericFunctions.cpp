#include "NumericFunctions.hpp"
#include <math.h>
#include <stdio.h>

/*
--------------------------------------------------------------------------------
AUTHOR  : TOM MERTENS
--------------------------------------------------------------------------------
VERSION 1.0 : FOR USING IN CFFI IBS PYTHON TEST EXAMPLE
DATE        : 16/02/2021
COPYRIGHT   : HELMHOLTZ ZENTRUM BERLIN

    DESCRIPTION:
        CALCULATES SIMGA_E FROM SIGMA_S GIVEN OMEGA0, LONGITUNAL TUNE
        AND ETA

--------------------------------------------------------------------------------
qs (double)     : longitunal tune
sigs (double)   : bunch length [m]
omega0 (double) : angular RF frequency
eta (double)    : phase slip factor
--------------------------------------------------------------------------------
RETURNS:
sige (double)   : dE/E (do not confuse with dp/p, there is a factor beta**2
difference)
--------------------------------------------------------------------------------
*/
double sigefromsigs(double omega0, double sigs, double qs, double eta) {
  return qs * omega0 * (sigs / (fabs(eta) * clight));
}

/*
--------------------------------------------------------------------------------
AUTHOR  : TOM MERTENS
--------------------------------------------------------------------------------
VERSION 1.0 : FOR USING IN CFFI IBS PYTHON TEST EXAMPLE
DATE        : 16/02/2021
COPYRIGHT   : HELMHOLTZ ZENTRUM BERLIN

    DESCRIPTION:
        CALCULATES ETA - PHASE SLIP FACTOR

--------------------------------------------------------------------------------
gamma (double)   : relativistic gamma
gammatr (double) : relativistic transition gamma of the lattice
--------------------------------------------------------------------------------
RETURNS:
eta (double)     : phase slip factor (1/ gammatr**2  - 1/ gamma**2)
--------------------------------------------------------------------------------
*/
double eta(double gamma, double gammatr) {
  return 1.0f / (gammatr * gammatr - 1.0f / (gamma * gamma));
}

/*
--------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
--------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES
(P-PB)
AUTHOR    : TOM MERTENS
DATE      : 05/02/2021 COPYRIGHT : CERN / HZB

    DESCRIPTION :
        CALCULATE CLASSICAL PARTICLE RADIUS

    REF:
        HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS

---------------------------------------------------------------------------------
a (double)      : parameter of fmohl
b (double)      : parameter of fmohl
q (double)      : parameter of fmohl
n (int)         : number of iterations (accuracy)
---------------------------------------------------------------------------------
RETURNS:
    fmohl (double)
---------------------------------------------------------------------------------
*/
double fmohl(double a, double b, double q, int n) {
  // constants
  double const euler = 0.577215664901533;
  double const pi = 3.141592653589793;

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