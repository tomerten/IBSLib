#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "NumericFunctions.hpp"
#include <iostream>
#include <math.h>
#include <stdio.h>

/* Freminder on how to pass arrays
double *arraytest(int array[][4]) {
  static double output[4];
  for (int i = 0; i < 2; ++i) {
    output[i] = 0;
    for (int j = 0; j < 4; ++j) {
      output[i] += array[i][j];
    }
  }
  return output;
}
*/
/*
-------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
-------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES
(P-PB)
AUTHOR    : TOM MERTENS
DATE      : 05/02/2021
COPYRIGHT : CERN / HZB

    DESCRIPTION :
        IBS ROUTINE PIWINSKI USING SMOOTH PARAMETERS

    REF:
        HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS P.126

    DETAILS:
    TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING MADX TWISS HEADER
        0 -> gamma
        1 -> length
        3 -> gammatr
        4 -> qx
        6 -> qy
-------------------------------------------------------------------------------
pnumber (double)        : number of real particles in the bunch
ex (double)             : horizontal emittance
ey (double)             : vertical emittance
sigs (double)           : bunch length
dponp (double)          : momentum spread
twiss (double[6])       : twiss header (see DETAILS above)
aatom (double)          : atomic number - for electrons this is
electron_mass_energy_MeV / proton_mass_energy_MeV
-------------------------------------------------------------------------------
output (double[3])      : growth rates
    0 -> ap
    1 -> ax
    2 -> ay
-------------------------------------------------------------------------------
USED VARIABLES:
pi       : constant
c        : speed of light
betar    : relativistic beta
betxAvg  : average betax derived from horizontal tune
betyAvg  : average betay derived from vertical tune
xdisp    : average horizontal dispersion derived from relativistic transition
gamma    : relativistic gamma
q1       : horizontal tune
q2       : vertical tune
rmsx     : RMS x
rmsy     : RMX y
r0       : classical radius
-------------------------------------------------------------------------------
*/
double *PiwinskiSmooth(double pnumber, double ex, double ey, double sigs,
                       double dponp, double twiss[5], double r0) {
  const double c = 299792458.0;
  const double pi = 3.141592653589793;

  static double output[3];

  double gamma = twiss[0];
  double len = twiss[1];
  double gammatr = twiss[2];
  double q1 = twiss[3];
  double q2 = twiss[4];

  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;
  double betaxAvg = len / (2.0 * pi * q1);             // avg betax
  double betayAvg = len / (2.0 * pi * q2);             // avgbety
  double betar = sqrt(1.0 - 1.0 / (gamma * gamma));    // relativistic beta
  double xdisp = len / (2.0 * pi * gammatr * gammatr); // avg dispersion

  // RMS transverse beam size
  double rmsx = sqrt(ex * betaxAvg); // rms x
  double rmsy = sqrt(ey * betayAvg); // rms y

  // longitudinal beam size - bunch height
  double sigh2inv = (1.0 / (dponp * dponp)) + (xdisp * xdisp / (rmsx * rmsx));
  double sigh = 1.0 / sqrt(sigh2inv);

  double atop = r0 * r0 * c * pnumber;
  double abot = 64.0 * pi * pi * betar * betar * betar * gamma * gamma * gamma *
                gamma * ex * ey * sigs * dponp;
  double ca = atop / abot;
  double a = sigh * betaxAvg / (gamma * rmsx);
  double b = sigh * betayAvg / (gamma * rmsy);
  double d = (rmsx <= rmsy) ? rmsx : rmsy;
  double q = sigh * betar * sqrt(2.0 * d / r0);

  // fmohl accuracy
  int npp = 1000;

  // calc fmohl values
  double fmohlp = fmohl(a, b, q, npp);
  double fmohlx = fmohl(1 / a, b / a, q / a, npp);
  double fmohly = fmohl(1 / b, a / b, q / b, npp);

  // calc IBS growth times ( AMPLITUDE - NOT EMITTANCE )
  output[0] = ca * fmohlp * (sigh * sigh / (dponp * dponp));
  output[1] =
      ca * (fmohlx + fmohlp * xdisp * xdisp * sigh * sigh / (rmsx * rmsx));
  output[2] = ca * fmohly;

  return output;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES
(P-PB)
AUTHOR    : TOM MERTENS
DATE      : 05/02/2021
COPYRIGHT : CERN / HZB

    DESCRIPTION :
        IBS ROUTINE PIWINSKI FOR INDIVIDUAL ELEMENTS

    REF:
        HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS

    DETAILS:
        TWISS HEADER
        TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
            0 -> gamma
            1 -> length

        TWISS TABLE
        TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
            0 -> length
            1 -> betax
            2 -> betay
            3 -> dx
---------------------------------------------------------------------------------------------------------------
pnumber (double)        : number of real particles in the bunch
ex (double)             : horizontal emittance
ey (double)             : vertical emittance
sigs (double)           : bunch length
dponp (double)          : momentum spread
twissheader (double[4]) : twiss header containing (see details above)
n (int)                 : number of columns in the twiss data table
twiss (double[])        : twiss table data for the required columns (see DETAILS
above) aatom (double)          : atomic number - for electrons this is
electron_mass_energy_MeV / proton_mass_energy_MeV
---------------------------------------------------------------------------------------------------------------
output (double[3])      : growth rates
    0 -> ap
    1 -> ax
    2 -> ay
---------------------------------------------------------------------------------------------------------------
USED VARIABLES:
pi       : constant
c        : speed of light
betar    : relativistic beta
betxAvg  : average betax derived from horizontal tune
betyAvg  : average betay derived from vertical tune
xdisp    : average horizontal dispersion derived from relativistic transition
gamma    : relativistic gamma
q1       : horizontal tune
q2       : vertical tune
rmsx     : RMS x
rmsy     : RMX y
r0       : classical radius
---------------------------------------------------------------------------------------------------------------
*/
double *PiwinskiLattice(double pnumber, double ex, double ey, double sigs,
                        double dponp, double twissheader[2], int n,
                        double (*twissdata)[4], // shape [4,n]
                        double r0) {
  const double c = clight;

  static double output[3];
  double gamma = twissheader[0];
  double len = twissheader[1];

  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;
  double betar = sqrt(1 - 1 / (gamma * gamma));

  double atop = r0 * r0 * c * pnumber;
  double abot = 64.0 * pi * pi * betar * betar * betar * gamma * gamma * gamma *
                gamma * ex * ey * sigs * dponp;
  double ca = atop / abot;

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;

#pragma omp parallel for shared(twissdata) reduction(+ : alfax0, alfay0, alfap0)
  for (int i = 0; i < n; i++) {
    // fmohl accuracy
    int npp = 1000;

    // local naming of twiss data
    // making code more readable (?efficiency)
    double L = twissdata[i][0];
    double bx = twissdata[i][1];
    double by = twissdata[i][2];
    double dx = twissdata[i][3];

    double rmsx = sqrt(bx * ex);
    double rmsy = sqrt(by * ey);
    double d = (rmsx <= rmsy) ? rmsx : rmsy;

    double sigh2inv = (1.0 / (dponp * dponp)) + (dx * dx / (rmsx * rmsx));
    double sigh = 1.0 / sqrt(sigh2inv);

    double a = sigh * bx / (gamma * rmsx);
    double b = sigh * by / (gamma * rmsy);
    double q = sigh * betar * sqrt(2.0 * d / r0);

    // calc fmohl values
    double fmohlp = fmohl(a, b, q, npp);
    double fmohlx = fmohl(1 / a, b / a, q / a, npp);
    double fmohly = fmohl(1 / b, a / b, q / b, npp);

    // calc IBS growth times ( AMPLITUDE - NOT EMITTANCE )
    alfap0 += ca * fmohlp * (sigh * sigh / (dponp * dponp)) * L;
    alfax0 +=
        ca * (fmohlx + fmohlp * dx * dx * sigh * sigh / (rmsx * rmsx)) * L;
    alfay0 += ca * fmohly * L;
  }

  output[0] = alfap0 / len;
  output[1] = alfax0 / len;
  output[2] = alfay0 / len;

  return output;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM
MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES
(P-PB)
AUTHOR    : TOM MERTENS
DATE      : 05/02/2021
COPYRIGHT : CERN / HZB

    DESCRIPTION :
        IBS ROUTINE PIWINSKI FOR INDIVIDUAL ELEMENTS

    REF:
        HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS

    DETAILS:
        TWISS HEADER
        TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
            0 -> gamma
            1 -> length

        TWISS TABLE
        TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
            0 -> length
            1 -> betax
            2 -> betay
            3 -> dx
            4 -> dpx
            5 -> alfx
---------------------------------------------------------------------------------------------------------------
pnumber (double)        : number of real particles in the bunch
ex (double)             : horizontal emittance
ey (double)             : vertical emittance
sigs (double)           : bunch length
dponp (double)          : momentum spread
twissheader (double[4]) : twiss header containing (see details above)
n (int)                 : number of columns in the twiss data table
twiss (double[])        : twiss table data for the required columns (see DETAILS
above) aatom (double)          : atomic number - for electrons this is
electron_mass_energy_MeV / proton_mass_energy_MeV
---------------------------------------------------------------------------------------------------------------
output (double[3])      : growth rates
    0 -> ap
    1 -> ax
    2 -> ay
---------------------------------------------------------------------------------------------------------------
USED VARIABLES:
pi       : constant
c        : speed of light
betar    : relativistic beta
betxAvg  : average betax derived from horizontal tune
betyAvg  : average betay derived from vertical tune
xdisp    : average horizontal dispersion derived from relativistic transition
gamma    : relativistic gamma
q1       : horizontal tune
q2       : vertical tune
rmsx     : RMS x
rmsy     : RMX y
r0       : classical radius
---------------------------------------------------------------------------------------------------------------

double *PiwinskiLatticeModified(double pnumber, double ex, double ey,
                                double sigs, double dponp,
                                double twissheader[5], int n,
                                double (*twissdata)[n], // shape [6,n]
                                double r0) {
  const double c = 299792458.0;
  const double pi = 3.141592653589793;

  static double output[3];
  double gamma = twissheader[0];
  double len = twissheader[1];

  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;
  double betar = sqrt(1 - 1 / (gamma * gamma));

  double atop = r0 * r0 * c * pnumber;
  double abot = 64.0 * pi * pi * betar * betar * betar * gamma * gamma * gamma *
                gamma * ex * ey * sigs * dponp;
  double ca = atop / abot;

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;

#pragma omp parallel for shared(twissdata) reduction(+ : alfax0, alfay0, alfap0)
  for (int i = 0; i < n; i++) {
    // fmohl accuracy
    int npp = 1000;

    double H0 = twissdata[3][i];
    double H1 =
        twissdata[1][i] * twissdata[4][i] + twissdata[5][i] * twissdata[3][i];
    double H = (H0 * H0 + H1 * H1) / twissdata[1][i];

    double rmsx = sqrt(twissdata[1][i] * ex);
    double rmsy = sqrt(twissdata[2][i] * ey);
    double d = (rmsx <= rmsy) ? rmsx : rmsy;

    double sigh2inv = (1.0f / (dponp * dponp)) + (H / ex);
    double sigh = 1.0f / sqrt(sigh2inv);

    double a = sigh * twissdata[1][i] / (gamma * rmsx);
    double b = sigh * twissdata[2][i] / (gamma * rmsy);
    double q = sigh * betar * sqrt(2.0f * d / r0);

    // calc fmohl values
    double fmohlp = fmohl(a, b, q, npp);
    double fmohlx = fmohl(1 / a, b / a, q / a, npp);
    double fmohly = fmohl(1 / b, a / b, q / b, npp);

    // calc IBS growth times ( AMPLITUDE - NOT EMITTANCE )
    alfap0 += ca * fmohlp * (sigh * sigh / (dponp * dponp)) * twissdata[0][i];
    alfax0 += ca *
              (fmohlx + fmohlp * twissdata[3][i] * twissdata[3][i] * sigh *
                            sigh / (rmsx * rmsx)) *
              twissdata[0][i];
    alfay0 += ca * fmohly * twissdata[0][i];
  }

  output[0] = alfap0 / len;
  output[1] = alfax0 / len;
  output[2] = alfay0 / len;

  return output;
}
*/
/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : RODERIK BRUCE, TOM MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES
(P-PB)
AUTHOR    : TOM MERTENS
DATE      : 06/02/2021
COPYRIGHT : CERN / HZB

    DESCRIPTION :
      IBS ROUTINE NAGAITSEV

    DETAILS:
        TWISS HEADER
        TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
            0 -> gamma
            1 -> charge
            2 -> length
            3 -> energy
            4 -> mass

        TWISS TABLE
        TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
            0 -> length
            1 -> betax
            2 -> betay
            3 -> dx
            4 -> dpx
            5 -> alfx

    REF:
        PRSTAB 8, 064403 (2005)

---------------------------------------------------------------------------------------------------------------
pnumber (double)        : number of real particles in the bunch
ex (double)             : horizontal emittance
ey (double)             : vertical emittance
sigs (double)           : bunch length
dponp (double)          : momentum spread
twissheader (double[4]) : twiss header containing (see details above)
n (int)                 : number of columns in the twiss data table
twiss (double[])        : twiss table data for the required columns (see DETAILS
above) aatom (double)          : atomic number - for electrons this is
electron_mass_energy_MeV / proton_mass_energy_MeV
---------------------------------------------------------------------------------------------------------------
output (double[3])      : growth rates
    0 -> ap
    1 -> ax
    2 -> ay

---------------------------------------------------------------------------------------------------------------
pnumber  : number of real particles in the bunch
epsx     : hor emit
epsy     : ver emit
sigs     : bunch len
dponp    : momentum spread
pi       : constant
circ     : accelerator circumference
clight   : speed of light
qatom    : charge
aatom    : atomic number
betar    : relativistic beta
betx     : BETX madx
bety     : BETY madx
dx       : DX madx
dxp      : DPX madx
l        : L madx
alfx     : ALFX madx
nelem    : number of elements
gamma0   : relativistic gamma
alfap0   : longitudinal growth rate
alfax0   : hor growth rate
alfay0   : ver growth rate
---------------------------------------------------------------------------------------------------------------

double *Nagaitsev(double pnumber, double ex, double ey, double sigs,
                  double dponp, double twissheader[5], int n,
                  double (*twissdata)[n], // shape [6,n]
                  double r0) {
  const double c = 299792458.0f;
  const double pi = 3.141592653589793f;

  static double output[3];
  double gamma = twissheader[0];
  double charge = twissheader[1];
  double len = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;

  // initialize
  double alfax0 = 0.0f;
  double alfay0 = 0.0f;
  double alfap0 = 0.0f;
  double betar = sqrt(1 - 1 / (gamma * gamma));
  double betar3 = betar * betar * betar;
  double gamma5 = gamma * gamma * gamma * gamma * gamma;

#pragma omp parallel for shared(twissdata,len) reduction(+: alfap0, alfax0,
alfay0) for (int i = 0; i < n; i++) { double bx = twissdata[1][i]; double by =
twissdata[2][i]; double dx = twissdata[3][i];

    double phi = twissdata[4][i] +
                 (twissdata[5][i] * (twissdata[3][i] / twissdata[1][i]));
    double axx = twissdata[1][i] / ex;
    double ayy = twissdata[2][i] / ey;

    double sigmax = sqrt(twissdata[3][i] * twissdata[3][i] * dponp * dponp +
                         ex * twissdata[1][i]);
    double sigmay = sqrt(ey * twissdata[2][i]);

    double as = axx * (twissdata[3][i] * twissdata[3][i] /
                           (twissdata[1][i] * twissdata[1][i]) +
                       phi * phi) +
                (1.0f / (dponp * dponp));
    double a1 = 0.5 * (axx + gamma * gamma * as);
    double a2 = 0.5 * (axx - gamma * gamma * as);
    double b1 = sqrt(a2 * a2 + gamma * gamma * axx * axx * phi * phi);

    double lambda1 = ayy;
    double lambda2 = a1 + b1;
    double lambda3 = a1 - b1;

    double R1 = (1.0f / lambda1) *
                rds((1.0f / lambda2), (1.0f / lambda3), (1.0f / lambda1));
    double R2 = (1.0f / lambda2) *
                rds((1.0f / lambda3), (1.0f / lambda1), (1.0f / lambda2));
    double R3 = 3.0f * sqrt((lambda1 * lambda2) / lambda3) -
                (lambda1 / lambda3) * R1 - (lambda2 / lambda3) * R2;

    double sp =
        (gamma * gamma / 2.0f) * (2.0f * R1 - R2 * (1.0f - 3.0f * a2 / b1) -
                                  R3 * (1.0f + 3.0f * a2 / b1));
    double sx = 0.5f * (2.0f * R1 - R2 * (1.0f + 3.0f * a2 / b1) -
                        R3 * (1.0f - 3.0f * a2 / b1));
    double sxp = (3.0f * gamma * gamma * phi * phi * axx) / b1 * (R3 - R2);

    double alfapp = sp / (sigmax * sigmay);
    double alfaxx = (twissdata[1][i] / (sigmax * sigmay)) *
                    (sx + sxp +
                     sp * (twissdata[3][i] * twissdata[3][i] /
                               (twissdata[1][i] * twissdata[1][i]) +
                           phi * phi));
    double alfayy =
        (twissdata[2][i] / (sigmax * sigmay)) * (-2.0f * R1 + R2 + R3);

    double clog[2];
    twclog(pnumber, bx, by, dx, 0.0, ex, ey, r0, gamma, charge, en0, amass,
           dponp, sigs, clog);
    alfap0 += (alfapp * twissdata[0][i] * clog[0]);
    alfax0 += (alfaxx * twissdata[0][i] * clog[0]);
    alfay0 += (alfayy * twissdata[0][i] * clog[0]);
  }

  output[0] = alfap0 / (dponp * dponp) * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;
  output[1] = alfax0 / ex * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;
  output[2] = alfay0 / ey * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;

  return output;
}
*/