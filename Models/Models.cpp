#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "NumericFunctions.hpp"
#include <iostream>
#include <math.h>
#include <stdio.h>

/* reminder on how to pass arrays
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
void printouts(double output[3]) {
  // factor 2 and 0.5 are for converting to emittance values from sigma values
  printf("\n");
  printf("    (Weighted) average rates - emittances (1/sec)\n");
  printf("    Longitudinal = %15.6f\n", 2.0 * output[0]);
  printf("    Horizontal   = %15.6f\n", 2.0 * output[1]);
  printf("    Vertical     = %15.6f\n", 2.0 * output[2]);

  printf("\n");
  printf("    (Weighted) average lifetimes - emittances (sec)\n");
  printf("    Longitudinal = %15.6f\n", 0.5 / output[0]);
  printf("    Horizontal   = %15.6f\n", 0.5 / output[1]);
  printf("    Vertical     = %15.6f\n", 0.5 / output[2]);
}

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
output (double[3])      : growth rates for emittances
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
  // factor two is to convert to emit growth rates
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
output (double[3])      : growth rates fro emittances
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
    double *L = &twissdata[i][0];
    double *bx = &twissdata[i][1];
    double *by = &twissdata[i][2];
    double *dx = &twissdata[i][3];

    double rmsx = sqrt(*bx * ex);
    double rmsy = sqrt(*by * ey);
    double d = (rmsx <= rmsy) ? rmsx : rmsy;

    double sigh2inv = (1.0 / (dponp * dponp)) + (*dx * *dx / (rmsx * rmsx));
    double sigh = 1.0 / sqrt(sigh2inv);

    double a = sigh * *bx / (gamma * rmsx);
    double b = sigh * *by / (gamma * rmsy);
    double q = sigh * betar * sqrt(2.0 * d / r0);

    // calc fmohl values
    double fmohlp = fmohl(a, b, q, npp);
    double fmohlx = fmohl(1 / a, b / a, q / a, npp);
    double fmohly = fmohl(1 / b, a / b, q / b, npp);

    // calc IBS growth times ( AMPLITUDE - NOT EMITTANCE )
    alfap0 += ca * fmohlp * (sigh * sigh / (dponp * dponp)) * *L;
    alfax0 +=
        ca * (fmohlx + fmohlp * *dx * *dx * sigh * sigh / (rmsx * rmsx)) * *L;
    alfay0 += ca * fmohly * *L;
  }
  // factor two is to convert to emittance growth rates
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
output (double[3])      : growth rates for emittances
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
double *PiwinskiLatticeModified(double pnumber, double ex, double ey,
                                double sigs, double dponp,
                                double twissheader[5], int n,
                                double (*twissdata)[6], // shape [6,n]
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
    // local naming
    double *L = &twissdata[i][0];
    double *bx = &twissdata[i][1];
    double *by = &twissdata[i][2];
    double *dx = &twissdata[i][3];
    double *dpx = &twissdata[i][4];
    double *alfx = &twissdata[i][5];

    // fmohl accuracy
    int npp = 1000;

    double H0 = *dx;
    double H1 = *bx * *dpx + *alfx * *dx;
    double H = (H0 * H0 + H1 * H1) / *bx;

    double rmsx = sqrt(*bx * ex);
    double rmsy = sqrt(*by * ey);
    double d = (rmsx <= rmsy) ? rmsx : rmsy;

    double sigh2inv = (1.0 / (dponp * dponp)) + (H / ex);
    double sigh = 1.0 / sqrt(sigh2inv);

    double a = sigh * *bx / (gamma * rmsx);
    double b = sigh * *by / (gamma * rmsy);
    double q = sigh * betar * sqrt(2.0 * d / r0);

    // calc fmohl values
    double fmohlp = fmohl(a, b, q, npp);
    double fmohlx = fmohl(1 / a, b / a, q / a, npp);
    double fmohly = fmohl(1 / b, a / b, q / b, npp);

    // calc IBS growth times ( AMPLITUDE - NOT EMITTANCE )
    alfap0 += ca * fmohlp * (sigh * sigh / (dponp * dponp)) * *L;
    alfax0 +=
        ca * (fmohlx + fmohlp * *dx * *dx * sigh * sigh / (rmsx * rmsx)) * *L;
    alfay0 += ca * fmohly * *L;
  }

  output[0] = alfap0 / len;
  output[1] = alfax0 / len;
  output[2] = alfay0 / len;

  return output;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : RODERIK BRUCE, TOM MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE
TYPES (P-PB) AUTHOR    : TOM MERTENS DATE      : 06/02/2021 COPYRIGHT : CERN
/ HZB

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
twiss (double[])        : twiss table data for the required columns (see
DETAILS above) aatom (double)          : atomic number - for electrons this
is electron_mass_energy_MeV / proton_mass_energy_MeV
---------------------------------------------------------------------------------------------------------------
output (double[3])      : growth rates for emittances
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
*/
double *Nagaitsev(double pnumber, double ex, double ey, double sigs,
                  double dponp, double twissheader[5], int n,
                  double (*twissdata)[6], // shape [6,n]
                  double r0) {
  const double c = clight;

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

#pragma omp parallel for shared(twissdata,len) reduction(+: alfap0, alfax0,alfay0)
  for (int i = 0; i < n; i++) {
    // local naming
    double *L = &twissdata[i][0];
    double *bx = &twissdata[i][1];
    double *by = &twissdata[i][2];
    double *dx = &twissdata[i][3];
    double *dpx = &twissdata[i][4];
    double *alfx = &twissdata[i][5];

    double phi = *dpx + (*alfx * (*dx / *bx));
    double axx = *bx / ex;
    double ayy = *by / ey;

    double sigmax = sqrt(*dx * *dx * dponp * dponp + ex * *bx);
    double sigmay = sqrt(ey * *by);

    double as =
        axx * (*dx * *dx / (*bx * *bx) + phi * phi) + (1.0 / (dponp * dponp));
    double a1 = 0.5 * (axx + gamma * gamma * as);
    double a2 = 0.5 * (axx - gamma * gamma * as);
    double b1 = sqrt(a2 * a2 + gamma * gamma * axx * axx * phi * phi);

    double lambda1 = ayy;
    double lambda2 = a1 + b1;
    double lambda3 = a1 - b1;

    double R1 = (1.0 / lambda1) *
                rds((1.0 / lambda2), (1.0 / lambda3), (1.0 / lambda1));
    double R2 = (1.0 / lambda2) *
                rds((1.0 / lambda3), (1.0 / lambda1), (1.0 / lambda2));
    double R3 = 3.0 * sqrt((lambda1 * lambda2) / lambda3) -
                (lambda1 / lambda3) * R1 - (lambda2 / lambda3) * R2;

    double sp = (gamma * gamma / 2.0) * (2.0 * R1 - R2 * (1.0 - 3.0 * a2 / b1) -
                                         R3 * (1.0 + 3.0 * a2 / b1));
    double sx = 0.5 * (2.0 * R1 - R2 * (1.0 + 3.0 * a2 / b1) -
                       R3 * (1.0 - 3.0 * a2 / b1));
    double sxp = (3.0 * gamma * gamma * phi * phi * axx) / b1 * (R3 - R2);

    double alfapp = sp / (sigmax * sigmay);
    double alfaxx = (*bx / (sigmax * sigmay)) *
                    (sx + sxp + sp * (*dx * *dx / (*bx * *bx) + phi * phi));
    double alfayy = (*by / (sigmax * sigmay)) * (-2.0 * R1 + R2 + R3);

    double clog[2];
    twclog(pnumber, *bx, *by, *dx, 0.0, ex, ey, r0, gamma, charge, en0, amass,
           dponp, sigs, clog);
    alfap0 += (alfapp * *L * clog[0]);
    alfax0 += (alfaxx * *L * clog[0]);
    alfay0 += (alfayy * *L * clog[0]);
    // printf("clog %12.6e\n", clog[0]);
  }

  // factor 2.0 is due to converstion to amplitudes from emittances
  output[0] = alfap0 / (dponp * dponp) * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;
  output[1] = alfax0 / ex * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;
  output[2] = alfay0 / ey * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;

  return output;
}

// doubling of code to avoid large amount of if evaluations when
// runnig long simulations
double *Nagaitsevtailcut(double pnumber, double ex, double ey, double sigs,
                         double dponp, double twissheader[5], int n,
                         double (*twissdata)[12], // shape [6,n]
                         double r0, double aatom) {
  const double c = clight;

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

#pragma omp parallel for shared(twissdata,len) reduction(+: alfap0, alfax0,alfay0)
  for (int i = 0; i < n; i++) {
    // local naming
    double *L = &(twissdata[i][0]);
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

    double phi = *dpx + (*ax * (*dx / *bx));
    double axx = *bx / ex;
    double ayy = *by / ey;

    double sigmax = sqrt(*dx * *dx * dponp * dponp + ex * *bx);
    double sigmay = sqrt(ey * *by);

    double as =
        axx * (*dx * *dx / (*bx * *bx) + phi * phi) + (1.0 / (dponp * dponp));
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

    double sp = (gamma * gamma / 2.0) * (2.0 * R1 - R2 * (1.0 - 3.0 * a2 / b1) -
                                         R3 * (1.0 + 3.0 * a2 / b1));
    double sx = 0.5 * (2.0 * R1 - R2 * (1.0 + 3.0 * a2 / b1) -
                       R3 * (1.0 - 3.0 * a2 / b1));
    double sxp = (3.0 * gamma * gamma * phi * phi * axx) / b1 * (R3 - R2);

    double alfapp = sp / (sigmax * sigmay);
    double alfaxx = (*bx / (sigmax * sigmay)) *
                    (sx + sxp + sp * (*dx * *dx / (*bx * *bx) + phi * phi));
    double alfayy = (*by / (sigmax * sigmay)) * (-2.0 * R1 + R2 + R3);

    double clog[2];
    twclogtail(pnumber, *L, *bx, *by, *dx, *dpx, *dy, *dpy, *ax, *ay, *angle,
               *k1l, *k1sl, ex, ey, r0, aatom, gamma, en0, len, amass, charge,
               dponp, sigs, clog);
    alfap0 += (alfapp * *L * clog[0]);
    alfax0 += (alfaxx * *L * clog[0]);
    alfay0 += (alfayy * *L * clog[0]);
    // printf("clog %12.6e\n", clog[0]);
    // twclog(pnumber, *bx, *by, *dx, 0.0, ex, ey, r0, gamma, charge, en0,
    // amass,
    //     dponp, sigs, clog);
    // printf("clog %12.6e\n", clog[0]);
  }

  output[0] = alfap0 / (dponp * dponp) * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;
  output[1] = alfax0 / ex * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;
  output[2] = alfay0 / ey * (pnumber * r0 * r0 * c) /
              (12.0 * pi * betar3 * gamma5 * sigs) / 2.0 / len;

  return output;
}

/*
================================================================================

IBS MODELS USING DEDICATED INTEGRATORS

NOTE:
-----
    USES MADX TWISS TABLE WITH TWISS AT CENTER OF ELEMENTS!!!!!!

================================================================================
*/

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : CERN MADX COPYRIGHT
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE
TYPES (P-PB)

    AUTHOR    :
        TOM MERTENS DATE      : 01/06/2021 COPYRIGHT : CERN/ HZB

    REF:
        CERN NOTE AB-2006--002

    DESCRIPTION :
      IBS ROUTINE MADX

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
            3 -> alfax
            4 -> alfay
            5 -> dx
            6 -> dpx
            7 -> dy
            8 -> dpy

---------------------------------------------------------------------------------------------------------------
pnumber (double)        : number of real particles in the bunch
ex (double)             : horizontal emittance
ey (double)             : vertical emittance
sigs (double)           : bunch length
dponp (double)          : momentum spread
twissheader (double[4]) : twiss header containing (see details above)
n (int)                 : number of rows in the twiss data table
twiss (double[])        : twiss table data for the required columns (see
DETAILS above) aatom (double)          : atomic number - for electrons this
is electron_mass_energy_MeV / proton_mass_energy_MeV
---------------------------------------------------------------------------------------------------------------
output (double[3])      : growth rates for emittances
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
*/
double *ibsmadx(double pnumber, double ex, double ey, double sigs, double sige,
                double twissheader[5], int n, double (*twissdata)[9], double r0,
                bool printout) {
  const double zero = 0.0;
  const double one = 1.0;
  const double two = 2.0;

  static double output[3];

  double gamma = twissheader[0];
  double charge = twissheader[1];
  double circ = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // relativistic beta
  double betar = sqrt(1 - 1 / (gamma * gamma));
  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;
  // NOTE:
  // ****************************************************************
  // Sige is the dE/E. dp/p needed as input for the IBS calculations
  //   dp/p= (dE/E)/beta**2
  // *****************************************************************

  // sige /= (betar * betar);

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;
  double alfas[3];

  double sbxb = zero;
  double sbxinv = zero;
  double sbyb = zero;
  double sbyinv = zero;
  double salxb = zero;
  double salyb = zero;
  double sdxb = zero;
  double sdpxb = zero;
  double sdyb = zero;
  double sdpyb = zero;
  double wnorm = zero;
  double dxwtd = zero;
  double dpxwtd = zero;
  double dywtd = zero;
  double dpywtd = zero;
  double bywtd = zero;
  double alxwtd = zero;
  double alywtd = zero;
  // double hscrpt  = zero;
  // double hscrpty = zero;
  double hscwtd = zero;
  double hscwtdy = zero;

  // CoulombLog(pnumber, ex, ey, twissheader, sige, sigs, aatom, 0, clog);

#pragma omp parallel for shared(twissdata,circ, alfas) reduction(+: alfap0, alfax0, alfay0, sbxb, sbyb, salxb, salyb,sdxb, sdyb )
  for (int i = 0; i < n; i++) {
    double dels = twissdata[i][0]; // l
    double betax = twissdata[i][1];
    double betay = twissdata[i][2];
    double alx = twissdata[i][3]; // alfx
    double aly = twissdata[i][4]; // alfy
    double dx = betar * twissdata[i][5];
    double dpx = betar * twissdata[i][6];
    double dy = betar * twissdata[i][7];
    double dpy = betar * twissdata[i][8];

    sbxb = sbxb + betax * dels;
    sbxinv = sbxinv + dels / betax;
    sbyb = sbyb + betay * dels;
    sbyinv = sbyinv + dels / betay;
    salxb = salxb + alx * dels;
    salyb = salyb + aly * dels;
    sdxb = sdxb + dx * dels;
    sdpxb = sdpxb + dpx * dels;
    sdyb = sdyb + dy * dels;
    sdpyb = sdpyb + dpy * dels;

    //---- Calculate weighted average in region of non-zero DX's.
    //     These values are used to calculate "average" ring lifetimes
    //     in TWSINT.
    if (dx > zero) {
      double wnorm = wnorm + dels;
      double dxwtd = dxwtd + dels * dx;
      double dpxwtd = dpxwtd + dels * dpx;
      double dywtd = dywtd + dels * dy;
      double dpywtd = dpywtd + dels * dpy;
      double bywtd = bywtd + dels / sqrt(betay);
      double alxwtd = alxwtd + dels * alx;
      double alywtd = alywtd + dels * aly;
      double hscrpt = betax * dpx * dpx + two * alx * dx * dpx +
                      (one + alx * alx) * dx * dx / betax;
      double hscrpty = betay * dpy * dpy + two * aly * dy * dpy +
                       (one + aly * aly) * dy * dy / betay;
      double hscwtd = hscwtd + dels * sqrt(hscrpt);
      double hscwtdy = hscwtdy + dels * sqrt(hscrpty);
    }

    //---- IBSIntegrator calculates the Bjorken/Mtingwa integral.
    twsint(pnumber, ex, ey, sigs, sige, gamma, betax, betay, alx, aly, dx, dpx,
           dy, dpy, alfas);
    double clog[2];
    twclog(pnumber, betax, betay, dx, dy, ex, ey, r0, gamma, charge, en0, amass,
           sige, sigs, clog);
    alfap0 += alfas[0] * dels * clog[1];
    alfax0 += alfas[1] * dels * clog[1];
    alfay0 += alfas[2] * dels * clog[1];
  }

  // ---- We have finished reading the lattice
  double bxbar = sbxb / circ;
  double bybar = sbyb / circ;
  double alxbar = salxb / circ;
  double alybar = salyb / circ;
  double dxbar = sdxb / circ;
  double dpxbr = sdpxb / circ;
  double dybar = sdyb / circ;
  double dpybr = sdpyb / circ;
  // double bxinv  = sbxinv / circ;
  // double byinv  = sbyinv / circ;

  dxwtd = dxwtd / wnorm;
  dpxwtd = dpxwtd / wnorm;
  dywtd = dywtd / wnorm;
  dpywtd = dpywtd / wnorm;
  bywtd = bywtd / wnorm;
  bywtd = one / (bywtd * bywtd);
  alxwtd = alxwtd / wnorm;
  alywtd = alywtd / wnorm;
  hscwtd = (hscwtd / wnorm) * (hscwtd / wnorm);

  double beteff = dxwtd * dxwtd / hscwtd;
  double beteffy = (hscwtdy != 0.0) ? dywtd * dywtd / hscwtdy : bywtd;
  double tbar[3], twtd[3];

  tbar[0] = zero;
  tbar[1] = zero;
  tbar[2] = zero;

  twtd[0] = zero;
  twtd[1] = zero;
  twtd[2] = zero;

  output[0] = alfap0 / circ / 2.0;
  output[1] = alfax0 / circ / 2.0;
  output[2] = alfay0 / circ / 2.0;

  if (printout) {
    //---- Integral for averaged quantities.
    twsint(pnumber, ex, ey, sigs, sige, gamma, bxbar, bybar, alxbar, alybar,
           dxbar, dpxbr, dybar, dpybr, tbar);

    // ---- Integral for effective quantities.
    twsint(pnumber, ex, ey, sigs, sige, gamma, beteff, beteffy, alxwtd, alywtd,
           dxwtd, dpxwtd, dywtd, dpywtd, twtd);

    // ---- Calculate the Coulomb logarithm.

    double clog2[2];
    twclog(pnumber, bxbar, bybar, dxbar, dybar, ex, ey, r0, gamma, charge, en0,
           amass, sige, sigs, clog2);

    printf("    r0                  ; %.8e m*rad\n", r0);
    printf("    X-emittance         ; %.8e m*rad\n", ex);
    printf("    Y-emittance         ; %.8e m*rad\n", ey);
    printf("    SIGS                ; %.8e m*rad\n", sigs);
    printf("    SIGE                ; %.8e m*rad\n", sige);
    printf("    ENERGY              : %.8e\n", en0);
    printf("    GAMMA               : %.8e\n", gamma);
    printf("    CONST               : %.8e\n", clog2[1]);
    printf("    COULOMB LOG         : %.8e\n", clog2[0]);
    printf("\n");
    printf("Ring average values \n");
    printf("    betx   = %.8e   bety   = %.8e   Dx  = %.8e   Dy  = %.8e \n",
           bxbar, bybar, dxbar, dybar);
    // CoulombLog(pnumber, ex, ey, avgtwissheader, sige, sigs, aatom, 1,
    // clog);

    // factors 2 and half are due to difference in
    // growth rates / times for emit or sig
    printf("\n");
    printf("    (Weighted) average rates (1/sec)\n");
    printf("    Longitudinal = %15.6f\n", 2.0 * output[0]);
    printf("    Horizontal   = %15.6f\n", 2.0 * output[1]);
    printf("    Vertical     = %15.6f\n", 2.0 * output[2]);

    printf("\n");
    printf("    (Weighted) average lifetimes (sec)\n");
    printf("    Longitudinal = %15.6f\n", 0.5 / output[0]);
    printf("    Horizontal   = %15.6f\n", 0.5 / output[1]);
    printf("    Vertical     = %15.6f\n", 0.5 / output[2]);
  }

  return output;
}

double *ibsmadxtailcut(double pnumber, double ex, double ey, double sigs,
                       double sige, double twissheader[5], int n,
                       double (*twissdata)[12], double r0, double aatom) {
  const double zero = 0.0;
  const double one = 1.0;
  const double two = 2.0;

  static double output[3];

  double gamma = twissheader[0];
  double charge = twissheader[1];
  double circ = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // relativistic beta
  double betar = sqrt(1 - 1 / (gamma * gamma));
  // NOTE:
  // ****************************************************************
  // Sige is the dE/E. dp/p needed as input for the IBS calculations
  //   dp/p= (dE/E)/beta**2
  // *****************************************************************

  // sige /= (betar * betar);

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;
  double alfas[3];

  double sbxb = zero;
  double sbxinv = zero;
  double sbyb = zero;
  double sbyinv = zero;
  double salxb = zero;
  double salyb = zero;
  double sdxb = zero;
  double sdpxb = zero;
  double sdyb = zero;
  double sdpyb = zero;
  double wnorm = zero;
  double dxwtd = zero;
  double dpxwtd = zero;
  double dywtd = zero;
  double dpywtd = zero;
  double bywtd = zero;
  double alxwtd = zero;
  double alywtd = zero;
  // double hscrpt  = zero;
  // double hscrpty = zero;
  double hscwtd = zero;
  double hscwtdy = zero;

  // CoulombLog(pnumber, ex, ey, twissheader, sige, sigs, aatom, 0, clog);

#pragma omp parallel for shared(twissdata,circ, alfas) reduction(+: alfap0, alfax0, alfay0, sbxb, sbyb, salxb, salyb,sdxb, sdyb )
  for (int i = 0; i < n; i++) {
    // local naming
    double *L = &(twissdata[i][0]);
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

    sbxb = sbxb + *bx * *L;
    sbxinv = sbxinv + *L / *bx;
    sbyb = sbyb + *by * *L;
    sbyinv = sbyinv + *L / *by;
    salxb = salxb + *ax * *L;
    salyb = salyb + *ay * *L;
    sdxb = sdxb + *dx * *L;
    sdpxb = sdpxb + *dpx * *L;
    sdyb = sdyb + *dy * *L;
    sdpyb = sdpyb + *dpy * *L;

    //---- Calculate weighted average in region of non-zero DX's.
    //     These values are used to calculate "average" ring lifetimes
    //     in TWSINT.
    if (*dx > zero) {
      double wnorm = wnorm + *L;
      double dxwtd = dxwtd + *L * *dx;
      double dpxwtd = dpxwtd + *L * *dpx;
      double dywtd = dywtd + *L * *dy;
      double dpywtd = dpywtd + *L * *dpy;
      double bywtd = bywtd + *L / sqrt(*by);
      double alxwtd = alxwtd + *L * *ax;
      double alywtd = alywtd + *L * *ay;
      double hscrpt = *bx * *dpx * *dpx + two * *ax * *dx * *dpx +
                      (one + *ax * *ax) * *dx * *dx / *bx;
      double hscrpty = *by * *dpy * *dpy + two * *ay * *dy * *dpy +
                       (one + *ay * *ay) * *dy * *dy / *by;
      double hscwtd = hscwtd + *L * sqrt(hscrpt);
      double hscwtdy = hscwtdy + *L * sqrt(hscrpty);
    }

    //---- IBSIntegrator calculates the Bjorken/Mtingwa integral.
    twsint(pnumber, ex, ey, sigs, sige, gamma, *bx, *by, *ax, *ay, *dx, *dpx,
           *dy, *dpy, alfas);
    double clog[2];
    twclogtail(pnumber, *L, *bx, *by, *dx, *dpx, *dy, *dpy, *ax, *ay, *angle,
               *k1l, *k1sl, ex, ey, r0, aatom, gamma, en0, circ, amass, charge,
               sige, sigs, clog);
    alfap0 += alfas[0] * *L * clog[1];
    alfax0 += alfas[1] * *L * clog[1];
    alfay0 += alfas[2] * *L * clog[1];
  }

  // factor 2 for converting to amplitudes from emit growth rates
  output[0] = alfap0 / circ / 2.0;
  output[1] = alfax0 / circ / 2.0;
  output[2] = alfay0 / circ / 2.0;

  return output;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : TOM MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE
TYPES (P-PB)

    AUTHOR    :
        TOM MERTENS DATE      : 01/06/2021 COPYRIGHT : HZB

    REF:
        CERN NOTE AB-2006--002

    DESCRIPTION :
      IBS ROUTINE MADX

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
            5 -> alfax
---------------------------------------------------------------------------------------------------------------

*/
double *BjorkenMtingwa2(double pnumber, double ex, double ey, double sigs,
                        double dponp, double twissheader[5], int n,
                        double (*twissdata)[6], double r0) {
  const double gamma = twissheader[0];
  const double charge = twissheader[1];
  const double circ = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // initialize
  double alfax0 = 0.0f;
  double alfay0 = 0.0f;
  double alfap0 = 0.0f;
  double gamma2 = gamma * gamma;
  double dponp2 = dponp * dponp;

  static double output[3];

#pragma omp parallel for shared(twissdata) reduction(+ : alfap0, alfax0, alfay0)
  for (int i = 0; i < n; i++) {
    double l = twissdata[i][0];
    double bx = twissdata[i][1];
    double by = twissdata[i][2];
    double dx = twissdata[i][3];
    double dpx = twissdata[i][4];
    double ax = twissdata[i][5];

    double dx2 = dx * dx;
    // double dy2 = dy * dy;
    double bx2 = bx * bx;
    // double by2 = by * by;
    double bxy = bx * by;
    double exy = ex * ey;

    double c1 = (bx / ex + by / ey);
    double gd2 = gamma2 / dponp2;

    double phix = dpx + (ax * (dx / bx));
    // double phiy = dpy + (ay * (dy / by));

    double phix2 = phix * phix;
    // double phiy2 = fabs(phiy * phiy);

    double hx = (dx2 + bx2 * phix2) / bx;
    // double hy = (dy2 + by2 * phiy2) / by;

    double a = (gamma2 * hx / ex) + gd2;
    double b = c1 * (gamma2 * dx2 / (ex * bx) + gd2) +
               (gamma2 * phix2 * bx * by / (ex * ey));
    double c = (bxy / exy) * (gamma2 * dx2 / (ex * bx) + gd2);

    double axx = 2.0 * gamma2 * hx / ex + 2.0 * gd2;
    double bxx =
        c1 * ((gamma2 * hx * hx / ex) + gd2) - bx2 / (ex * ex) * gamma2 * phix2;

    double al = 2.0 * gamma2 * (hx / ex + 1.0 / dponp2);
    double bl = b;

    double ayy = -gamma2 * hx / ex - gd2;
    double byy = b - 3.0 * bx / ex * (gamma2 * dx2 / (ex * bx) + gd2);

    double integrals[3];

    intSimpson(IBSIntegralIntegrand, axx, bxx, ayy, byy, al, bl, a, b, c,
               integrals);

    double clog[2];
    twclog(pnumber, bx, by, dx, 0.0, ex, ey, r0, gamma, charge, en0, amass,
           dponp, sigs, clog);

    alfap0 += l * gd2 * integrals[0] * clog[1];
    alfax0 += l * hx * integrals[1] * clog[1];
    alfay0 += l * by * integrals[2] * clog[1];
  }
  alfap0 /= circ;
  alfax0 *= gamma2 / ex / circ;
  alfay0 *= 1.0 / ey / circ;

  output[0] = alfap0 / 2.0;
  output[1] = alfax0 / 2.0;
  output[2] = alfay0 / 2.0;

  return output;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : TOM MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE
TYPES (P-PB)

    AUTHOR    :
        TOM MERTENS DATE      : 01/06/2021 COPYRIGHT : HZB

    REF:
        CERN NOTE AB-2006--002

    DESCRIPTION :
      IBS ROUTINE MADX

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
            3 -> alfx
            4 -> alfy
            5 -> dx
            6 -> dpx
            7 -> dy
            8 -> dpy
---------------------------------------------------------------------------------------------------------------
*/
double *BjorkenMtingwa(double pnumber, double ex, double ey, double sigs,
                       double dponp, double twissheader[], int n,
                       double (*twissdata)[9], double r0) {
  // constants
  const double gamma = twissheader[0];
  const double charge = twissheader[1];
  const double circ = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;
  double betar = sqrt(1 - 1 / (gamma * gamma));

  static double output[3];

#pragma omp parallel for shared(twissdata) reduction(+ : alfap0, alfax0, alfay0)
  for (int i = 0; i < n; i++) {
    double l = twissdata[i][0];
    double bx = twissdata[i][1];
    double by = twissdata[i][2];
    double ax = twissdata[i][3];
    double ay = twissdata[i][4];
    double dx = betar * twissdata[i][5];
    double dpx = betar * twissdata[i][6];
    double dy = betar * twissdata[i][7];
    double dpy = betar * twissdata[i][8];

    double integrals[3];

    BjorkenMtingwaInt(pnumber, ex, ey, sigs, betar * betar * dponp, gamma, bx,
                      by, ax, ay, dx, dpx, dy, dpy, integrals);

    double clog[2];
    twclog(pnumber, bx, by, dx, 0.0, ex, ey, r0, gamma, charge, en0, amass,
           dponp, sigs, clog);

    alfap0 += l * integrals[0] * clog[1];
    alfax0 += l * integrals[1] * clog[1];
    alfay0 += l * integrals[2] * clog[1];
  }

  alfap0 /= circ;
  alfax0 /= circ;
  alfay0 /= circ;

  output[0] = alfap0 / 2.0;
  output[1] = alfax0 / 2.0;
  output[2] = alfay0 / 2.0;

  return output;
}
double *BjorkenMtingwatailcut(double pnumber, double ex, double ey, double sigs,
                              double dponp, double twissheader[], int n,
                              double (*twissdata)[12], double r0,
                              double aatom) {
  // constants
  const double gamma = twissheader[0];
  const double charge = twissheader[1];
  const double circ = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;
  double betar = sqrt(1 - 1 / (gamma * gamma));

  static double output[3];

#pragma omp parallel for shared(twissdata) reduction(+ : alfap0, alfax0, alfay0)
  for (int i = 0; i < n; i++) {
    double *L = &(twissdata[i][0]);
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

    double integrals[3];

    BjorkenMtingwaInt(pnumber, ex, ey, sigs, betar * betar * dponp, gamma, *bx,
                      *by, *ax, *ay, *dx, *dpx, *dy, *dpy, integrals);

    double clog[2];
    twclogtail(pnumber, *L, *bx, *by, *dx, *dpx, *dy, *dpy, *ax, *ay, *angle,
               *k1l, *k1sl, ex, ey, r0, aatom, gamma, en0, circ, amass, charge,
               dponp, sigs, clog);

    alfap0 += *L * integrals[0] * clog[1];
    alfax0 += *L * integrals[1] * clog[1];
    alfay0 += *L * integrals[2] * clog[1];
  }

  // factor 2 for converting to amplitudes from emit growth rates
  alfap0 /= circ;
  alfax0 /= circ;
  alfay0 /= circ;

  output[0] = alfap0 / 2.0;
  output[1] = alfax0 / 2.0;
  output[2] = alfay0 / 2.0;

  return output;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : TOM MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE
TYPES (P-PB)

    AUTHOR    :
        TOM MERTENS DATE      : 01/06/2021 COPYRIGHT : HZB

    REF:
        CERN NOTE AB-2006--002

    DESCRIPTION :
      IBS ROUTINE MADX

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
            3 -> alfx
            4 -> alfy
            5 -> dx
            6 -> dpx
            7 -> dy
            8 -> dpy

---------------------------------------------------------------------------------------------------------------
*/
double *ConteMartini(double pnumber, double ex, double ey, double sigs,
                     double dponp, double twissheader[], int n,
                     double (*twissdata)[9], double r0) {
  // constants
  const double gamma = twissheader[0];
  const double charge = twissheader[1];
  const double circ = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;
  double betar = sqrt(1 - 1 / (gamma * gamma));
  static double output[3];

#pragma omp parallel for shared(twissdata) reduction(+ : alfap0, alfax0, alfay0)
  for (int i = 0; i < n; i++) {
    double l = twissdata[i][0]; // l
    double bx = twissdata[i][1];
    double by = twissdata[i][2];
    double ax = twissdata[i][3]; // alfx
    double ay = twissdata[i][4]; // alfy
    double dx = betar * twissdata[i][5];
    double dpx = betar * twissdata[i][6];
    double dy = betar * twissdata[i][7];
    double dpy = betar * twissdata[i][8];

    double integrals[3];

    ConteMartiniInt(pnumber, ex, ey, sigs, betar * betar * dponp, gamma, bx, by,
                    ax, ay, dx, dpx, dy, dpy, integrals);

    double clog[2];
    twclog(pnumber, bx, by, dx, 0.0, ex, ey, r0, gamma, charge, en0, amass,
           dponp, sigs, clog);

    alfap0 += l * integrals[0] * clog[1];
    alfax0 += l * integrals[1] * clog[1];
    alfay0 += l * integrals[2] * clog[1];
  }

  alfap0 /= circ;
  alfax0 /= circ;
  alfay0 /= circ;

  output[0] = alfap0 / 2.0;
  output[1] = alfax0 / 2.0;
  output[2] = alfay0 / 2.0;

  return output;
}

double *ConteMartinitailcut(double pnumber, double ex, double ey, double sigs,
                            double dponp, double twissheader[], int n,
                            double (*twissdata)[12], double r0, double aatom) {
  // constants
  const double gamma = twissheader[0];
  const double charge = twissheader[1];
  const double circ = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;
  double betar = sqrt(1 - 1 / (gamma * gamma));
  static double output[3];

#pragma omp parallel for shared(twissdata) reduction(+ : alfap0, alfax0, alfay0)
  for (int i = 0; i < n; i++) {
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

    double integrals[3];

    ConteMartiniInt(pnumber, ex, ey, sigs, betar * betar * dponp, gamma, *bx,
                    *by, *ax, *ay, *dx, *dpx, *dy, *dpy, integrals);

    double clog[2];
    twclogtail(pnumber, *l, *bx, *by, *dx, *dpx, *dy, *dpy, *ax, *ay, *angle,
               *k1l, *k1sl, ex, ey, r0, aatom, gamma, en0, circ, amass, charge,
               dponp, sigs, clog);

    alfap0 += *l * integrals[0] * clog[1];
    alfax0 += *l * integrals[1] * clog[1];
    alfay0 += *l * integrals[2] * clog[1];
  }

  alfap0 /= circ;
  alfax0 /= circ;
  alfay0 /= circ;

  output[0] = alfap0 / 2.0;
  output[1] = alfax0 / 2.0;
  output[2] = alfay0 / 2.0;

  return output;
}

/*
---------------------------------------------------------------------------------------------------------------
ORIGINAL AUTHORS : TOM MERTENS
---------------------------------------------------------------------------------------------------------------
VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE
TYPES (P-PB)

    AUTHOR    :
        TOM MERTENS DATE      : 01/06/2021 COPYRIGHT : HZB

    REF:
        CERN NOTE AB-2006--002

    DESCRIPTION :
      IBS ROUTINE MADX

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
            3 -> alfx
            4 -> alfy
            5 -> dx
            6 -> dpx
            7 -> dy
            8 -> dpy

---------------------------------------------------------------------------------------------------------------
*/
double *MadxIBS(double pnumber, double ex, double ey, double sigs, double dponp,
                double twissheader[5], int n, double (*twissdata)[9],
                double r0) {
  // constants
  const double gamma = twissheader[0];
  const double charge = twissheader[1];
  const double circ = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];

  // initialize
  double alfax0 = 0.0;
  double alfay0 = 0.0;
  double alfap0 = 0.0;
  double betar = sqrt(1 - 1 / (gamma * gamma));

  static double output[3];

#pragma omp parallel for shared(twissdata) reduction(+ : alfap0, alfax0, alfay0)
  for (int i = 0; i < n; i++) {
    double l = twissdata[i][0];
    double bx = twissdata[i][1];
    double by = twissdata[i][2];
    double ax = twissdata[i][3];
    double ay = twissdata[i][4];
    double dx = twissdata[i][5];
    double dpx = twissdata[i][6];
    double dy = twissdata[i][7];
    double dpy = twissdata[i][8];

    double integrals[3];

    MadxInt(pnumber, ex, ey, sigs, betar * betar * dponp, gamma, bx, by, ax, ay,
            dx, dpx, dy, dpy, integrals);

    double clog[2];
    twclog(pnumber, bx, by, dx, 0.0, ex, ey, r0, gamma, charge, en0, amass,
           dponp, sigs, clog);

    alfap0 += l * integrals[0] * clog[1];
    alfax0 += l * integrals[1] * clog[1];
    alfay0 += l * integrals[2] * clog[1];
  }

  alfap0 /= circ;
  alfax0 /= circ;
  alfay0 /= circ;

  output[0] = alfap0 / 2.0;
  output[1] = alfax0 / 2.0;
  output[2] = alfay0 / 2.0;

  return output;
}