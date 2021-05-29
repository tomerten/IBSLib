#include "CoulombLogFunctions.hpp"
#include "NumericFunctions.hpp"
#include <math.h>
#include <stdio.h>

/*
------------------------------------------------------------------------------
ORIGINAL AUTHORS : MADX AUTHORS COPYRIGHT CERN
------------------------------------------------------------------------------
VERSION 1.0 : SIMPSON INTEGRATOR WITH DECADE SPLITTING TO BE ABLE TO ADD
              CONTE-MARTINI AND BJORKEN-MTINGWA - NORMAL SIMPSON GIVES TOO
              LARGE DEVIATIONS COMPARED TO SIMPSON DECADE.

AUTHOR    : TOM MERTENS
DATE      : 16/02/2021
OPYRIGHT  : CERN / HZB

    DESCRIPTION :
        SIMPSON INTEGRATOR SPLITTING THE INTEGRAL IN
        DECADES FOR FASTER CONVERGENCE

    DETAILS:
        MAX DECADES = 30 (10**30)
        INTERVAL SPLITS (NS) = 50

    REF:
        MADX ORIGINAL SOURCE CODE IN TWSINT FUNCTION
        CERN NOTE CERN-AB-2006-002 EQ 8

-------------------------------------------------------------------------------
a (double)      : lambda**2 coefficient integral denominator
b (double)      : lambda coefficient integral denominator
c (double)      : constant integral denominator
cl (double)     : longitudinal growth time factor
cx (double)     : horizontal growth time factor
cy (double)     : vertical growht time factor
cprime (double) : scaling factor
cyy (double)    : scaling factor adapted to sqrt denominator
tl1 (double)    : lambda coefficient integral numerator
tl2 (double)    : constant term integral numerator
tx1 (double)    : lambda coefficient integral numerator
tx2 (double)    : constant term integral numerator
ty1 (double)    : lambda coefficient integral numerator
ty2 (double)    : constant term integral numerator
-------------------------------------------------------------------------------
test (double)   : iteration convergence test value
tau (double[3]) : {long growth factor, hor growth factor, ver growth factor}
if no convergence all are set to 0.0
-------------------------------------------------------------------------------
*/
void SimpsonDecade(double a, double b, double c, double cl, double cx,
                   double cy, double cprime, double cyy, double tl1, double tl2,
                   double tx1, double tx2, double ty1, double ty2,
                   double *tau) {
  const int maxdec = 30, ns = 50;

  const double ten = 10.0;
  const double three = 3.0;
  const double test = 1.0e-7;
  const double coeff[2] = {2.0, 4.0};

  double ar[31], br[30];

  double h, aloop, alam, term, cof, f;
  double suml, sumx, sumy, tmpl, tmpx, tmpy;
  double polyl, polyx, polyy, func;
  bool flag = 0;

  double zintl = 0.0;
  double zintx = 0.0;
  double zinty = 0.0;

  ar[0] = 0.0;

  for (int iloop = 0; iloop < maxdec; iloop++) {
    br[iloop] = pow(ten, iloop);
    ar[iloop + 1] = br[iloop];
    h = (br[iloop] - ar[iloop]) / ns;
    aloop = ar[iloop];

    term = sqrt(cyy * cyy *
                (aloop * aloop * aloop + a * aloop * aloop + b * aloop + c));

    func = sqrt(aloop) / (term * term * term);
    polyl = tl1 * aloop + tl2;
    polyx = tx1 * aloop + tx2;
    polyy = ty1 * aloop + ty2;
    suml = func * polyl;
    sumx = func * polyx;
    sumy = func * polyy;

    for (int iiz = 0; iiz < ns; iiz++) {
      alam = aloop + iiz * h;
      cof = coeff[iiz % 2];
      term = sqrt(cyy * cyy *
                  (alam * alam * alam + a * alam * alam + b * alam + c));

      f = sqrt(alam) / (term * term * term);
      polyl = tl1 * alam + tl2;
      polyx = tx1 * alam + tx2;
      polyy = ty1 * alam + ty2;

      suml = suml + cof * f * polyl;
      sumx = sumx + cof * f * polyx;
      sumy = sumy + cof * f * polyy;
    }

    suml = suml - f * polyl;
    sumx = sumx - f * polyx;
    sumy = sumy - f * polyy;
    tmpl = (suml / three) * h;
    tmpx = (sumx / three) * h;
    tmpy = (sumy / three) * h;
    zintl = zintl + tmpl;
    zintx = zintx + tmpx;
    zinty = zinty + tmpy;

    if ((fabs(tmpl / zintl) < test) & (fabs(tmpx / zintx) < test) &
        (fabs(tmpy / zinty) < test)) {
      // ---- Divide answers by cprime to account for scaling.
      tau[1] = cx * (zintx / cprime);
      tau[0] = cl * (zintl / cprime);
      tau[2] = cy * (zinty / cprime);

      flag = 1;
      break;
    }
  }

  if (flag == 0) {
    tau[0] = 0.0;
    tau[1] = 0.0;
    tau[2] = 0.0;
  }
}
