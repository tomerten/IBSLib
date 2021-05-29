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
/*
------------------------------------------------------------------------------
VERSION 1.0 : IBS INTEGRALS INTEGRAND AS FUNCTION TO BE ABLE TO ADD
              CONTE-MARTINI AND BJORKEN-MTINGWA - NORMAL SIMPSON GIVES TOO
              LARGE DEVIATIONS COMPARED TO SIMPSON DECADE.

AUTHOR    : TOM MERTENS
DATE      : 16/02/2021
OPYRIGHT  : HZB

    DESCRIPTION :
        IBS INTEGRALS INTEGRAND

-------------------------------------------------------------------------------
a (double)       : lambda**2 coefficient integral denominator
b (double)       : lambda coefficient integral denominator
c (double)       : constant integral denominator
ax (double)      : lambda coefficient integral numerator
bx (double)      : constant term integral numerator
lambda (double)  : integration variable
-------------------------------------------------------------------------------
*/

double IBSIntegralIntegrand(double lambda, double ax, double bx, double a,
                            double b, double c) {
  double num = sqrt(lambda) * (ax * lambda + bx);
  double term = lambda * lambda * lambda + a * lambda * lambda + b * lambda + c;
  double denom = sqrt(term * term * term);
  return num / denom;
}

double simpson(double ibsintegrand(double, double, double, double, double,
                                   double),
               double ax, double bx, double a, double b, double c, double al,
               double bl, int n) {
  double h, integral, x, sum = 0.0;
  int i;

  h = fabs(bl - al) / n;
  for (i = 0; i < n; i++) {
    x = al + i * h;
    if (i % 2 == 0) {
      sum += 2.0 * ibsintegrand(x, ax, bx, a, b, c);
    } else {
      sum += 4.0 * ibsintegrand(x, ax, bx, a, b, c);
    }
  }

  double fa = ibsintegrand(al, ax, bx, a, b, c);
  double fb = ibsintegrand(bl, ax, bx, a, b, c);

  integral = (h / 3.0) * (fa + fb + sum);

  return integral;
}

void intSimpson(double BjorkenMtingwaIntegrand(double, double, double, double,
                                               double, double),
                double ax, double bx, double ay, double by, double as,
                double bs, double a, double b, double c, double *integral) {
  double al[31], bl[30], aloop, bloop;

  int maxdec = 30, ns = 50;
  double test = 1.0e-7;
  bool flag = 0;
  // double integral[3] = {0.0,0.0,0.0};

  al[0] = 0.0;
  for (int iloop = 0; iloop < maxdec; iloop++) {
    bl[iloop] = pow(10.0, iloop);
    al[iloop + 1] = bl[iloop];
    aloop = al[iloop];
    bloop = bl[iloop];

    double nintx =
        simpson(BjorkenMtingwaIntegrand, ax, bx, a, b, c, aloop, bloop, ns);
    double ninty =
        simpson(BjorkenMtingwaIntegrand, ay, by, a, b, c, aloop, bloop, ns);
    double nints =
        simpson(BjorkenMtingwaIntegrand, as, bs, a, b, c, aloop, bloop, ns);

    integral[0] += nints;
    integral[1] += nintx;
    integral[2] += ninty;

    if ((fabs(nints / integral[0]) < test) &
        (fabs(nintx / integral[1]) < test) & (fabs(ninty / integral[2]) < test)

    ) {
      flag = 1;
      break;
    }
  }
  if (flag == 0) {
    integral[0] = 0.0;
    integral[1] = 0.0;
    integral[2] = 0.0;
  }
  // return integral;
}