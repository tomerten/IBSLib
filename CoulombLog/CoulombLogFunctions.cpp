#include "NumericFunctions.hpp"
#include "RadiationDamping.hpp"
#include <math.h>
#include <stdio.h>

void twclog(double pnumber, double bx, double by, double dx, double dy,
            double ex, double ey, double r0, double gamma, double charge,
            double en0, double amass, double sige, double sigt, double *clog) {
  // constants
  const double ot2 = 1.0e2;
  const double ft8 = 5.0e8;
  const double ot5 = 1.0e5;
  const double ttm3 = 2.0e-3;
  const double fac1 = 743.4;
  const double fac2 = 1.44e-7;
  const double c = clight;

  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // IMPORTANT: HBAR USED HERE IS REDUCED PLANCK CONSTANT IN GEV!!!!
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // const double hbar = 6.582119569e-25; // 1.0545718176461565e-34;

  const double two = 2.0;
  const double eight = 8.0f;

  double betar = sqrt(1.0 - 1.0 / (gamma * gamma));
  double qion = fabs(charge);
  double etrans = ft8 * (gamma * en0 - amass) * (ex / bx);
  double tempev = two * etrans;

  //---- Calculate beam volume to get density (in cm**-3).
  double sigxcm = ot2 * sqrt(ex * bx + (dx * dx * sige * sige));
  double sigycm = ot2 * sqrt(ey * by + (dy * dy * sige * sige));
  double sigtcm = ot2 * sigt;

  double vol = eight * sqrt(pi * pi * pi) * sigxcm * sigycm * sigtcm;
  double densty = pnumber / vol;

  //---- Calculate RMAX as smaller of SIGXCM and DEBYE length.
  double debyel = fac1 * sqrt(tempev / densty) / qion;
  double rmax = min(sigxcm, debyel);

  //---- Calculate RMIN as larger of classical distance of closest approach
  //     or quantum mechanical diffraction limit from nuclear radius.
  double rmincl = fac2 * qion * qion / tempev;
  double rminqm = hbar * c * ot5 / (two * sqrt(ttm3 * etrans * amass));
  double rmin = max(rmincl, rminqm);
  double coulog = log(rmax / rmin);
  // calculate coulomb log constant pre-factor
  double constt = pnumber * coulog * r0 * r0 * c /
                  (eight * pi * betar * betar * betar * gamma * gamma * gamma *
                   gamma * ex * ey * sige * sigt);

  clog[0] = coulog;
  clog[1] = constt;
}

void twclogtail(double pnumber, double l, double bx, double by, double dx,
                double dpx, double dy, double dpy, double ax, double ay,
                double angle, double k1l, double k1sl, double ex, double ey,
                double r0, double aatom, double gamma, double en0, double len,
                double amass, double charge, double sige, double sigt,
                double *clog) {
  // constants
  const double ot2 = 1.0e2;
  const double ft8 = 5.0e8;
  const double ot5 = 1.0e5;
  const double ttm3 = 2.0e-3;
  const double fac1 = 743.4;
  const double fac2 = 1.44e-7;
  const double c = clight;
  // const double hbar = 1.0545718176461565e-34;
  const double electron_volt_joule_relationship = 1.602176634e-19;

  // For tailcut we need local radiation integrals
  double *radint = RadiationDampingElement(l, bx, by, dx, dpx, dy, dpy, ax, ay,
                                           angle, k1l, k1sl);

  double p0 = en0 * 1.0e9;
  double restE = amass * 1.0e9;
  double particle_radius = charge * charge / aatom * 1.54e-18;
  double CalphaEC = particle_radius * c / (3.0 * restE * restE * restE) *
                    (p0 * p0 * p0 / len);

  // transverse partition numbers
  double jx = 1.0 - radint[3] / radint[0];
  double jy = 1.0 - radint[4] / radint[0];
  double alphax = 2.0 * CalphaEC * radint[0] * jx;
  double alphay = 2.0 * CalphaEC * radint[0] * jy;
  double alphas = 2.0 * CalphaEC * radint[0] * (jx + jy);
  double tauradxbunch = (1.0 / alphax) / gamma;
  double tauradybunch = (1.0 / alphay) / gamma;
  double tauradsbunch = (1.0 / alphas) / gamma;
  double tauradmax = max(max(tauradxbunch, tauradybunch), tauradsbunch);
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // IMPORTANT: HBAR USED HERE IS REDUCED PLANCK CONSTANT IN GEV!!!!
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // const double hbar = 6.582119569e-25; // 1.0545718176461565e-34;

  const double two = 2.0;
  const double eight = 8.0;

  double betar = sqrt(1.0 - 1.0 / (gamma * gamma));
  double qion = fabs(charge);
  double etrans = ft8 * (gamma * en0 - amass) * (ex / bx);
  double tempev = two * etrans;

  //---- Calculate beam volume to get density (in cm**-3).
  double sigxcm = ot2 * sqrt(ex * bx + (dx * dx * sige * sige));
  double sigycm = ot2 * sqrt(ey * by + (dy * dy * sige * sige));
  double sigtcm = ot2 * sigt;

  double vol = eight * sqrt(pi * pi * pi) * sigxcm * sigycm * sigtcm;
  double densty = pnumber / vol;

  //---- Calculate RMAX as smaller of SIGXCM and DEBYE length.
  double debyel = fac1 * sqrt(tempev / densty) / qion;
  double rmax = min(sigxcm, debyel);

  //---- Calculate RMIN as larger of classical distance of closest approach
  //     or quantum mechanical diffraction limit from nuclear radius.
  double rmincl = fac2 * qion * qion / tempev;
  double rminqm = hbar * c * ot5 / (two * sqrt(ttm3 * etrans * amass));
  double rmin;
  if (tauradmax > 0.0) {
    double rmintailcut =
        1.0f / sqrt(pnumber * pi * tauradmax * c * gamma * sqrt(ex / bx));
    rmin = max(max(rmincl, rminqm), rmintailcut);
    // printf("rmax %12.6e   rmin %12.6e rmincl %12.6e rminqm %12.6e rmintailcut
    // "
    //       "%12.6e\n",
    //       rmax, rmin, rmincl, rminqm, rmintailcut);
  } else {
    rmin = max(rmincl, rminqm);
    // printf("rmax %12.6e   rmin %12.6e rmincl %12.6e rminqm %12.6e\n", rmax,
    //       rmin, rmincl, rminqm);
  }
  double coulog = log(rmax / rmin);

  // calculate coulomb log constant pre-factor
  double constt = pnumber * coulog * r0 * r0 * c /
                  (eight * pi * betar * betar * betar * gamma * gamma * gamma *
                   gamma * ex * ey * sige * sigt);

  clog[0] = coulog;
  clog[1] = constt;
}

/*
-------------------------------------------------------------------------------
ORIGINAL AUTHORS : MADX AUTHORS COPYRIGHT CERN
-------------------------------------------------------------------------------
VERSION 1.0 : INTEGRATE COULOMB LOG CALCULATION IN IBS MODELS AT EACH STEP
(BEFORE FIXED VALUE FOR ENTIRE SIM)
AUTHOR    : TOM MERTENS
DATE      : 19/12/2016
COPYRIGHT : CERN / HZB

    DESCRIPTION :
        FUNCTIONS TO CALCULATE COULOMB LOG
        NOTE: DEPENDS ON CURRENT EMIT

    DETAILS:
        TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
            0 -> gamma
            1 -> charge
            3 -> length
            4 -> energy
            5 -> mass
            6 -> q1
            7 -> q2
            8 -> dxrms
            9 -> dyrms
            10-> gammatr

    REF:
        Calculation of Coulomb logarithm (and print)
        based on the formulae in AIP physics vade mecum p.264 (1981)

-------------------------------------------------------------------------------
pnumber (double)        : number of particles in the bunch
ex (double)             : horizontal emittance
ey (double)             : vertical emittance
twissheader (double[])  : twiss header (madx)
sige (double)           : dE / E
sigt (double)           : simga_s [m]
aatom (double)          : atomic mass
printout (_Bool)        : flag to print out data to screen
-------------------------------------------------------------------------------
clog (double*)          : output (0 -> coulomblog, 1 -> IBS constant factor)
with coulomblog (double)    Constant in eq. (IV.9.1), ZAP user's manual.
-------------------------------------------------------------------------------
*/
void CoulombLog(double pnumber, double ex, double ey, double twissheader[10],
                double sige, double sigt, double r0, bool printout,
                double *clog) {
  // const double one   = 1.0f;
  const double two = 2.0;
  // const double four  = 4.0f;
  const double eight = 8.0f;
  const double ot2 = 1.0e2;
  const double ft8 = 5.0e8;
  const double ot5 = 1.0e5;
  const double ttm3 = 2.0e-3;
  const double fac1 = 743.4;
  const double fac2 = 1.44e-7;

  const double c = clight;
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // IMPORTANT: HBAR USED HERE IS REDUCED PLANCK CONSTANT IN GEV!!!!
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // static double output[2];

  double gamma = twissheader[0];
  double charge = twissheader[1];
  double len = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];
  double q1 = twissheader[5];
  double q2 = twissheader[6];
  double dxbar = twissheader[7];
  double dybar = twissheader[8];
  // double gammatr = twissheader[9];

  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;
  double betar = sqrt(1 - 1 / (gamma * gamma));
  double bxbar = len / (2.0 * pi * q1); // avg betax
  double bybar = len / (2.0 * pi * q2); // avgbety

  //---- Calculate transverse temperature as 2*P*X',
  //     i.e., assume the transverse energy is temperature/2.
  double qion = fabs(charge);
  double etrans = ft8 * (gamma * en0 - amass) * (ex / bxbar);
  double tempev = two * etrans;

  //---- Calculate beam volume to get density (in cm**-3).
  double sigxcm = ot2 * sqrt(ex * bxbar + (dxbar * dxbar * sige * sige));
  double sigycm = ot2 * sqrt(ey * bybar + (dybar * dybar * sige * sige));
  double sigtcm = ot2 * sigt;

  double vol = eight * sqrt(pi * pi * pi) * sigxcm * sigycm * sigtcm;
  double densty = pnumber / vol;

  //---- Calculate RMAX as smaller of SIGXCM and DEBYE length.
  double debyel = fac1 * sqrt(tempev / densty) / qion;
  double rmax = min(sigxcm, debyel);

  //---- Calculate RMIN as larger of classical distance of closest approach
  //     or quantum mechanical diffraction limit from nuclear radius.
  double rmincl = fac2 * qion * qion / tempev;
  double rminqm = hbar * c * ot5 / (two * sqrt(ttm3 * etrans * amass));
  double rmin = max(rmincl, rminqm);
  double coulog = log(rmax / rmin);

  double constt = pnumber * coulog * r0 * r0 * c /
                  (eight * pi * betar * betar * betar * gamma * gamma * gamma *
                   gamma * ex * ey * sige * sigt);
  double cbunch = qion * pnumber * ec * betar * c / len;

  if (printout) {
    printf("Radius max          : %.8e\n", rmax);
    printf("Density             : %.8e\n", pow(densty, -1.0 / 3.0));
    printf("CONST               ; %.8e\n", constt);
    printf("ENERGY              ; %.8e GeV\n", en0);
    printf("BETA                ; %.8e\n", betar);
    printf("GAMMA               ; %.8e\n", gamma);
    printf("COULOMB LOG         ; %.8e\n", coulog);
    printf("X-emittance         ; %.8e m*rad\n", ex);
    printf("Y-emittance         ; %.8e m*rad\n", ey);
    printf("Momentum spread     ; %.8e\n", sige);
    printf("Bunch length        ; %.8e m \n", sigt);
    printf("Particles per bunch ; %.8e\n", pnumber);
    printf("Bunch current       ; %.8e A\n", cbunch);
  }

  clog[0] = coulog;
  clog[1] = constt;
}

/*
------------------------------------------------------------------------------
ORIGINAL AUTHORS : MADX AUTHORS COPYRIGHT CERN
------------------------------------------------------------------------------
VERSION 1.0 : INTEGRATE COULOMB LOG CALCULATION IN IBS MODELS AT EACH STEP
(BEFORE FIXED VALUE FOR ENTIRE SIM)
AUTHOR    : TOM MERTENS
DATE      : 19/12/2016
OPYRIGHT : CERN / HZB

    DESCRIPTION :
        FUNCTIONS TO CALCULATE TAIL CUT COULOMB LOG
        NOTE: DEPENDS ON CURRENT EMIT

    DETAILS:
        TAKES AN ARRAY AS INPUT USING THE FOLLOWING INDEX MAPPING
            0 -> gamma
            1 -> charge
            3 -> length
            4 -> energy
            5 -> mass
            6 -> q1
            7 -> q2
            8 -> dxrms
            9 -> dyrms
            10-> gammatr

    REF:
        Calculation of Tail cut Coulomb logarithm (and print)
        based on the formulae in AIP physics vade mecum p.264 (1981)

-------------------------------------------------------------------------------
pnumber (double)        : number of particles in the bunch
ex (double)             : horizontal emittance
ey (double)             : vertical emittance
twissheader (double[])  : twiss header (madx)
sige (double)           : dE / E
sigt (double)           : simga_s [m]
aatom (double)          : atomic mass
printout (_Bool)        : flag to print out data to screen
-------------------------------------------------------------------------------
clog (double*)          : output (0 -> coulomblog, 1 -> IBS constant factor)
with coulomblog (double)    Constant in eq. (IV.9.1), ZAP user's manual.
-------------------------------------------------------------------------------
*/
void TailCutCoulombLog(double pnumber, double ex, double ey,
                       double twissheader[10], double sige, double sigt,
                       double tauradx, double taurady, double taurads,
                       double r0, bool printout, double *clog) {
  // const double one   = 1.0f;
  const double two = 2.0f;
  // const double four  = 4.0f;
  const double eight = 8.0f;
  const double ot2 = 1.0e2;
  const double ft8 = 5.0e8;
  const double ot5 = 1.0e5;
  const double ttm3 = 2.0e-3;
  const double fac1 = 743.4;
  const double fac2 = 1.44e-7;

  const double c = clight;
  // const double pi = 3.141592653589793f;
  // const double ec = 1.602176634e-19;
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // IMPORTANT: HBAR USED HERE IS REDUCED PLANCK CONSTANT IN GEV!!!!
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  const double hbar = 6.582119569e-25; // 1.0545718176461565e-34;

  // static double output[2];

  double gamma = twissheader[0];
  double charge = twissheader[1];
  double len = twissheader[2];
  double en0 = twissheader[3];
  double amass = twissheader[4];
  double q1 = twissheader[5];
  double q2 = twissheader[6];
  double dxbar = twissheader[7];
  double dybar = twissheader[8];
  // double gammatr = twissheader[9];

  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;
  double betar = sqrt(1 - 1 / (gamma * gamma));
  double bxbar = len / (2.0 * pi * q1); // avg betax
  double bybar = len / (2.0 * pi * q2); // avgbety

  double tauradxbunch = tauradx / gamma;
  double tauradybunch = taurady / gamma;
  double tauradsbunch = taurads / gamma;
  double tauradmax = max(max(tauradxbunch, tauradybunch), tauradsbunch);

  //---- Calculate transverse temperature as 2*P*X',
  //     i.e., assume the transverse energy is temperature/2.
  double qion = fabs(charge);
  double etrans = ft8 * (gamma * en0 - amass) * (ex / bxbar);
  double tempev = two * etrans;

  //---- Calculate beam volume to get density (in cm**-3).
  double sigxcm = ot2 * sqrt(ex * bxbar + (dxbar * dxbar * sige * sige));
  double sigycm = ot2 * sqrt(ey * bybar + (dybar * dybar * sige * sige));
  double sigtcm = ot2 * sigt;

  double vol = eight * sqrt(pi * pi * pi) * sigxcm * sigycm * sigtcm;
  double densty = pnumber / vol;

  //---- Calculate RMAX as smaller of SIGXCM and DEBYE length.
  double debyel = fac1 * sqrt(tempev / densty) / qion;
  double rmax = min(sigxcm, debyel);

  //---- Calculate RMIN as larger of classical distance of closest approach
  //     or quantum mechanical diffraction limit from nuclear radius.
  double rmincl = fac2 * qion * qion / tempev;
  double rminqm = hbar * c * ot5 / (two * sqrt(ttm3 * etrans * amass));
  double rmintailcut =
      1.0 / sqrt(pnumber * pi * tauradmax * c * gamma * sqrt(ex / bxbar));
  double rmin = max(max(rmincl, rminqm), rmintailcut);
  double coulog = log(rmax / rmin);

  double constt = pnumber * coulog * r0 * r0 * c /
                  (eight * pi * betar * betar * betar * gamma * gamma * gamma *
                   gamma * ex * ey * sige * sigt);
  double cbunch = qion * pnumber * ec * betar * c / len;

  if (printout) {
    printf("Radius max          : %.8e\n", rmax);
    printf("Density             : %.8e\n", pow(densty, -1.0f / 3.0f));
    printf("CONST               ; %.8e\n", constt);
    printf("ENERGY              ; %.8e GeV\n", en0);
    printf("BETA                ; %.8e\n", betar);
    printf("GAMMA               ; %.8e\n", gamma);
    printf("COULOMB LOG         ; %.8e\n", coulog);
    printf("X-emittance         ; %.8e m*rad\n", ex);
    printf("Y-emittance         ; %.8e m*rad\n", ey);
    printf("Momentum spread     ; %.8e\n", sige);
    printf("Bunch length        ; %.8e m \n", sigt);
    printf("Particles per bunch ; %.8e\n", pnumber);
    printf("Bunch current       ; %.8e A\n", cbunch);
  }

  clog[0] = coulog;
  clog[1] = constt;
}