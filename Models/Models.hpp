void printouts(double output[3]);

double *PiwinskiSmooth(double pnumber, double ex, double ey, double sigs,
                       double dponp, double twiss[5], double r0);

double *PiwinskiLattice(double pnumber, double ex, double ey, double sigs,
                        double dponp, double twissheader[5], int n,
                        double (*twissdata)[4], // shape [4,n]
                        double r0);

double *PiwinskiLatticeModified(double pnumber, double ex, double ey,
                                double sigs, double dponp,
                                double twissheader[5], int n,
                                double (*twissdata)[6], // shape [6,n]
                                double r0);

double *Nagaitsev(double pnumber, double ex, double ey, double sigs,
                  double dponp, double twissheader[5], int n,
                  double (*twissdata)[6], // shape [6,n]
                  double r0);
/*
================================================================================

IBS MODELS USING DEDICATED INTEGRATORS (SEE INTEGRATORS)

================================================================================
*/
double *ibsmadx(double pnumber, double ex, double ey, double sigs, double sige,
                double twissheader[5], int n, double (*twissdata)[9], double r0,
                bool printout);
double *BjorkenMtingwa2(double pnumber, double ex, double ey, double sigs,
                        double dponp, double twissheader[5], int n,
                        double (*twissdata)[6], double r0);

double *BjorkenMtingwa(double pnumber, double ex, double ey, double sigs,
                       double dponp, double twissheader[5], int n,
                       double (*twissdata)[9], double r0);

double *ConteMartini(double pnumber, double ex, double ey, double sigs,
                     double dponp, double twissheader[5], int n,
                     double (*twissdata)[9], double r0);

double *MadxIBS(double pnumber, double ex, double ey, double sigs, double dponp,
                double twissheader[5], int n, double (*twissdata)[9],
                double r0);