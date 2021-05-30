double *PiwinskiSmooth(double pnumber, double ex, double ey, double sigs,
                       double dponp, double twiss[5], double r0);

double *PiwinskiLattice(double pnumber, double ex, double ey, double sigs,
                        double dponp, double twissheader[5], int n,
                        double (*twissdata)[4], // shape [4,n]
                        double r0);
/*

double *PiwinskiLatticeModified(double pnumber, double ex, double ey,
                                double sigs, double dponp,
                                double twissheader[5], int n,
                                double (*twissdata)[n], // shape [6,n]
                                double r0);

double *Nagaitsev(double pnumber, double ex, double ey, double sigs,
                  double dponp, double twissheader[5], int n,
                  double (*twissdata)[n], // shape [6,n]
                  double r0);
*/