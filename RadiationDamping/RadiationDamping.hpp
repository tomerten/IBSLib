void printradint(double out[6]);
double *RadiationDampingApprox(double latticeLength, double gammaTransition,
                               double dipoleBendingRadius, double betax,
                               double betay);
double *RadiationDampingElement(double l, double bx, double by, double dx,
                                double dpx, double dy, double dpy, double ax,
                                double ay, double angle, double k1l,
                                double k1sl);
double *RadiationDampingLattice(int rows, double (*twissdata)[12]);
double *RadiationDampingGrowthRatesAndEquilibriumEmittances(
    double twiss[5], double radiationIntegrals[6], double aatom);
double RadiationLossesPerTurn(double twiss[5], double I2, double aatom);
double *RadiationCriticalEnergy(double rho, double gamma, double omega);