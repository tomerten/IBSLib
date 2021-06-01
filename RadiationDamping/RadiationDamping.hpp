void printradint(double out[6]);
double *RadiationDampingApprox(double latticeLength, double gammaTransition,
                               double dipoleBendingRadius, double betax,
                               double betay);
double *RadiationDampingLattice(int rows, double (*twissdata)[12]);
double *RadiationDampingGrowthRatesAndEquilibriumEmittances(
    double twiss[5], double radiationIntegrals[6], double aatom);
/*
double twiss[], double radiationIntegrals[6], double aatom);
double RadiationLossesPerTurn(double twiss[], double I2, double aatom);
double *RadiationCriticalEnergy(double rho, double gamma, double omega);
*/