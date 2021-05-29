void twclog(double pnumber, double bx, double by, double dx, double dy,
            double ex, double ey, double r0, double gamma, double charge,
            double en0, double amass, double sige, double sigt, double *clog);

void CoulombLog(double pnumber, double ex, double ey, double twissheader[10],
                double sige, double sigt, double r0, bool printout,
                double *clog);

void TailCutCoulombLog(double pnumber, double ex, double ey,
                       double twissheader[10], double sige, double sigt,
                       double tauradx, double taurady, double taurads,
                       double r0, bool printout, double *clog);