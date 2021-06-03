#include <vector>
using namespace std;
void ODE(double twiss[], int nrows, double (*twissdata)[12], double harmon[],
         double voltages[], double stepsize, int maxsteps, vector<double> &ex,
         vector<double> &ey, vector<double> &sigs, vector<double> sige);