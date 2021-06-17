#include <algorithm>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

void WriteToFile(string filename, vector<double> &ex, vector<double> &ey,
                 vector<double> &sigs);

void ODE(map<string, double> &twiss, map<string, vector<double>> &twissdata,
         int nrf, double harmon[], double voltages[], double stepsize,
         int maxsteps, vector<double> &ex, vector<double> &ey,
         vector<double> &sigs, vector<double> sige, int model, double pnumber);
