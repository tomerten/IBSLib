#ifndef NUMERIC_FUNCTIONS_HPP
#define NUMERIC_FUNCTIONS_HPP
#include <map>
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;

const double clight = 299792458.0;
const double hbar = 6.582119569e-25; // Reduced Planck Constant in GeV !!!!
const double emass = 0.51099895000e-3;
const double pmass = 0.93827208816;
const double nmass = 0.93956542052;          // GeV CODATA 2018
const double mumass = 0.1056583755;          //     ! GeV CODATA 2018
const double atomicmassunit = 0.93149410242; // GeV scipy constants
const double pi = 3.141592653589793;
const double ec = 1.602176634e-19;

// Classical radius [m]
const double erad = 2.8179403262e-15;

// proton radius
double prad = erad * emass / pmass;

double sigefromsigs(double omega0, double sigs, double qs, double eta);
double sigsfromsige(double sige, double gamma, double gammatr, double omega);

double eta(double gamma, double gammatr);

double fmohl(double a, double b, double q, int n);

double particle_radius(double charge, double aatom);

double BetaRelativisticFromGamma(double gamma);

double rds(double x, double y, double z);

double VeffRFeV(double phi, double charge, int nrf, double harmon[],
                double voltages[]);

double VeffRFeVPrime(double phi, double charge, int nrf, double harmon[],
                     double voltages[]);

double VeffRFeVRadlosses(double phi, double U0, double charge, int nrf,
                         double harmon[], double voltages[]);

double synchronuousphase(double target, double init_phi, double U0,
                         double charge, int nrf, double harmon[],
                         double voltages[], double epsilon);

double VeffRFeVPotentialWellDistortion(double phi, double U0, double charge,
                                       int nrf, double harmon[],
                                       double voltages[], double L, double N,
                                       double sigs, double pc);

double VeffRFeVPotentialWellDistortionPrime(double phi, double U0,
                                            double charge, int nrf,
                                            double harmon[], double voltages[],
                                            double L, double N, double sigs,
                                            double pc);

double synchronuousphasewithPWD(double target, double init_phi, double U0,
                                double charge, int nrf, double harmon[],
                                double voltages[], double L, double N,
                                double sigs, double pc, double epsilon);

double synchrotronTune(double omega0, double U0, double charge, int nrf,
                       double harmon[], double voltages[], double phis,
                       double eta, double pc);

double synchrotronTunePWD(double omega0, double U0, double charge, int nrf,
                          double harmon[], double voltages[], double L,
                          double N, double sigs, double phis, double eta,
                          double pc);
double csige(double v0, double h0, double sigs, double U0, double gamma,
             double gammatr, double pc, double circ, double phis,
             bool printout);
void updateTwiss(map<string, vector<double>> &table);
void printTwissMap(string key, map<string, vector<double>> &table);
#endif