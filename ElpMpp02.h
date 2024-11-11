// Header file for ElpMpp02.cpp
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define PI 3.14159265358979323846

// Arguments for ELP/MPP02 series
struct Elp_args {
  double W1, D, F, L, Lp, zeta, Me, Ve, EM, Ma, Ju, Sa, Ur, Ne;
};

// Factors multiplied by B1-B5 for longitude and latitude
struct Elp_facs {
  double fA, fB1, fB2, fB3, fB4, fB5;
};

// parameters for adjusting the lunar and planetary arguments
struct Elp_paras {
  // parameters adjusted to fit data
  double Dw1_0, Dw2_0, Dw3_0, Deart_0, Dperi, Dw1_1, Dgam, De, Deart_1, Dep,
         Dw2_1, Dw3_1, Dw1_2, Dw1_3, Dw1_4, Dw2_2, Dw2_3, Dw3_2, Dw3_3;
  // parameters derived from the previous parameters
  double Cw2_1, Cw3_1;
};

// Coefficients for the ELP/MPP02 series
struct Elp_coefs {
  // Main problem
  int n_main_long, n_main_lat, n_main_dist;
  int **i_main_long, **i_main_lat, **i_main_dist;
  double *A_main_long, *A_main_lat, *A_main_dist;

  // Perturbation, longitude
  int n_pert_longT0, n_pert_longT1, n_pert_longT2, n_pert_longT3;
  int **i_pert_longT0, **i_pert_longT1, **i_pert_longT2, **i_pert_longT3;
  double *A_pert_longT0, *A_pert_longT1, *A_pert_longT2, *A_pert_longT3;
  double *ph_pert_longT0, *ph_pert_longT1, *ph_pert_longT2, *ph_pert_longT3;

  // Perturbation, latitude
  int n_pert_latT0, n_pert_latT1, n_pert_latT2;
  int **i_pert_latT0, **i_pert_latT1, **i_pert_latT2;
  double *A_pert_latT0, *A_pert_latT1, *A_pert_latT2;
  double *ph_pert_latT0, *ph_pert_latT1, *ph_pert_latT2;

  // Perturbation, distance
  int n_pert_distT0, n_pert_distT1, n_pert_distT2, n_pert_distT3;
  int **i_pert_distT0, **i_pert_distT1, **i_pert_distT2, **i_pert_distT3;
  double *A_pert_distT0, *A_pert_distT1, *A_pert_distT2, *A_pert_distT3;
  double *ph_pert_distT0, *ph_pert_distT1, *ph_pert_distT2, *ph_pert_distT3;
};

// Calculate the Moon's geocentric X,Y,Z coordinates with respect to
// J2000.0 mean ecliptic and equinox.
// T is the TDB Julian century from J2000.0 = (TBD JD - 2451545)/36525
void getX2000(double T,  Elp_paras &paras, Elp_coefs &coefs,
              double &X, double &Y, double &Z);

// Calculate the Moon's geocentric X,Y,Z coordinates and their time derivatives
// with respect to J2000.0 mean ecliptic and equinox.
// T is the TDB Julian century from J2000.0 = (TBD JD - 2451545)/36525
// Xvec[0] = X, Xvec[1] = Y, Xvec[2] = Z, Xvec[3] = Xdot, Xvec[4] = Ydot, Xvec[5] = Zdot
void getX2000_Xdot2000(double T,  Elp_paras &paras, Elp_coefs &coefs,
              double (&Xvec)[6]);

// restrict x to [-pi, pi)
double mod2pi(double x);

// Set up adjustable parameters
// corr=0: fit to LLR data, corr=1: fit to DE405
void setup_parameters(int corr, Elp_paras &paras, Elp_facs &facs);

// Read main problem file
// n is the number of terms in the series, which is stored in the
//   first line of the data file
void read_main_problem_file(const char *infile, int &n, int ** &i_main, double * &A_main,
                            double fA, Elp_facs facs);

// Read perturbation file
// n is the number of terms in the series, which is stored in the
//   first line of the data file
void read_perturbation_file(const char *infile, int &n, int ** &i_pert, double * &A_pert,
                            double * &phase);

// set up coefficients for the ELP/MPP02 series
void setup_Elp_coefs(Elp_coefs &coefs, Elp_facs facs);

// Compute the lunar and planetary arguments used in the ELP/MPP02 series
void compute_Elp_arguments(double T, Elp_paras paras, Elp_args &args);

// Compute the time derivatives of the lunar and planetary arguments used in the ELP/MPP02 series
// return the time derivatives in units of rad/day
void compute_Elp_arguments_dot(double T, Elp_paras paras, Elp_args &args_dot);

// Compute the time derivatives of the lunar and planetary arguments used in the ELP/MPP02 series
// return the time derivatives in units of rad/day
void compute_Elp_arguments_dot(double T, Elp_paras paras, Elp_args &args_dot);

// Sum the ELP/MPP02 series for the main problem
// dist = 0: sine series; dist != 0: cosine series
double Elp_main_sum(int n, int ** &i_main, double * &A_main, Elp_args &args, int dist);

// Sum the ELP/MPP02 series for perturbations
double Elp_perturbation_sum(int n, int ** &i_pert, double * &A_pert, double * &ph_pert,
                            Elp_args &args);

// Sum the ELP/MPP02 series and its time derivative for the main problem
// dist = 0: sine series; dist != 0: cosine series
void Elp_main_sum_and_derv(int n, int ** &i_main, double * &A_main, Elp_args &args,
               Elp_args &args_dot, int dist, double &sum, double &sum_dot);

// Sum the ELP/MPP02 series and its time derivative for perturbations
void Elp_perturbation_sum_and_derv(int n, int ** &i_pert, double * &A_pert, double * &ph_pert,
                            Elp_args &args, Elp_args &args_dot, double &sum, double &sum_dot);

