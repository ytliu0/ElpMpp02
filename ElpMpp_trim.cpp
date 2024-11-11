// -------------------------------------------------------------------------
//  Generate C++ code and data files for a truncated ELP/MPP02 series.
//
//  Usage: 
//    1. Set the parameter corr: corr=0 uses parameters fitted
//       to the lunar laser ranging (LLR) observation data, corr=1 uses
//       parameters fitted to JPL's DE405/DE406 ephemerides.
//    2. Call the function setup_parameters() to set up parameters
//       corresponding to the choice of corr. There are two sets of
//       parameters: a) parameters for adjusting the lunar and
//       planetary arguments, stored in the struct Elp_paras;
//       b) parameters for adjusting the coefficeients in
//       the ELP/MPP02 series for the main problem, stored in the struct
//       Elp_facs.
//    3. Call the function setup_Elp_coefs() to set up the coefficients
//       for the ELP/MPP02 series. The coefficients are stored in the struct
//       Elp_coefs.
//    4. Set the parameters AthU, AthV, AthR and tau for the truncated series. 
//       See the pdf documentation for how the truncated series 
//       is created based on these 4 parameters.
//    5. Call the function trim_Elp_coefs() to compute the coefficients 
//       for the truncated ELP/MPP02 series. 
//    6. Call output_data_files() to output the coefficients for the truncated 
//       ELP/MPP02 series to 14 data files.
//    7. Call generate_cpp_code() to generate a C++ code that containing functions
//       to calculate the lunar position using the truncated ELP/MPP02 series. 
//    8. (Optional) Call error_est() to estimate the accuracy of the truncated 
//       series by performing a Monte Carlo simulation.
//
//  See example_usingElpMpp_trim.cpp for an example of using this code.
// -------------------------------------------------------------------------

#include "ElpMpp02.h"
#include <iomanip>

// Parameters providing the statistics of the deviation between the truncated 
// series and the full series.
struct devStats{
  double absmax_long, rms_long;
  double absmax_lat, rms_lat;
  double absmax_dist, rms_dist;
};

// Compute the coefficients of the truncated ELP/MPP02 series for the main problem.
void trim_main_coefficients(int n, int ** &i_full, double * &A_full, double Ath, 
                            int &n_trim, int ** &i_trim, double * &A_trim) {
   int i,k;
   // calculate n_trim:
   n_trim = 0;
   for (i=0; i<n; i++) {
      if (fabs(A_full[i]) > Ath) n_trim++;
   }

   if (n_trim > 0) { 
      // Allocate arrays
      i_trim = new int *[n_trim];
      for (i=0; i<n_trim; i++) i_trim[i] = new int[4];
      A_trim = new double[n_trim];
      // set coefficients 
      int j=0;
      for (i=0; i<n; i++) { 
         if (fabs(A_full[i]) > Ath) {
           for (k=0; k<4; k++) {
              i_trim[j][k] = i_full[i][k];
           }
           A_trim[j] = A_full[i];
           j++;
         }
      }
   }
}

// Same as trim_main_coefficients() but includes parameters for a crude error estimate.
void trim_main_coefficients_errorEst(int n, int ** &i_full, double * &A_full, double Ath,
                            int &n_trim, int ** &i_trim, double * &A_trim, 
                            int &n_drop, int ** &i_drop, double * &A_drop,
                            double &absmax, double &sqsum) {
   int i,k;
   // calculate n_trim and n_drop:
   n_trim = 0; n_drop = 0;
   for (i=0; i<n; i++) {
      if (fabs(A_full[i]) > Ath) { 
        n_trim++; 
      } else {
        n_drop++;
      }
   }

   absmax = 0.0;
   sqsum = 0.0; 
   if (n_trim > 0) {
      // Allocate arrays
      i_trim = new int *[n_trim];
      for (i=0; i<n_trim; i++) i_trim[i] = new int[4];
      A_trim = new double[n_trim];
      // set coefficients
      int j=0;
      for (i=0; i<n; i++) {
         if (fabs(A_full[i]) > Ath) {
           for (k=0; k<4; k++) {
              i_trim[j][k] = i_full[i][k];
           }
           A_trim[j] = A_full[i];
           j++;
         } else {
           absmax += fabs(A_full[i]);
           sqsum += A_full[i]*A_full[i];
         }
      }
   }

   // dropped coefficients
   if (n_drop > 0) {
      // Allocate arrays
      i_drop = new int *[n_drop];
      for (i=0; i<n_drop; i++) i_drop[i] = new int[4];
      A_drop = new double[n_drop];
      // set coefficients
      int j=0;
      for (i=0; i<n; i++) {
         if (fabs(A_full[i]) <= Ath) {
           for (k=0; k<4; k++) {
              i_drop[j][k] = i_full[i][k];
           }
           A_drop[j] = A_full[i];
           j++;
         }
      }
   }
}

// Compute the coefficients of the truncated ELP/MPP02 series for perturbatioins.
void trim_pert_coefficients(int n, int ** &i_full, double * &A_full, double * &ph_full,
                            double Ath, int &n_trim, int ** &i_trim, double * &A_trim, 
                            double * &ph_trim) {
   int i,k;
   // calculate n_trim:
   n_trim = 0;
   for (i=0; i<n; i++) {
      if (fabs(A_full[i]) > Ath) n_trim++;
   }

   if (n_trim > 0) {
      // Allocate arrays
      i_trim = new int *[n_trim];
      for (i=0; i<n_trim; i++) i_trim[i] = new int[13];
      A_trim = new double[n_trim];
      ph_trim = new double[n_trim];
      // set coefficients
      int j=0;
      for (i=0; i<n; i++) {
         if (fabs(A_full[i]) > Ath) {
           for (k=0; k<13; k++) {
              i_trim[j][k] = i_full[i][k];
           }
           A_trim[j] = A_full[i];
           ph_trim[j] = ph_full[i];
           j++;
         }
      }
   }
}

// Same as trim_pert_coefficients() but includes parameters for a crude error estimate.
void trim_pert_coefficients_errorEst(int n, int ** &i_full, double * &A_full, 
                            double * &ph_full, double Ath, 
                            int &n_trim, int ** &i_trim, double * &A_trim, double * &ph_trim,
                            int &n_drop, int ** &i_drop, double * &A_drop, double * &ph_drop,
                            double &absmax, double &sqsum) {
   int i,k;
   // calculate n_trim and n_drop:
   n_trim = 0;
   n_drop = 0;
   for (i=0; i<n; i++) {
      if (fabs(A_full[i]) > Ath) { 
        n_trim++;
      } else {
        n_drop++;
      }
   }

   absmax = 0.0;
   sqsum = 0.0;
   if (n_trim > 0) {
      // Allocate arrays
      i_trim = new int *[n_trim];
      for (i=0; i<n_trim; i++) i_trim[i] = new int[13];
      A_trim = new double[n_trim];
      ph_trim = new double[n_trim];
      // set coefficients
      int j=0;
      for (i=0; i<n; i++) {
         if (fabs(A_full[i]) > Ath) {
           for (k=0; k<13; k++) {
              i_trim[j][k] = i_full[i][k];
           }
           A_trim[j] = A_full[i];
           ph_trim[j] = ph_full[i];
           j++;
         } else {
           absmax += fabs(A_full[i]);
           sqsum += A_full[i]*A_full[i];
         }
      }
   }

   // dropped coefficients 
   if (n_drop > 0) {
      // Allocate arrays
      i_drop = new int *[n_drop];
      for (i=0; i<n_drop; i++) i_drop[i] = new int[13];
      A_drop = new double[n_drop];
      ph_drop = new double[n_drop];
      // set coefficients
      int j=0;
      for (i=0; i<n; i++) {
         if (fabs(A_full[i]) <= Ath) {
           for (k=0; k<13; k++) {
              i_drop[j][k] = i_full[i][k];
           }
           A_drop[j] = A_full[i];
           ph_drop[j] = ph_full[i];
           j++;
         }
      }
   }
}

// Compute the coefficients for the truncated ELP/MPP02 series
void trim_Elp_coefs(Elp_coefs &coefs, Elp_coefs &coefs_trim, 
                    double AthU, double AthV, double AthR, double tau) {
   double tau2 = tau*tau;
   double tau3 = tau*tau2;

   // main problem
   trim_main_coefficients(coefs.n_main_long, coefs.i_main_long, coefs.A_main_long, 
                          AthV, coefs_trim.n_main_long, coefs_trim.i_main_long, 
                          coefs_trim.A_main_long);
   trim_main_coefficients(coefs.n_main_lat, coefs.i_main_lat, coefs.A_main_lat,
                          AthU, coefs_trim.n_main_lat, coefs_trim.i_main_lat,
                          coefs_trim.A_main_lat);
   trim_main_coefficients(coefs.n_main_dist, coefs.i_main_dist, coefs.A_main_dist,
                          AthR, coefs_trim.n_main_dist, coefs_trim.i_main_dist,
                          coefs_trim.A_main_dist);

   // perturbation, longitude
   trim_pert_coefficients(coefs.n_pert_longT0, coefs.i_pert_longT0, 
                          coefs.A_pert_longT0, coefs.ph_pert_longT0, AthV,
                          coefs_trim.n_pert_longT0, coefs_trim.i_pert_longT0,
                          coefs_trim.A_pert_longT0, coefs_trim.ph_pert_longT0);
   trim_pert_coefficients(coefs.n_pert_longT1, coefs.i_pert_longT1,
                          coefs.A_pert_longT1, coefs.ph_pert_longT1, AthV/tau,
                          coefs_trim.n_pert_longT1, coefs_trim.i_pert_longT1,
                          coefs_trim.A_pert_longT1, coefs_trim.ph_pert_longT1);
   trim_pert_coefficients(coefs.n_pert_longT2, coefs.i_pert_longT2,
                          coefs.A_pert_longT2, coefs.ph_pert_longT2, AthV/tau2,
                          coefs_trim.n_pert_longT2, coefs_trim.i_pert_longT2,
                          coefs_trim.A_pert_longT2, coefs_trim.ph_pert_longT2);
   trim_pert_coefficients(coefs.n_pert_longT3, coefs.i_pert_longT3,
                          coefs.A_pert_longT3, coefs.ph_pert_longT3, AthV/tau3,
                          coefs_trim.n_pert_longT3, coefs_trim.i_pert_longT3,
                          coefs_trim.A_pert_longT3, coefs_trim.ph_pert_longT3);

   // perturbation, latitude
   trim_pert_coefficients(coefs.n_pert_latT0, coefs.i_pert_latT0,
                          coefs.A_pert_latT0, coefs.ph_pert_latT0, AthU,
                          coefs_trim.n_pert_latT0, coefs_trim.i_pert_latT0,
                          coefs_trim.A_pert_latT0, coefs_trim.ph_pert_latT0);
   trim_pert_coefficients(coefs.n_pert_latT1, coefs.i_pert_latT1,
                          coefs.A_pert_latT1, coefs.ph_pert_latT1, AthU/tau,
                          coefs_trim.n_pert_latT1, coefs_trim.i_pert_latT1,
                          coefs_trim.A_pert_latT1, coefs_trim.ph_pert_latT1);
   trim_pert_coefficients(coefs.n_pert_latT2, coefs.i_pert_latT2,
                          coefs.A_pert_latT2, coefs.ph_pert_latT2, AthU/tau2,
                          coefs_trim.n_pert_latT2, coefs_trim.i_pert_latT2,
                          coefs_trim.A_pert_latT2, coefs_trim.ph_pert_latT2);

   // perturbation, distance
   trim_pert_coefficients(coefs.n_pert_distT0, coefs.i_pert_distT0,
                          coefs.A_pert_distT0, coefs.ph_pert_distT0, AthR,
                          coefs_trim.n_pert_distT0, coefs_trim.i_pert_distT0,
                          coefs_trim.A_pert_distT0, coefs_trim.ph_pert_distT0);
   trim_pert_coefficients(coefs.n_pert_distT1, coefs.i_pert_distT1,
                          coefs.A_pert_distT1, coefs.ph_pert_distT1, AthR/tau,
                          coefs_trim.n_pert_distT1, coefs_trim.i_pert_distT1,
                          coefs_trim.A_pert_distT1, coefs_trim.ph_pert_distT1);
   trim_pert_coefficients(coefs.n_pert_distT2, coefs.i_pert_distT2,
                          coefs.A_pert_distT2, coefs.ph_pert_distT2, AthR/tau2,
                          coefs_trim.n_pert_distT2, coefs_trim.i_pert_distT2,
                          coefs_trim.A_pert_distT2, coefs_trim.ph_pert_distT2);
   trim_pert_coefficients(coefs.n_pert_distT3, coefs.i_pert_distT3,
                          coefs.A_pert_distT3, coefs.ph_pert_distT3, AthR/tau3,
                          coefs_trim.n_pert_distT3, coefs_trim.i_pert_distT3,
                          coefs_trim.A_pert_distT3, coefs_trim.ph_pert_distT3);
}

// Same as trim_Elp_coefs() but include a crude error estimate.
void trim_Elp_coefs_errorEst(Elp_coefs &coefs, Elp_coefs &coefs_trim, Elp_coefs &coefs_drop,
                    double AthU, double AthV, double AthR, double tau, 
                    double T1, double T2, devStats &stats) {
   double tau2 = tau*tau;
   double tau3 = tau*tau2;
   double Tmax = fabs(T2) > fabs(T1) ? fabs(T2):fabs(T1);
   double f1 = (pow(T2,3) - pow(T1,3))/(6.0*(T2-T1));
   double f2 = 0.1*(pow(T2,5) - pow(T1,5))/(T2-T1);
   double f3 = (pow(T2,7) - pow(T1,7))/(14.0*(T2-T1));
   stats.absmax_long = 0.0;
   stats.absmax_lat  = 0.0;
   stats.absmax_dist = 0.0;
   stats.rms_long = 0.0;
   stats.rms_lat  = 0.0;
   stats.rms_dist = 0.0;

   double absmax, sqsum;
   // main problem
   trim_main_coefficients_errorEst(coefs.n_main_long, coefs.i_main_long, coefs.A_main_long, 
                          AthV, coefs_trim.n_main_long, coefs_trim.i_main_long, 
                          coefs_trim.A_main_long, coefs_drop.n_main_long, 
                          coefs_drop.i_main_long, coefs_drop.A_main_long,
                          absmax, sqsum);
   stats.absmax_long += absmax;
   stats.rms_long += 0.5*sqsum;
   trim_main_coefficients_errorEst(coefs.n_main_lat, coefs.i_main_lat, coefs.A_main_lat,
                          AthU, coefs_trim.n_main_lat, coefs_trim.i_main_lat,
                          coefs_trim.A_main_lat, coefs_drop.n_main_lat, 
                          coefs_drop.i_main_lat, coefs_drop.A_main_lat,
                          absmax, sqsum);
   stats.absmax_lat += absmax;
   stats.rms_lat += 0.5*sqsum;
   trim_main_coefficients_errorEst(coefs.n_main_dist, coefs.i_main_dist, coefs.A_main_dist,
                          AthR, coefs_trim.n_main_dist, coefs_trim.i_main_dist,
                          coefs_trim.A_main_dist, coefs_drop.n_main_dist, 
                          coefs_drop.i_main_dist, coefs_drop.A_main_dist,
                          absmax, sqsum);
   stats.absmax_dist += absmax;
   stats.rms_dist += 0.5*sqsum;

   // perturbation, longitude
   trim_pert_coefficients_errorEst(coefs.n_pert_longT0, coefs.i_pert_longT0, 
                          coefs.A_pert_longT0, coefs.ph_pert_longT0, AthV,
                          coefs_trim.n_pert_longT0, coefs_trim.i_pert_longT0,
                          coefs_trim.A_pert_longT0, coefs_trim.ph_pert_longT0, 
                          coefs_drop.n_pert_longT0, coefs_drop.i_pert_longT0,
                          coefs_drop.A_pert_longT0, coefs_drop.ph_pert_longT0,
                          absmax, sqsum);
   stats.absmax_long += absmax;
   stats.rms_long += 0.5*sqsum;
   trim_pert_coefficients_errorEst(coefs.n_pert_longT1, coefs.i_pert_longT1,
                          coefs.A_pert_longT1, coefs.ph_pert_longT1, AthV/tau,
                          coefs_trim.n_pert_longT1, coefs_trim.i_pert_longT1,
                          coefs_trim.A_pert_longT1, coefs_trim.ph_pert_longT1, 
                          coefs_drop.n_pert_longT1, coefs_drop.i_pert_longT1,
                          coefs_drop.A_pert_longT1, coefs_drop.ph_pert_longT1,
                          absmax, sqsum);
   stats.absmax_long += absmax*Tmax;
   stats.rms_long += sqsum*f1;
   trim_pert_coefficients_errorEst(coefs.n_pert_longT2, coefs.i_pert_longT2,
                          coefs.A_pert_longT2, coefs.ph_pert_longT2, AthV/tau2,
                          coefs_trim.n_pert_longT2, coefs_trim.i_pert_longT2,
                          coefs_trim.A_pert_longT2, coefs_trim.ph_pert_longT2, 
                          coefs_drop.n_pert_longT2, coefs_drop.i_pert_longT2,
                          coefs_drop.A_pert_longT2, coefs_drop.ph_pert_longT2,
                          absmax, sqsum);
   stats.absmax_long += absmax*Tmax*Tmax;
   stats.rms_long += sqsum*f2;
   trim_pert_coefficients_errorEst(coefs.n_pert_longT3, coefs.i_pert_longT3,
                          coefs.A_pert_longT3, coefs.ph_pert_longT3, AthV/tau3,
                          coefs_trim.n_pert_longT3, coefs_trim.i_pert_longT3,
                          coefs_trim.A_pert_longT3, coefs_trim.ph_pert_longT3, 
                          coefs_drop.n_pert_longT3, coefs_drop.i_pert_longT3,
                          coefs_drop.A_pert_longT3, coefs_drop.ph_pert_longT3,
                          absmax, sqsum);
   stats.absmax_long += absmax*Tmax*Tmax*Tmax;
   stats.rms_long += sqsum*f3;

   // perturbation, latitude
   trim_pert_coefficients_errorEst(coefs.n_pert_latT0, coefs.i_pert_latT0,
                          coefs.A_pert_latT0, coefs.ph_pert_latT0, AthU,
                          coefs_trim.n_pert_latT0, coefs_trim.i_pert_latT0,
                          coefs_trim.A_pert_latT0, coefs_trim.ph_pert_latT0, 
                          coefs_drop.n_pert_latT0, coefs_drop.i_pert_latT0,
                          coefs_drop.A_pert_latT0, coefs_drop.ph_pert_latT0,
                          absmax, sqsum);
   stats.absmax_lat += absmax;
   stats.rms_lat += 0.5*sqsum;
   trim_pert_coefficients_errorEst(coefs.n_pert_latT1, coefs.i_pert_latT1,
                          coefs.A_pert_latT1, coefs.ph_pert_latT1, AthU/tau,
                          coefs_trim.n_pert_latT1, coefs_trim.i_pert_latT1,
                          coefs_trim.A_pert_latT1, coefs_trim.ph_pert_latT1, 
                          coefs_drop.n_pert_latT1, coefs_drop.i_pert_latT1,
                          coefs_drop.A_pert_latT1, coefs_drop.ph_pert_latT1,
                          absmax, sqsum);
   stats.absmax_lat += absmax*Tmax;
   stats.rms_lat += sqsum*tau2;
   trim_pert_coefficients_errorEst(coefs.n_pert_latT2, coefs.i_pert_latT2,
                          coefs.A_pert_latT2, coefs.ph_pert_latT2, AthU/tau2,
                          coefs_trim.n_pert_latT2, coefs_trim.i_pert_latT2,
                          coefs_trim.A_pert_latT2, coefs_trim.ph_pert_latT2, 
                          coefs_drop.n_pert_latT2, coefs_drop.i_pert_latT2,
                          coefs_drop.A_pert_latT2, coefs_drop.ph_pert_latT2,
                          absmax, sqsum);
   stats.absmax_lat += absmax*Tmax*Tmax;
   stats.rms_lat += sqsum*f2;

   // perturbation, distance
   trim_pert_coefficients_errorEst(coefs.n_pert_distT0, coefs.i_pert_distT0,
                          coefs.A_pert_distT0, coefs.ph_pert_distT0, AthR,
                          coefs_trim.n_pert_distT0, coefs_trim.i_pert_distT0,
                          coefs_trim.A_pert_distT0, coefs_trim.ph_pert_distT0, 
                          coefs_drop.n_pert_distT0, coefs_drop.i_pert_distT0,
                          coefs_drop.A_pert_distT0, coefs_drop.ph_pert_distT0,
                          absmax, sqsum);
   stats.absmax_dist += absmax;
   stats.rms_dist += sqsum*0.5;
   trim_pert_coefficients_errorEst(coefs.n_pert_distT1, coefs.i_pert_distT1,
                          coefs.A_pert_distT1, coefs.ph_pert_distT1, AthR/tau,
                          coefs_trim.n_pert_distT1, coefs_trim.i_pert_distT1,
                          coefs_trim.A_pert_distT1, coefs_trim.ph_pert_distT1, 
                          coefs_drop.n_pert_distT1, coefs_drop.i_pert_distT1,
                          coefs_drop.A_pert_distT1, coefs_drop.ph_pert_distT1,
                          absmax, sqsum);
   stats.absmax_dist += absmax*Tmax;
   stats.rms_dist += sqsum*f1;
   trim_pert_coefficients_errorEst(coefs.n_pert_distT2, coefs.i_pert_distT2,
                          coefs.A_pert_distT2, coefs.ph_pert_distT2, AthR/tau2,
                          coefs_trim.n_pert_distT2, coefs_trim.i_pert_distT2,
                          coefs_trim.A_pert_distT2, coefs_trim.ph_pert_distT2, 
                          coefs_drop.n_pert_distT2, coefs_drop.i_pert_distT2,
                          coefs_drop.A_pert_distT2, coefs_drop.ph_pert_distT2,
                          absmax, sqsum);
   stats.absmax_dist += absmax*Tmax*Tmax;
   stats.rms_dist += sqsum*f2;
   trim_pert_coefficients_errorEst(coefs.n_pert_distT3, coefs.i_pert_distT3,
                          coefs.A_pert_distT3, coefs.ph_pert_distT3, AthR/tau3,
                          coefs_trim.n_pert_distT3, coefs_trim.i_pert_distT3,
                          coefs_trim.A_pert_distT3, coefs_trim.ph_pert_distT3, 
                          coefs_drop.n_pert_distT3, coefs_drop.i_pert_distT3,
                          coefs_drop.A_pert_distT3, coefs_drop.ph_pert_distT3,
                          absmax, sqsum);
   stats.absmax_dist += absmax*Tmax*Tmax*Tmax;
   stats.rms_dist += sqsum*f3;

   stats.rms_long = sqrt(stats.rms_long);
   stats.rms_lat  = sqrt(stats.rms_lat);
   stats.rms_dist = sqrt(stats.rms_dist);
}

// Compute the total number of terms in an ELP/MPP series
int compute_ntot(Elp_coefs &coefs) {
  int ntot = coefs.n_main_long + coefs.n_main_lat + coefs.n_main_dist +
         coefs.n_pert_longT0 + coefs.n_pert_longT1 + coefs.n_pert_longT2 +
         coefs.n_pert_longT3 + coefs.n_pert_latT0 + coefs.n_pert_latT1 +
         coefs.n_pert_latT2 + coefs.n_pert_distT0 + coefs.n_pert_distT1 +
         coefs.n_pert_distT2 + coefs.n_pert_distT3;
  return ntot;
}

// Calculate the differences in Moon's geocentric ecliptic longituide, latitude and distance
// between the full and truncated series by summing the dropped ELP/MPP02 terms. 
// T is the TDB Julian century from J2000.0 = (TBD JD - 2451545)/36525
void calculate_dVUr(double T,  Elp_paras &paras, Elp_coefs &coefs_drop,
              double &dV, double &dU, double &dr) {
  double T2 = T*T;
  double T3 = T*T2;
  Elp_args args;
  compute_Elp_arguments(T, paras, args);

  // Sum the dropped ELP/MPP02 series
  // main problem series
  double main_long = Elp_main_sum(coefs_drop.n_main_long, coefs_drop.i_main_long,
                                  coefs_drop.A_main_long, args, 0);
  double main_lat = Elp_main_sum(coefs_drop.n_main_lat, coefs_drop.i_main_lat,
                                  coefs_drop.A_main_lat, args, 0);
  double main_dist = Elp_main_sum(coefs_drop.n_main_dist, coefs_drop.i_main_dist,
                                  coefs_drop.A_main_dist, args, 1);
  // perturbation, longitude
  double pert_longT0 = Elp_perturbation_sum(coefs_drop.n_pert_longT0, coefs_drop.i_pert_longT0,
                                            coefs_drop.A_pert_longT0, coefs_drop.ph_pert_longT0, args);
  double pert_longT1 = Elp_perturbation_sum(coefs_drop.n_pert_longT1, coefs_drop.i_pert_longT1,
                                            coefs_drop.A_pert_longT1, coefs_drop.ph_pert_longT1, args);
  double pert_longT2 = Elp_perturbation_sum(coefs_drop.n_pert_longT2, coefs_drop.i_pert_longT2,
                                            coefs_drop.A_pert_longT2, coefs_drop.ph_pert_longT2, args);
  double pert_longT3 = Elp_perturbation_sum(coefs_drop.n_pert_longT3, coefs_drop.i_pert_longT3,
                                            coefs_drop.A_pert_longT3, coefs_drop.ph_pert_longT3, args);
  // perturbation, latitude
  double pert_latT0 = Elp_perturbation_sum(coefs_drop.n_pert_latT0, coefs_drop.i_pert_latT0,
                                            coefs_drop.A_pert_latT0, coefs_drop.ph_pert_latT0, args);
  double pert_latT1 = Elp_perturbation_sum(coefs_drop.n_pert_latT1, coefs_drop.i_pert_latT1,
                                            coefs_drop.A_pert_latT1, coefs_drop.ph_pert_latT1, args);
  double pert_latT2 = Elp_perturbation_sum(coefs_drop.n_pert_latT2, coefs_drop.i_pert_latT2,
                                            coefs_drop.A_pert_latT2, coefs_drop.ph_pert_latT2, args);
  // perturbation, distance
  double pert_distT0 = Elp_perturbation_sum(coefs_drop.n_pert_distT0, coefs_drop.i_pert_distT0,
                                            coefs_drop.A_pert_distT0, coefs_drop.ph_pert_distT0, args);
  double pert_distT1 = Elp_perturbation_sum(coefs_drop.n_pert_distT1, coefs_drop.i_pert_distT1,
                                            coefs_drop.A_pert_distT1, coefs_drop.ph_pert_distT1, args);
  double pert_distT2 = Elp_perturbation_sum(coefs_drop.n_pert_distT2, coefs_drop.i_pert_distT2,
                                            coefs_drop.A_pert_distT2, coefs_drop.ph_pert_distT2, args);
  double pert_distT3 = Elp_perturbation_sum(coefs_drop.n_pert_distT3, coefs_drop.i_pert_distT3,
                                            coefs_drop.A_pert_distT3, coefs_drop.ph_pert_distT3, args);

  // Moon's longitude, latitude and distance
  dV = main_long + pert_longT0 + mod2pi(pert_longT1*T) +
                 mod2pi(pert_longT2*T2) + mod2pi(pert_longT3*T3);
  dU  = main_lat + pert_latT0 + mod2pi(pert_latT1*T) + mod2pi(pert_latT2*T2);
  const double ra0 = 384747.961370173/384747.980674318;
  dr = ra0*(main_dist +  pert_distT0 + pert_distT1*T + pert_distT2*T2 + pert_distT3*T3);
}

// Estimate the accuracy of the truncated series by comparing it with 
// the untruncated one at n randomly chosen times T1 <= T <= T2.
void error_est(Elp_coefs &coefs_drop, Elp_paras &paras, 
               int n, double T1, double T2, devStats &errStats) {
  errStats.absmax_long = 0.0;
  errStats.rms_long = 0.0;
  errStats.absmax_lat = 0.0;
  errStats.rms_lat = 0.0;
  errStats.absmax_dist = 0.0;
  errStats.rms_dist = 0.0;

  double dlong, dlat, ddist;
  double dT = T2-T1;
  srand(4852618); // random number seed 
  for (int i=0; i<n; i++) {
    double T = rand()*dT/RAND_MAX + T1;
    calculate_dVUr(T, paras, coefs_drop, dlong,dlat,ddist);
    errStats.absmax_long = fabs(dlong) > errStats.absmax_long ?  fabs(dlong):errStats.absmax_long;
    errStats.absmax_lat = fabs(dlat) > errStats.absmax_lat ?  fabs(dlat):errStats.absmax_lat;
    errStats.absmax_dist = fabs(ddist) > errStats.absmax_dist ?  fabs(ddist):errStats.absmax_dist;
    errStats.rms_long += dlong*dlong;
    errStats.rms_lat += dlat*dlat;
    errStats.rms_dist += ddist*ddist;
    if (i % 1000==0) { 
      cout << "  " << setprecision(3) << i*100.0/n << "% done..." << "\r" << flush;
    }
  }
  cout << "                " << "\r" << flush;
  cout << endl << endl;
  errStats.rms_long = sqrt(errStats.rms_long/n);
  errStats.rms_lat = sqrt(errStats.rms_lat/n);
  errStats.rms_dist = sqrt(errStats.rms_dist/n);
}

// Generate a C++ code that computes the Moon's position using the truncated 
// ELP/MPP02 series. The C++ code will be outputted to the file specified by outfile. 
// The parameter dataFileSuffix must be exactly the same as the one in output_data_files().
void generate_cpp_code(const char* outfile, const char *dataFileSuffix, int corr, 
                       double AthU, double AthV, double AthR, double tau, 
                       Elp_paras &paras) {
     ofstream fin(outfile);
     fin << "// ----------------------------------------------------------------" << endl;
     fin << "//  This code computes a truncated ELP/MPP02 series. " << endl;
     fin << "//" << endl;
     fin << "//  ELP/MPP02 is a semi-analytic solution for the lunar motion developed by" << endl;
     fin << "//  J. Chapront and G. Francou in 2002. It is an improvement of the ELP2000-82B" << endl;
     fin << "//  lunar theory." << endl;
     fin << "//" << endl;
     fin << "//  ELP/MPP02 source paper:" << endl;
     fin << "//    The lunar theory ELP revisited. Introduction of new planetary perturbations " << endl;
     fin << "//     by J. Chapront and G. Francou, Astronomy and Astrophysics, v.404, p.735-742 (2003)" << endl;
     fin << "//     http://adsabs.harvard.edu/abs/2003A%26A...404..735C" << endl;
     fin << "//" << endl;
     fin << "//    Original FORTRAN code and data files:" << endl;
     fin << "//    ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/" << endl;
     fin << "//" << endl;
     fin << "//  This code is generated by \"Elp_Mpp_trim.cpp\" using the following parameters: " << endl;
     fin << "//  corr = " << corr << "," << endl;
     fin << "//  AthU = " << setprecision(17) << AthU << "," << endl;
     fin << "//  AthV = " << setprecision(17) << AthV << "," << endl;
     fin << "//  AthR = " << setprecision(17) << AthR << "," << endl;
     fin << "//  tau = " << setprecision(17) << tau << "," << endl;
     fin << "//" << endl;
     fin << "//  It requires the following 14 data files also generated by the code" << endl;
     fin << "//  \"Elp_Mpp_trim.cpp\" using the same parameters: " << endl;
     fin << "//    elp_main.long" << dataFileSuffix << ", elp_main.lat" << dataFileSuffix 
         << ", elp_main.dist" << dataFileSuffix << "," << endl;
     fin << "//    elp_pert.longT0" << dataFileSuffix << ", elp_pert.longT1" << dataFileSuffix
         << ", elp_pert.longT2" << dataFileSuffix << ", elp_pert.longT3" << dataFileSuffix 
         << ", " << endl;
     fin << "//    elp_pert.latT0" << dataFileSuffix << ", elp_pert.latT1" << dataFileSuffix
         << ", elp_pert.latT2" << dataFileSuffix << ", " << endl;
     fin << "//    elp_pert.distT0" << dataFileSuffix << ", elp_pert.distT1" << dataFileSuffix
         << ", elp_pert.distT2" << dataFileSuffix << ", elp_pert.distT3" << dataFileSuffix
         << endl;
     fin << "//" << endl;
     fin << "//  Usage:" << endl;
     fin << "//    1. Call the function setup_Elp_coefs() to set up the coefficients" << endl;
     fin << "//       for the ELP/MPP02 series. The coefficients are stored in the struct" << endl;
     fin << "//       Elp_coefs." << endl;
     fin << "//    2. Call getX2000() to compute the rectangular geocentric coordinates" << endl;
     fin << "//       of the Moon's position with respect to the mean ecliptic and" << endl;
     fin << "//       equinox of J2000.0." << endl;
     fin << "// ---------------------------------------------------------------- " << endl;
     fin << endl;
     fin << "#include <cmath>" << endl;
     fin << "#include <cstdlib>" << endl;
     fin << "#include <iostream>" << endl;
     fin << "#include <fstream>" << endl;
     fin << "#include <string>" << endl;
     fin << "" << endl;
     fin << "using namespace std;" << endl;
     fin << "" << endl;
     fin << "#define PI 3.14159265358979323846" << endl;
     fin << "" << endl;
     fin << "// Arguments for the ELP/MPP02 series" << endl;
     fin << "struct Elp_args {" << endl;
     fin << "  double W1, D, F, L, Lp, zeta, Me, Ve, EM, Ma, Ju, Sa, Ur, Ne;" << endl;
     fin << "};" << endl;
     fin << "" << endl;
     fin << "// coefficients for the ELP/Mpp02 series" << endl;
     fin << "struct Elp_coefs {" << endl;
     fin << "  // Main problem" << endl;
     fin << "  int n_main_long, n_main_lat, n_main_dist;" << endl;
     fin << "  int **i_main_long, **i_main_lat, **i_main_dist;" << endl;
     fin << "  double *A_main_long, *A_main_lat, *A_main_dist;" << endl;
     fin << "" << endl;
     fin << "  // Perturbation, longitude" << endl;
     fin << "  int n_pert_longT0, n_pert_longT1, n_pert_longT2, n_pert_longT3;" << endl;
     fin << "  int **i_pert_longT0, **i_pert_longT1, **i_pert_longT2, **i_pert_longT3;" << endl;
     fin << "  double *A_pert_longT0, *A_pert_longT1, *A_pert_longT2, *A_pert_longT3;" << endl;
     fin << "  double *ph_pert_longT0, *ph_pert_longT1, *ph_pert_longT2, *ph_pert_longT3;" << endl;
     fin << "" << endl;
     fin << "  // Perturbation, latitude" << endl;
     fin << "  int n_pert_latT0, n_pert_latT1, n_pert_latT2;" << endl;
     fin << "  int **i_pert_latT0, **i_pert_latT1, **i_pert_latT2;" << endl;
     fin << "  double *A_pert_latT0, *A_pert_latT1, *A_pert_latT2;" << endl;
     fin << "  double *ph_pert_latT0, *ph_pert_latT1, *ph_pert_latT2;" << endl;
     fin << "" << endl;
     fin << "  // Perturbation, distance" << endl;
     fin << "  int n_pert_distT0, n_pert_distT1, n_pert_distT2, n_pert_distT3;" << endl;
     fin << "  int **i_pert_distT0, **i_pert_distT1, **i_pert_distT2, **i_pert_distT3;" << endl;
     fin << "  double *A_pert_distT0, *A_pert_distT1, *A_pert_distT2, *A_pert_distT3;" << endl;
     fin << "  double *ph_pert_distT0, *ph_pert_distT1, *ph_pert_distT2, *ph_pert_distT3;" << endl;
     fin << "};" << endl;
     fin << "" << endl;
     fin << "// restrict x to [-pi,pi) " << endl;
     fin << "double mod2pi(double x) {" << endl;
     fin << "  const double tpi = 2.0*PI;" << endl;
     fin << "  return x - tpi*floor((x + PI)/tpi);" << endl;
     fin << "}" << endl;
     fin << "" << endl;
     fin << "// Read main problem file" << endl;
     fin << "// n is the number of terms in the series, which is stored in the" << endl;
     fin << "//   first line of the data file" << endl;
     fin << "void read_main_problem_file(const char *infile, int &n, int ** &i_main, double * &A_main) {" << endl;
     fin << "  int i;" << endl;
     fin << "  ifstream file(infile, ios::in);" << endl;
     fin << "  if (!file) {" << endl;
     fin << "    cerr << \"Error in opening \" << infile << endl;" << endl;
     fin << "    exit(1);" << endl;
     fin << "  }" << endl;
     fin << "  file >> n;" << endl;
     fin << "  if (n > 0) {" << endl;
     fin << "     i_main = new int *[n];" << endl;
     fin << "     A_main = new double[n];" << endl;
     fin << "     for (i=0; i<n; i++) i_main[i] = new int[4];" << endl;
     fin << "     for (i=0; i<n; i++) {" << endl;
     fin << "        if (file.eof()) {" << endl;
     fin << "          cerr << \"Reached the end of the file \" << infile " << endl;
     fin << "               << \" before reading all data!\" << endl;" << endl;
     fin << "          exit(1);" << endl;
     fin << "        }" << endl;
     fin << "        file >> i_main[i][0] >> i_main[i][1] >> i_main[i][2] >> i_main[i][3]" << endl;
     fin << "             >> A_main[i];" << endl;
     fin << "     }" << endl;
     fin << "  }" << endl;
     fin << "  file.close();" << endl;
     fin << "}" << endl;
     fin << "" << endl;
     fin << "// Read perturbation file" << endl;
     fin << "// n is the number of terms in the series, which is stored in the" << endl;
     fin << "//   first line of the data file" << endl;
     fin << "void read_perturbation_file(const char *infile, int &n, int ** &i_pert, double * &A_pert," << endl;
     fin << "                            double * &phase) {" << endl;
     fin << "  int i;" << endl;
     fin << "  ifstream file(infile);" << endl;
     fin << "  if (!file) {" << endl;
     fin << "    cerr << \"Error in opening \" << infile << endl;" << endl;
     fin << "    exit(1);" << endl;
     fin << "  }" << endl;
     fin << "  file >> n;" << endl;
     fin << "  if (n > 0) {" << endl;
     fin << "     i_pert = new int *[n];" << endl;
     fin << "     A_pert = new double[n];" << endl;
     fin << "     phase = new double[n];" << endl;
     fin << "     for (i=0; i<n; i++) i_pert[i] = new int[13];" << endl;
     fin << "     for (i=0; i<n; i++) {" << endl;
     fin << "        if (file.eof()) {" << endl;
     fin << "          cerr << \"Reached the end of the file \" << infile" << endl;
     fin << "               << \" before reading all data!\" << endl;" << endl;
     fin << "          exit(1);" << endl;
     fin << "        }" << endl;
     fin << "        file >> i_pert[i][0] >> i_pert[i][1] >> i_pert[i][2] >> i_pert[i][3]" << endl;
     fin << "             >> i_pert[i][4] >> i_pert[i][5] >> i_pert[i][6] >> i_pert[i][7]" << endl;
     fin << "             >> i_pert[i][8] >> i_pert[i][9] >> i_pert[i][10] >> i_pert[i][11]" << endl;
     fin << "             >> i_pert[i][12] >> A_pert[i] >> phase[i];" << endl;
     fin << "     }" << endl;
     fin << "  }" << endl;
     fin << "  file.close();" << endl;
     fin << "}" << endl;
     fin << "" << endl;
     fin << "// set up coefficients for the ELP/MPP02 series" << endl;
     fin << "void setup_Elp_coefs(Elp_coefs &coefs) {" << endl;
     fin << "  string infile;" << endl;
     fin << "" << endl;
     fin << "  // Main problem" << endl;
     fin << "  infile = \"elp_main.long" << dataFileSuffix << "\";" << endl;
     fin << "  read_main_problem_file(infile.c_str(), coefs.n_main_long, coefs.i_main_long," << endl;
     fin << "                         coefs.A_main_long);" << endl;
     fin << "  infile = \"elp_main.lat" << dataFileSuffix << "\";" << endl;
     fin << "  read_main_problem_file(infile.c_str(), coefs.n_main_lat, coefs.i_main_lat," << endl;
     fin << "                         coefs.A_main_lat);" << endl;
     fin << "  infile = \"elp_main.dist" << dataFileSuffix << "\";" << endl;
     fin << "  read_main_problem_file(infile.c_str(), coefs.n_main_dist, coefs.i_main_dist, " << endl;
     fin << "                         coefs.A_main_dist);" << endl;
     fin << "" << endl;
     fin << "  // perturbation, longitude" << endl;
     fin << "  infile = \"elp_pert.longT0" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_longT0, coefs.i_pert_longT0, " << endl;
     fin << "                         coefs.A_pert_longT0, coefs.ph_pert_longT0);" << endl;
     fin << "  infile = \"elp_pert.longT1" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_longT1, coefs.i_pert_longT1," << endl;
     fin << "                         coefs.A_pert_longT1, coefs.ph_pert_longT1);" << endl;
     fin << "  infile = \"elp_pert.longT2" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_longT2, coefs.i_pert_longT2," << endl;
     fin << "                         coefs.A_pert_longT2, coefs.ph_pert_longT2);" << endl;
     fin << "  infile = \"elp_pert.longT3" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_longT3, coefs.i_pert_longT3," << endl;
     fin << "                         coefs.A_pert_longT3, coefs.ph_pert_longT3);" << endl;
     fin << "" << endl;
     fin << "  // perturbation, latitude" << endl;
     fin << "  infile = \"elp_pert.latT0" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_latT0, coefs.i_pert_latT0," << endl;
     fin << "                         coefs.A_pert_latT0, coefs.ph_pert_latT0);" << endl;
     fin << "  infile = \"elp_pert.latT1" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_latT1, coefs.i_pert_latT1," << endl;
     fin << "                         coefs.A_pert_latT1, coefs.ph_pert_latT1);" << endl;
     fin << "  infile = \"elp_pert.latT2" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_latT2, coefs.i_pert_latT2," << endl;
     fin << "                         coefs.A_pert_latT2, coefs.ph_pert_latT2);" << endl;
     fin << "" << endl;
     fin << "  // perturbation, distance" << endl;
     fin << "  infile = \"elp_pert.distT0" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_distT0, coefs.i_pert_distT0," << endl;
     fin << "                         coefs.A_pert_distT0, coefs.ph_pert_distT0);" << endl;
     fin << "  infile = \"elp_pert.distT1" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_distT1, coefs.i_pert_distT1," << endl;
     fin << "                         coefs.A_pert_distT1, coefs.ph_pert_distT1);" << endl;
     fin << "  infile = \"elp_pert.distT2" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_distT2, coefs.i_pert_distT2," << endl;
     fin << "                         coefs.A_pert_distT2, coefs.ph_pert_distT2);" << endl;
     fin << "  infile = \"elp_pert.distT3" << dataFileSuffix << "\";" << endl;
     fin << "  read_perturbation_file(infile.c_str(), coefs.n_pert_distT3, coefs.i_pert_distT3," << endl;
     fin << "                         coefs.A_pert_distT3, coefs.ph_pert_distT3);" << endl;
     fin << "}" << endl;
     fin << "" << endl;
     fin << "// Compute the lunar and planetary arguments used in the ELP/MPP02 series" << endl;
     fin << "void compute_Elp_arguments(double T, Elp_args &args) {" << endl;
     fin << "  double T2 = T*T;" << endl;
     fin << "  double T3 = T*T2;" << endl;
     fin << "  double T4 = T2*T2;" << endl << endl;
     const double deg = PI/180.0; // degrees -> radians
     const double sec = PI/648000.0; // arcsecs -> radians
     
     double w10 = (-142.0 + 18.0/60.0 +(59.95571 + paras.Dw1_0)/3600.0)*deg;
     double w11 = (1732559343.73604 + paras.Dw1_1)*sec;
     double w12 = (-6.8084 + paras.Dw1_2)*sec;
     double w13 = (0.006604 + paras.Dw1_3)*sec;
     double w14 = (-3.169e-5 + paras.Dw1_4)*sec;
     fin << "  double W1 = " << setprecision(17) << w10 << ";" << endl;
     fin << "  W1 += mod2pi(" << setprecision(17) << w11 << "*T);" << endl;
     fin << "  W1 += mod2pi(" << setprecision(17) << w12 << "*T2);" << endl;
     fin << "  W1 += mod2pi(" << setprecision(17) << w13 << "*T3);" << endl;
     fin << "  W1 += mod2pi(" << setprecision(17) << w14 << "*T4);" << endl;

     double w20 = (83.0 + 21.0/60.0 + (11.67475 + paras.Dw2_0)/3600.0)*deg;
     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*sec;
     double w22 = (-38.2631 + paras.Dw2_2)*sec;
     double w23 = (-0.045047+ paras.Dw2_3)*sec;
     double w24 = 0.00021301*sec;
     fin << "  double W2 = " << setprecision(17) << w20 << ";" << endl;
     fin << "  W2 += mod2pi(" << setprecision(17) << w21 << "*T);" << endl;
     fin << "  W2 += mod2pi(" << setprecision(17) << w22 << "*T2);" << endl;
     fin << "  W2 += mod2pi(" << setprecision(17) << w23 << "*T3);" << endl;
     fin << "  W2 += mod2pi(" << setprecision(17) << w24 << "*T4);" << endl;

     double w30 = (125.0 + 2.0/60.0 + (40.39816 + paras.Dw3_0)/3600.0)*deg;
     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*sec;
     double w32 = (6.359 + paras.Dw3_2)*sec;
     double w33 = (0.007625 + paras.Dw3_3)*sec;
     double w34 = -3.586e-5*sec;
     fin << "  double W3 = " << setprecision(17) << w30 << ";" << endl;
     fin << "  W3 += mod2pi(" << setprecision(17) << w31 << "*T);" << endl;
     fin << "  W3 += mod2pi(" << setprecision(17) << w32 << "*T2);" << endl;
     fin << "  W3 += mod2pi(" << setprecision(17) << w33 << "*T3);" << endl;
     fin << "  W3 += mod2pi(" << setprecision(17) << w34 << "*T4);" << endl;

     double Ea0 = (100.0 + 27.0/60.0 + (59.13885 + paras.Deart_0)/3600.0)*deg;
     double Ea1 = (129597742.293 + paras.Deart_1)*sec;
     double Ea2 = -0.0202*sec;
     double Ea3 = 9e-6*sec;
     double Ea4 = 1.5e-7*sec;
     fin << "  double Ea = " << setprecision(17) << Ea0 << ";" << endl;
     fin << "  Ea += mod2pi(" << setprecision(17) << Ea1 << "*T);" << endl;
     fin << "  Ea += mod2pi(" << setprecision(17) << Ea2 << "*T2);" << endl;
     fin << "  Ea += mod2pi(" << setprecision(17) << Ea3 << "*T3);" << endl;
     fin << "  Ea += mod2pi(" << setprecision(17) << Ea4 << "*T4);" << endl;

     double p0 = (102.0 + 56.0/60.0 + (14.45766 + paras.Dperi)/3600.0)*deg;
     double p1 = 1161.24342*sec;
     double p2 = 0.529265*sec;
     double p3 = -1.1814e-4*sec;
     double p4 = 1.1379e-5*sec;
     fin << "  double pomp = " << setprecision(17) << p0 << ";" << endl;
     fin << "  pomp += mod2pi(" << setprecision(17) << p1 << "*T);" << endl;
     fin << "  pomp += mod2pi(" << setprecision(17) << p2 << "*T2);" << endl;
     fin << "  pomp += mod2pi(" << setprecision(17) << p3 << "*T3);" << endl;
     fin << "  pomp += mod2pi(" << setprecision(17) << p4 << "*T4);" << endl;

     double Me = (-108.0 + 15.0/60.0 + 3.216919/3600.0)*deg;
     double Me1 = 538101628.66888*sec;
     double Ve = (-179.0 + 58.0/60.0 + 44.758419/3600.0)*deg;
     double Ve1 = 210664136.45777*sec;
     double EM = (100.0 + 27.0/60.0 + 59.13885/3600.0)*deg;
     double EM1 = 129597742.293*sec;
     double Ma = (-5.0 + 26.0/60.0 + 3.642778/3600.0)*deg;
     double Ma1 = 68905077.65936*sec;
     double Ju = (34.0 + 21.0/60.0 + 5.379392/3600.0)*deg;
     double Ju1 = 10925660.57335*sec;
     double Sa = (50.0 + 4.0/60.0 + 38.902495/3600.0)*deg;
     double Sa1 = 4399609.33632*sec;
     double Ur = (-46.0 + 3.0/60.0 + 4.354234/3600.0)*deg;
     double Ur1 = 1542482.57845*sec;
     double Ne = (-56.0 + 20.0/60.0 + 56.808371/3600.0)*deg;
     double Ne1 = 786547.897*sec;
     fin << "  double Me = " << setprecision(17) << Me << ";" << endl;
     fin << "  Me += mod2pi(" << setprecision(17) << Me1 << "*T);" << endl;
     fin << "  double Ve = " << setprecision(17) << Ve << ";" << endl;
     fin << "  Ve += mod2pi(" << setprecision(17) << Ve1 << "*T);" << endl;
     fin << "  double EM = " << setprecision(17) << EM << ";" << endl;
     fin << "  EM += mod2pi(" << setprecision(17) << EM1 << "*T);" << endl;
     fin << "  double Ma = " << setprecision(17) << Ma << ";" << endl;
     fin << "  Ma += mod2pi(" << setprecision(17) << Ma1 << "*T);" << endl;
     fin << "  double Ju = " << setprecision(17) << Ju << ";" << endl;
     fin << "  Ju += mod2pi(" << setprecision(17) << Ju1 << "*T);" << endl;
     fin << "  double Sa = " << setprecision(17) << Sa << ";" << endl;
     fin << "  Sa += mod2pi(" << setprecision(17) << Sa1 << "*T);" << endl;
     fin << "  double Ur = " << setprecision(17) << Ur << ";" << endl;
     fin << "  Ur += mod2pi(" << setprecision(17) << Ur1 << "*T);" << endl;
     fin << "  double Ne = " << setprecision(17) << Ne << ";" << endl;
     fin << "  Ne += mod2pi(" << setprecision(17) << Ne1 << "*T);" << endl;
     fin << endl;
     fin << "  // Mean longitude of the Moon" << endl;
     fin << "  args.W1 = mod2pi(W1);" << endl;
     fin << "  // Arguments of Delaunay" << endl;
     fin << "  args.D = mod2pi(W1-Ea + PI);" << endl;
     fin << "  args.F = mod2pi(W1-W3);" << endl;
     fin << "  args.L = mod2pi(W1-W2);" << endl;
     fin << "  args.Lp = mod2pi(Ea-pomp);" << endl;
     fin << endl;
     fin << "  //zeta" << endl;
     fin << "  args.zeta = mod2pi(W1 + 0.02438029560881907*T);" << endl;
     fin << endl;
     fin << "  // Planetary arguments (mean longitudes and mean motions)" << endl;
     fin << "  args.Me = mod2pi(Me);" << endl;
     fin << "  args.Ve = mod2pi(Ve);" << endl;
     fin << "  args.EM = mod2pi(EM);" << endl;
     fin << "  args.Ma = mod2pi(Ma);" << endl;
     fin << "  args.Ju = mod2pi(Ju);" << endl;
     fin << "  args.Sa = mod2pi(Sa);" << endl;
     fin << "  args.Ur = mod2pi(Ur);" << endl;
     fin << "  args.Ne = mod2pi(Ne);" << endl;
     fin << "}" << endl;
     fin << "" << endl;
     fin << "// Sum the ELP/MPP02 series for the main problem" << endl;
     fin << "// dist = 0: sine series; dist != 0: cosine series" << endl;
     fin << "double Elp_main_sum(int n, int ** &i_main, double * &A_main, Elp_args &args, int dist) {" << endl;
     fin << "    int i;" << endl;
     fin << "    double sum = 0.0;" << endl;
     fin << "    double phase;" << endl;
     fin << "    if (dist==0) {" << endl;
     fin << "       // sine series " << endl;
     fin << "       for (i=0; i<n; i++) {" << endl;
     fin << "          phase = i_main[i][0]*args.D + i_main[i][1]*args.F + i_main[i][2]*args.L + " << endl;
     fin << "                  i_main[i][3]*args.Lp;" << endl;
     fin << "          sum += A_main[i]*sin(phase);" << endl;
     fin << "       }" << endl;
     fin << "    } else {" << endl;
     fin << "       // cosine series" << endl;
     fin << "       for (i=0; i<n; i++) {" << endl;
     fin << "          phase = i_main[i][0]*args.D + i_main[i][1]*args.F + i_main[i][2]*args.L +" << endl;
     fin << "                  i_main[i][3]*args.Lp;" << endl;
     fin << "          sum += A_main[i]*cos(phase);" << endl;
     fin << "       }" << endl;
     fin << "    }" << endl;
     fin << "    return sum;" << endl;
     fin << "}" << endl;
     fin << "" << endl;
     fin << "// Sum the ELP/MPP02 series for perturbations" << endl;
     fin << "double Elp_perturbation_sum(int n, int ** &i_pert, double * &A_pert, double * &ph_pert," << endl;
     fin << "                            Elp_args &args) {" << endl;
     fin << "    int i;" << endl;
     fin << "    double sum = 0.0;" << endl;
     fin << "    double phase;" << endl;
     fin << "    for (i=0; i<n; i++) {" << endl;
     fin << "       phase = ph_pert[i] + i_pert[i][0]*args.D + i_pert[i][1]*args.F + " << endl;
     fin << "               i_pert[i][2]*args.L + i_pert[i][3]*args.Lp + i_pert[i][4]*args.Me + " << endl;
     fin << "               i_pert[i][5]*args.Ve + i_pert[i][6]*args.EM + i_pert[i][7]*args.Ma + " << endl;
     fin << "               i_pert[i][8]*args.Ju + i_pert[i][9]*args.Sa + i_pert[i][10]*args.Ur +" << endl;
     fin << "               i_pert[i][11]*args.Ne + i_pert[i][12]*args.zeta;" << endl;
     fin << "       sum += A_pert[i]*sin(phase);" << endl;
     fin << "    }" << endl;
     fin << "    return sum;" << endl;
     fin << "}" << endl;
     fin << "" << endl;
     fin << "// Calculate the Moon's geocentric X,Y,Z coordinates with respect to " << endl;
     fin << "// J2000.0 mean ecliptic and equinox." << endl;
     fin << "// T is the TDB Julian century from J2000.0 = (TBD JD - 2451545)/36525" << endl;
     fin << "void getX2000(double T, Elp_coefs &coefs, double &X, double &Y, double &Z) {" << endl;
     fin << "  double T2 = T*T;" << endl;
     fin << "  double T3 = T*T2;" << endl;
     fin << "  double T4 = T2*T2;" << endl;
     fin << "  double T5 = T2*T3;" << endl;
     fin << "  Elp_args args;" << endl;
     fin << "  compute_Elp_arguments(T, args);" << endl;
     fin << "" << endl;
     fin << "  // Sum the ELP/MPP02 series" << endl;
     fin << "  // main problem series" << endl;
     fin << "  double main_long = Elp_main_sum(coefs.n_main_long, coefs.i_main_long, " << endl;
     fin << "                                  coefs.A_main_long, args, 0);" << endl;
     fin << "  double main_lat = Elp_main_sum(coefs.n_main_lat, coefs.i_main_lat," << endl;
     fin << "                                  coefs.A_main_lat, args, 0);" << endl;
     fin << "  double main_dist = Elp_main_sum(coefs.n_main_dist, coefs.i_main_dist," << endl;
     fin << "                                  coefs.A_main_dist, args, 1);" << endl;
     fin << "  // perturbation, longitude" << endl;
     fin << "  double pert_longT0 = Elp_perturbation_sum(coefs.n_pert_longT0, coefs.i_pert_longT0, " << endl;
     fin << "                                            coefs.A_pert_longT0, coefs.ph_pert_longT0, args);" << endl;
     fin << "  double pert_longT1 = Elp_perturbation_sum(coefs.n_pert_longT1, coefs.i_pert_longT1," << endl;
     fin << "                                            coefs.A_pert_longT1, coefs.ph_pert_longT1, args);" << endl;
     fin << "  double pert_longT2 = Elp_perturbation_sum(coefs.n_pert_longT2, coefs.i_pert_longT2," << endl;
     fin << "                                            coefs.A_pert_longT2, coefs.ph_pert_longT2, args);" << endl;
     fin << "  double pert_longT3 = Elp_perturbation_sum(coefs.n_pert_longT3, coefs.i_pert_longT3," << endl;
     fin << "                                            coefs.A_pert_longT3, coefs.ph_pert_longT3, args);" << endl;
     fin << "  // perturbation, latitude" << endl;
     fin << "  double pert_latT0 = Elp_perturbation_sum(coefs.n_pert_latT0, coefs.i_pert_latT0," << endl;
     fin << "                                            coefs.A_pert_latT0, coefs.ph_pert_latT0, args);" << endl;
     fin << "  double pert_latT1 = Elp_perturbation_sum(coefs.n_pert_latT1, coefs.i_pert_latT1," << endl;
     fin << "                                            coefs.A_pert_latT1, coefs.ph_pert_latT1, args);" << endl;
     fin << "  double pert_latT2 = Elp_perturbation_sum(coefs.n_pert_latT2, coefs.i_pert_latT2," << endl;
     fin << "                                            coefs.A_pert_latT2, coefs.ph_pert_latT2, args);" << endl;
     fin << "  // perturbation, distance" << endl;
     fin << "  double pert_distT0 = Elp_perturbation_sum(coefs.n_pert_distT0, coefs.i_pert_distT0," << endl;
     fin << "                                            coefs.A_pert_distT0, coefs.ph_pert_distT0, args);" << endl;
     fin << "  double pert_distT1 = Elp_perturbation_sum(coefs.n_pert_distT1, coefs.i_pert_distT1," << endl;
     fin << "                                            coefs.A_pert_distT1, coefs.ph_pert_distT1, args);" << endl;
     fin << "  double pert_distT2 = Elp_perturbation_sum(coefs.n_pert_distT2, coefs.i_pert_distT2," << endl;
     fin << "                                            coefs.A_pert_distT2, coefs.ph_pert_distT2, args);" << endl;
     fin << "  double pert_distT3 = Elp_perturbation_sum(coefs.n_pert_distT3, coefs.i_pert_distT3," << endl;
     fin << "                                            coefs.A_pert_distT3, coefs.ph_pert_distT3, args);" << endl;
     fin << "" << endl;
     fin << "  // Moon's longitude, latitude and distance" << endl;
     fin << "  double longM = args.W1 + main_long + pert_longT0 + mod2pi(pert_longT1*T) + " << endl;
     fin << "                 mod2pi(pert_longT2*T2) + mod2pi(pert_longT3*T3);" << endl;
     fin << "  double latM  = main_lat + pert_latT0 + mod2pi(pert_latT1*T) + mod2pi(pert_latT2*T2);" << endl;
     fin << "  // ra0 = a0(DE405)/a0(ELP) = 384747.961370173/384747.980674318 = 0.99999994982652041." << endl;
     fin << "  double r = " << setprecision(17) << 384747.961370173/384747.980674318 
         << "*(main_dist +  pert_distT0 + pert_distT1*T + pert_distT2*T2 + pert_distT3*T3);" 
         << endl;
     fin << "  double x0 = r*cos(longM)*cos(latM);" << endl;
     fin << "  double y0 = r*sin(longM)*cos(latM);" << endl;
     fin << "  double z0 = r*sin(latM);" << endl;
     fin << "" << endl;
     fin << "  // Precession matrix" << endl;
     fin << "  double P = 0.10180391e-4*T + 0.47020439e-6*T2 - 0.5417367e-9*T3 " << endl;
     fin << "             - 0.2507948e-11*T4 + 0.463486e-14*T5;" << endl;
     fin << "  double Q = -0.113469002e-3*T + 0.12372674e-6*T2 + 0.12654170e-8*T3 " << endl;
     fin << "             - 0.1371808e-11*T4 - 0.320334e-14*T5;" << endl;
     fin << "  double sq = sqrt(1 - P*P - Q*Q);" << endl;
     fin << "  double p11 = 1 - 2*P*P;" << endl;
     fin << "  double p12 = 2*P*Q;" << endl;
     fin << "  double p13 = 2*P*sq;" << endl;
     fin << "  double p21 = 2*P*Q;" << endl;
     fin << "  double p22 = 1-2*Q*Q;" << endl;
     fin << "  double p23 = -2*Q*sq;" << endl;
     fin << "  double p31 = -2*P*sq;" << endl;
     fin << "  double p32 = 2*Q*sq;" << endl;
     fin << "  double p33 = 1 - 2*P*P - 2*Q*Q;" << endl;
     fin << "" << endl;
     fin << "  // Finally, components of position vector wrt J2000.0 mean ecliptic and equinox" << endl;
     fin << "  X = p11*x0 + p12*y0 + p13*z0;" << endl;
     fin << "  Y = p21*x0 + p22*y0 + p23*z0;" << endl;
     fin << "  Z = p31*x0 + p32*y0 + p33*z0;" << endl;
     fin << "}" << endl;
     fin.close();
}

void write_main_problem_file(const char *outfile, int &n, int ** &i_main, double * &A_main) {
     ofstream fin(outfile);
     fin << n << endl;
     for (int i=0; i<n; i++) {
        fin << i_main[i][0] << "  " << i_main[i][1] << "  " << i_main[i][2] 
            << "  " << i_main[i][3] << "  " << setprecision(17) << A_main[i] << endl;
     }
     fin.close();
}

void write_perturbation_file(const char *outfile, int &n, int ** &i_pert, 
                             double * &A_pert, double * &phase) { 
     ofstream fin(outfile);
     fin << n << endl;
     for (int i=0; i<n; i++) {
        fin << i_pert[i][0] << "  " << i_pert[i][1] << "  " << i_pert[i][2] << "  " 
            << i_pert[i][3] << "  " << i_pert[i][4] << "  " << i_pert[i][5] << "  "
            << i_pert[i][6] << "  " << i_pert[i][7] << "  " << i_pert[i][8] << "  "
            << i_pert[i][9] << "  " << i_pert[i][10] << "  " << i_pert[i][11] << "  "
            << i_pert[i][12] << "  " << setprecision(17) << A_pert[i] << "  " 
            << phase[i] << endl;
     }
     fin.close();
}

// Output the coefficients of the truncated ELP/MPP02 series to 14 data files. 
// All file names have the form elp_MMMM.[X]dataFileSuffix, where MMMM is 
// main or pert, [X] can be long, lat, dist, longT0, longT1, longT2, longT3, 
// latT0, latT1, latT2, distT0, distT1, distT2, distT3. dataFileSuffix are 
// characters appended to the file name. It provides a shorthand notation 
// for the truncated series.
void output_data_files(const char *dataFileSuffix, Elp_coefs &coefs) {
  char outfile[80];

  // Main problem files 
  sprintf(outfile, "elp_main.long%s",dataFileSuffix);
  write_main_problem_file(outfile, coefs.n_main_long, coefs.i_main_long,
                         coefs.A_main_long);
  sprintf(outfile, "elp_main.lat%s",dataFileSuffix);
  write_main_problem_file(outfile, coefs.n_main_lat, coefs.i_main_lat,
                         coefs.A_main_lat);
  sprintf(outfile, "elp_main.dist%s",dataFileSuffix);
  write_main_problem_file(outfile, coefs.n_main_dist, coefs.i_main_dist,
                         coefs.A_main_dist);

  // perturbation, longitude
  sprintf(outfile, "elp_pert.longT0%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_longT0, coefs.i_pert_longT0,
                         coefs.A_pert_longT0, coefs.ph_pert_longT0);
  sprintf(outfile, "elp_pert.longT1%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_longT1, coefs.i_pert_longT1,
                         coefs.A_pert_longT1, coefs.ph_pert_longT1);
  sprintf(outfile, "elp_pert.longT2%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_longT2, coefs.i_pert_longT2,
                         coefs.A_pert_longT2, coefs.ph_pert_longT2);
  sprintf(outfile, "elp_pert.longT3%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_longT3, coefs.i_pert_longT3,
                         coefs.A_pert_longT3, coefs.ph_pert_longT3);

  // perturbation, latitude
  sprintf(outfile, "elp_pert.latT0%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_latT0, coefs.i_pert_latT0,
                         coefs.A_pert_latT0, coefs.ph_pert_latT0);
  sprintf(outfile, "elp_pert.latT1%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_latT1, coefs.i_pert_latT1,
                         coefs.A_pert_latT1, coefs.ph_pert_latT1);
  sprintf(outfile, "elp_pert.latT2%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_latT2, coefs.i_pert_latT2,
                         coefs.A_pert_latT2, coefs.ph_pert_latT2);

  // perturbation, distance
  sprintf(outfile, "elp_pert.distT0%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_distT0, coefs.i_pert_distT0,
                         coefs.A_pert_distT0, coefs.ph_pert_distT0);
  sprintf(outfile, "elp_pert.distT1%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_distT1, coefs.i_pert_distT1,
                         coefs.A_pert_distT1, coefs.ph_pert_distT1);
  sprintf(outfile, "elp_pert.distT2%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_distT2, coefs.i_pert_distT2,
                         coefs.A_pert_distT2, coefs.ph_pert_distT2);
  sprintf(outfile, "elp_pert.distT3%s", dataFileSuffix);
  write_perturbation_file(outfile, coefs.n_pert_distT3, coefs.i_pert_distT3,
                         coefs.A_pert_distT3, coefs.ph_pert_distT3);
}
