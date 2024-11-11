#include "ElpMpp_trim.cpp"

using namespace std;

int main() 
{ 
  const double sec = 3.14159265358979323846/648000;
  Elp_paras paras;
  Elp_facs facs;
  Elp_coefs coefs, coefs_trim, coefs_drop;
  devStats stats;
  int corr = 1; // using parameters fitted to DE 405/406
  setup_parameters(corr, paras, facs);
  setup_Elp_coefs(coefs, facs);
  int ntot = compute_ntot(coefs); // total number of terms in the original series

  // Set up parameters to create a truncated series
  double AthU = 1.e-3*sec, AthV = 1.e-3*sec, AthR = 0.1, tau = 50.0;
  double T1 = -50.0, T2 = 10.0;

  // Compute the trimmed coefficients
  //trim_Elp_coefs(coefs, coefs_trim, AthU, AthV, AthR, tau);
  trim_Elp_coefs_errorEst(coefs, coefs_trim, coefs_drop, AthU, AthV, AthR, tau, T1,T2, stats);

  int ntot_trim = compute_ntot(coefs_trim); // total number of terms in the truncated series

  printf("----------------------------------------------------\n");
  printf("Creating a truncated ELP/MPP02 series using the following parameters:\n");
  printf("corr = %d\n",corr);
  printf("AthU = %e radians = %f arcsecs\n",AthU,AthU/sec);
  printf("AthV = %e radians = %f arcsecs\n",AthV,AthV/sec);
  printf("AthR = %f km\n",AthR);
  printf("tau = %f\n",tau);
  printf("----------------------------------------------------\n\n");
  printf("Total number of terms in the original series: %d\n",ntot);
  printf("Total number of terms in the truncated series: %d\n\n",ntot_trim);

  printf("Estimation of errors based on the amplitudes of the truncated coefficients\n");
  printf("between T = %f and T = %f...\n", T1,T2);
  printf("Maximum possible error in ecliptic longitude = %f arcsecs\n",stats.absmax_long/sec);
  printf("Estimated root mean square error in ecliptic longitude = %f arcsecs\n",stats.rms_long/sec);
  printf("Maximum possible error in ecliptic latitude = %f arcsecs\n",stats.absmax_lat/sec);
  printf("Estimated root mean square error in ecliptic latitude = %f arcsecs\n",stats.rms_lat/sec);  
  printf("Maximum possible error in geocentric distance = %f km\n",stats.absmax_dist);
  printf("Estimated root mean square error in geocentric distance = %f km\n\n",stats.rms_dist);

  char outfile[] = "ElpMpp02C.cpp";
  char dataFileSuffix[] = "C";
  generate_cpp_code(outfile, dataFileSuffix, 
                    corr, AthU, AthV, AthR, tau, paras);
  output_data_files(dataFileSuffix, coefs_trim);
  printf("The following files are created:\n");
  printf("C++ code: %s\n",outfile);
  printf("14 data files:\n");
  printf("elp_main.long%s, elp_main.lat%s, elp_main.dist%s\n",dataFileSuffix,
           dataFileSuffix,dataFileSuffix);
  printf("elp_pert.longT0%s, elp_pert.longT1%s, elp_pert.longT2%s, elp_pert.longT3%s,\n", 
            dataFileSuffix,dataFileSuffix,dataFileSuffix,dataFileSuffix);
  printf("elp_pert.latT0%s, elp_pert.latT1%s, elp_pert.latT2%s,\n",
            dataFileSuffix,dataFileSuffix,dataFileSuffix);
  printf("elp_pert.distT0%s, elp_pert.distT1%s, elp_pert.distT2%s, elp_pert.distT3%s.\n\n",
            dataFileSuffix,dataFileSuffix,dataFileSuffix,dataFileSuffix);

  int n = 10000;
  printf("Estimate the error of the truncated series by comparing it\n");
  printf("to the original series at %d randomly chosen T between %f and %f...\n",n,T1,T2);
  printf("Please wait...\n");

  error_est(coefs_drop, paras, n, T1, T2, stats);

  printf("Maximum deviation found in ecliptic longitude: %f arcsecs\n",stats.absmax_long/sec);
  printf("Root mean square deviation in ecliptic longitude: %f arcsecs\n",stats.rms_long/sec);
  printf("Maximum deviation found in ecliptic latitude: %f arcsecs\n",stats.absmax_lat/sec);
  printf("Root mean square deviation in ecliptic latitude: %f arcsecs\n",stats.rms_lat/sec);
  printf("Maximum deviation found in geocentric distance: %f km\n",stats.absmax_dist);
  printf("Root mean square deviation in geocentric distance: %f km\n",stats.rms_dist);
  return 0;
}
