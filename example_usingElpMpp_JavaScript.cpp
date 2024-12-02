#include "ElpMpp_JavaScript.cpp"

int main() {
  const double sec = 3.14159265358979323846/648000;
  Elp_paras paras;
  Elp_facs facs;
  Elp_coefs coefs, coefs_trim;

  int corr = 1; // use parameters fitted to DE405/DE406
  setup_parameters(corr, paras, facs);
  setup_Elp_coefs(coefs, facs);
  int ntot = compute_ntot(coefs); // total number of terms in the original series

  double AthU = 1.e-3*sec, AthV = 1.e-3*sec, AthR = 0.1, tau = 50.0;

  // Compute the trimmed coefficients
  trim_Elp_coefs(coefs, coefs_trim, AthU, AthV, AthR, tau);
  int ntot_trim = compute_ntot(coefs_trim);

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

  char outfile[] = "ElpMpp02C.js", outfile_v[] = "ElpMpp02Cv.js";
  char funSuffix[] = "_C", funSuffix_v[] = "_Cv";
  generate_javascript_code(outfile, corr, AthU, AthV, AthR, tau, paras, coefs_trim, funSuffix);
  generate_javascript_code_with_velocity(outfile_v, corr, AthU, AthV, AthR, tau, paras, coefs_trim, funSuffix_v);
  char outfmin[] = "ElpMpp02C_min.js", outfmin_v[] = "ElpMpp02Cv_min.js";
  generate_javascript_code_min(outfmin, paras, coefs_trim, funSuffix);
  generate_javascript_code_with_velocity_min(outfmin_v, paras, coefs_trim, funSuffix_v);

  printf("These 4 files are created: %s, %s, %s, %s\n", outfile, outfile_v, outfmin, outfmin_v);
  printf("%s and %s are human-readable js files,\n", outfile, outfile_v);
  printf("whereas %s and %s are the minified versions.\n", outfmin, outfmin_v);
  printf("The minified versions provide exactly the same functions,\n");
  printf("but white spaces and comments have been removed to optimize performance.\n");
  printf("%s and %s only provide a function to calculate Moon's position,\n", outfile, outfmin);
  printf("whereas %s and %s provide a function to calculation Moon's position and velocity.\n", outfile_v, outfmin_v);
  printf("Note that the calculation of position and velocity doubles the computation time,\n");
  printf("and should be avoided if velocity is not needed.\n\n");
  return 0;
}
