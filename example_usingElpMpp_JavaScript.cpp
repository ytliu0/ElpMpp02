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

  char outfile[] = "ElpMpp02C.js";
  char funSuffix[] = "_C";
  generate_javascript_code(outfile, corr, AthU, AthV, AthR, tau, paras, coefs_trim, funSuffix);
  char outfmin[] = "ElpMpp02C_min.js";
  generate_javascript_code_min(outfmin, corr, AthU, AthV, AthR, tau, paras, coefs_trim, funSuffix);

  printf("These two files are created: %s, %s\n\n",outfile, outfmin);
  printf("%s is a human-readable js file,\n", outfile);
  printf("whereas %s is the minified version.\n", outfmin);
  printf("These two files provide exactly the same functions,\n");
  printf("but unnecessary white spaces and comments have been removed\n");
  printf("in %s to optimize performance.\n",outfmin);
  return 0;
}
