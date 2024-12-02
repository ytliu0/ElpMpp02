// --------------------------------------------------------------------------------
//   Generate a JavaScript code for a truncated ELP/MPP02 series.
//   
//   Usage:
//      1. Set the parameter corr: corr=0 uses parameters fitted
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
//    6. Call the function generate_javascript_code() to generate a Javascript code 
//       that computes the truncated series.
//
//   **Note that no data file will be created. Terms in the ELP/MPP02 series are 
//     written explicitly in the JavaScript code.**
//
//   See example_usingElpMpp_JavaScript.cpp for an example of using this code.
// --------------------------------------------------------------------------------
#include "ElpMpp_trim.cpp"

// write function for the main problem
void write_main_problem_function(ostream &fin, int n, int ** &i_main, double * &A_main, 
              int istart) {
  string strig;
  if (istart==0) {
    // sine series 
    strig = "Math.sin";
  } else {
    // cosine series. The first term is a constant, so starts at i=1
    strig = "Math.cos";
  }
  const char *trig = strig.c_str();
  char args[4][10] = {"args.D", "args.F", "args.L", "args.Lp"};
  for (int i=istart; i<n; i++) {
     int firstnonzero = 1;
     fin << "  s += " << setprecision(17) << A_main[i] << " * " << trig << "(";
     for (int j=0; j<4; j++) {
        int k = i_main[i][j];
        if (k != 0) {
          if (firstnonzero==1) {
            firstnonzero = 0;
            if (k==1) {
               fin << args[j];
            } else if (k==-1) {
               fin << "-" << args[j];
            } else { 
               fin << k << "*" << args[j];
            }
          } else {
            if (k > 0) { 
              fin << " + ";
            } else {
              fin << " - ";
            }
            if (k==1 || k==-1) {
               fin << args[j];
            } else {
               fin << abs(k) << "*" << args[j];
            }
          }
        }
     }
     fin << ");" << endl;
  }
  fin << "  return s;" << endl; 
  fin << "}" << endl << endl;
}

// write function for the main problem and its time derivative
void write_main_problem_and_derv_function(ostream &fin, int n, int ** &i_main, double * &A_main, 
              int istart) {
  string strig, dstrig;
  if (istart==0) {
    // sine series 
    strig = "Math.sin";
    dstrig = "Math.cos";
  } else {
    // cosine series. The first term is a constant, so starts at i=1
    strig = "Math.cos";
    dstrig = "Math.sin";
  }
  string args[4] = {"args.D", "args.F", "args.L", "args.Lp"};
  string args_dot[4] = {"args_dot.D", "args_dot.F", "args_dot.L", "args_dot.Lp"};
  for (int i=istart; i<n; i++) {
     for (int derv=0; derv<2; derv++) {
        int firstnonzero = 1;
        if (derv==0) {
           fin << "  phase = ";
        } else {
           fin << "  phase_dot = ";
        }
        for (int j=0; j<4; j++) {
           int k = i_main[i][j];
           if (k != 0) {
             string arg = (derv==0 ? args[j]:args_dot[j]);
             if (firstnonzero==1) {
               firstnonzero = 0;
               if (k==1) {
                  fin << arg;
               } else if (k==-1) {
                  fin << "-" << arg;
               } else { 
                  fin << k << "*" << arg;
               }
             } else {
               if (k > 0) { 
                 fin << " + ";
               } else {
                 fin << " - ";
               }
               if (k==1 || k==-1) {
                  fin << arg;
               } else {
                  fin << abs(k) << "*" << arg;
               }
             }
           }
        }
        fin << ";" << endl;
     }
     fin << "  s += " << setprecision(17) << A_main[i] << " * " << strig << "(phase);" << endl;
     fin << (istart==0 ? "  sdot += ":"  sdot -= ") << setprecision(17) << A_main[i] << " * " << dstrig << "(phase)*phase_dot;" << endl;
  }
  fin << "  return {sum:s, sum_dot:sdot};" << endl; 
  fin << "}" << endl << endl;
}

// write function for perturbation 
void write_perturbation_function(ostream &fin, int n, int ** &i_pert, double * &A_pert, 
                           double * &ph_pert) {
  fin << "  let s = 0.0;" << endl;

  char args[13][10] = {"args.D", "args.F", "args.L", "args.Lp", "args.Me", "args.Ve", 
                       "args.EM", "args.Ma", "args.Ju", "args.Sa", "args.Ur", "args.Ne", 
                       "args.zeta"};
  for (int i=0; i<n; i++) {
     fin << "  s += " << setprecision(17) << A_pert[i] << " * Math.sin(";
     int firstnonzero = 1;
     for (int j=0; j<13; j++) {
        int k = i_pert[i][j];
        if (k != 0) {
          if (firstnonzero==1) {
            firstnonzero = 0;
            if (k==1) {
               fin << args[j];
            } else if (k==-1) {
               fin << "-" << args[j];
            } else {
               fin << k << "*" << args[j];
            }
          } else {
            if (k > 0) {
              fin << " + ";
            } else {
              fin << " - ";
            }
            if (k==1 || k==-1) {
               fin << args[j];
            } else {
               fin << abs(k) << "*" << args[j];
            }
          }
        }
     }
     if (ph_pert[i] != 0 ) {
       if (ph_pert[i] > 0) {
         fin << " +";
       } else {
         fin << " -";
       }
       fin << fabs(ph_pert[i]);
     }
     fin << ");" << endl;
  }
  fin << "  return s;" << endl;
  fin << "}" << endl << endl;
}

// write function for perturbation and its derivative
void write_perturbation_and_derv_function(ostream &fin, int n, int ** &i_pert, double * &A_pert, 
                           double * &ph_pert) {
  fin << "  let s = 0.0, sdot = 0.0, phase, phase_dot;" << endl;

  string args[13] = {"args.D", "args.F", "args.L", "args.Lp", "args.Me", "args.Ve", 
                       "args.EM", "args.Ma", "args.Ju", "args.Sa", "args.Ur", "args.Ne", 
                       "args.zeta"};
  string args_dot[13] = {"args_dot.D", "args_dot.F", "args_dot.L", "args_dot.Lp", "args_dot.Me", 
                         "args_dot.Ve", "args_dot.EM", "args_dot.Ma", "args_dot.Ju", "args_dot.Sa",
                         "args_dot.Ur", "args_dot.Ne", "args_dot.zeta"};
  for (int i=0; i<n; i++) {
     for (int derv=0; derv < 2; derv++) {
        int firstnonzero = 1;
        if (derv==0) {
           fin << "  phase = ";
        } else {
           fin << "  phase_dot = ";
        }
        for (int j=0; j<13; j++) {
           int k = i_pert[i][j];
           if (k != 0) {
             string arg = (derv==0 ? args[j]:args_dot[j]);
             if (firstnonzero==1) {
               firstnonzero = 0;
               if (k==1) {
                  fin << arg;
               } else if (k==-1) {
                  fin << "-" << arg;
               } else {
                  fin << k << "*" << arg;
               }
             } else {
               if (k > 0) {
                 fin << " + ";
               } else {
                 fin << " - ";
               }
               if (k==1 || k==-1) {
                  fin << arg;
               } else {
                  fin << abs(k) << "*" << arg;
               }
             }
           }
        }
        if (ph_pert[i] != 0 && derv==0) {
          if (ph_pert[i] > 0) {
            fin << " +";
          } else {
            fin << " -";
          }
          fin << fabs(ph_pert[i]);
        }
        if (derv==1 && firstnonzero==1) { 
          // only constant term in the phase, set phase_dot = 0.
          fin << 0;
        }
        fin << ";" << endl;
     }
     fin << "  s += " << setprecision(17) << A_pert[i] << "*Math.sin(phase);" << endl;
     fin << "  sdot += " << setprecision(17) << A_pert[i] << "*Math.cos(phase)*phase_dot;" << endl;
  }
  fin << "  return {sum:s, sum_dot:sdot};" << endl;
  fin << "}" << endl << endl;
}

// Generate JavaScript code: create header and define the function mod2pi
void generate_javascript_code_header(ostream &fin, int corr,
                       double AthU, double AthV, double AthR, double tau,
                       const char* funSuffix, int derv) {
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
     fin << "//  This code is generated by \"ElpMpp_JavaScript.cpp\" using the following parameters: " << endl;
     fin << "//  corr = " << corr << "," << endl;
     fin << "//  AthU = " << setprecision(17) << AthU << "," << endl;
     fin << "//  AthV = " << setprecision(17) << AthV << "," << endl;
     fin << "//  AthR = " << setprecision(17) << AthR << "," << endl;
     fin << "//  tau = " << setprecision(17) << tau << "," << endl;
     fin << "//" << endl;
     if (derv==0) {
       fin << "//  Usage: Simply call the function getX2000" << funSuffix << "() to" << endl;
       fin << "//         compute the rectangular geocentric coordinators of the Moon" << endl;
       fin << "//         with respect to the mean ecliptic and equinox of J2000.0." << endl;
     } else {
       fin << "//  Usage: Simply call the function getX2000_Xdot2000" << funSuffix << "() to" << endl;
       fin << "//         compute the rectangular geocentric position and velocity of the Moon" << endl;
       fin << "//         with respect to the mean ecliptic and equinox of J2000.0." << endl;
     }
     fin << "//" << endl;
     fin << "// **Note: You should use the minified version of the JavaScript code" << endl;
     fin << "//         instead of this file to optimize performance.**" << endl;
     fin << "//" << endl;
     fin << "// ---------------------------------------------------------------- " << endl;
     fin << endl;
     fin << "\"use strict\";" << endl << endl;
     fin << "// restrict x to [-pi,pi) " << endl;
     fin << "let mod2pi" << funSuffix << " = function(x) {" << endl;
     fin << "  return x - 6.283185307179586*Math.floor(0.5*(x*0.3183098861837907 + 1));" << endl;
     fin << "};" << endl << endl;
}

void generate_javascript_code_getX2000(ostream &fin, const char* funSuffix) {
     fin << "// Calculate the Moon's geocentric X,Y,Z coordinates with respect to " << endl;
     fin << "// J2000.0 mean ecliptic and equinox." << endl;
     fin << "function getX2000" << funSuffix << "(T) {" << endl;
     fin << "  let T2 = T*T;" << endl;
     fin << "  let T3 = T*T2;" << endl;
     fin << "  let T4 = T2*T2;" << endl;
     fin << "  let T5 = T2*T3;" << endl << endl;
     fin << "  // Moon's longitude, latitude and distance" << endl;
     fin << "  let args = compute_Elp_arguments" << funSuffix << "(T);" << endl;
     fin << "  let longM = args.W1 + Elp_main_long" << funSuffix << "(args) + Elp_pert_longT0" << funSuffix << "(args) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_longT1" << funSuffix << "(args)*T) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_longT2" << funSuffix << "(args)*T2) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_longT3" << funSuffix << "(args)*T3);" << endl;
     fin << "  let latM =  Elp_main_lat" << funSuffix << "(args) + Elp_pert_latT0" << funSuffix << "(args) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_latT1" << funSuffix << "(args)*T) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_latT2" << funSuffix << "(args)*T2);" << endl;
     fin << "  let r = " << setprecision(17) << 384747.961370173/384747.980674318 << "*" 
         << "(Elp_main_dist" << funSuffix << "(args) + Elp_pert_distT0" << funSuffix << "(args) +" << endl;
     fin << "              Elp_pert_distT1" << funSuffix << "(args)*T +" << endl;
     fin << "              Elp_pert_distT2" << funSuffix << "(args)*T2 +" << endl;
     fin << "              Elp_pert_distT3" << funSuffix << "(args)*T3);" << endl << endl;
     fin << "  let x0 = r*Math.cos(longM)*Math.cos(latM);" << endl;
     fin << "  let y0 = r*Math.sin(longM)*Math.cos(latM);" << endl;
     fin << "  let z0 = r*Math.sin(latM);" << endl << endl;
     fin << "  // Precession matrix" << endl;
     fin << "  let P = 0.10180391e-4*T + 0.47020439e-6*T2 - 0.5417367e-9*T3 " << endl;
     fin << "             - 0.2507948e-11*T4 + 0.463486e-14*T5;" << endl;
     fin << "  let Q = -0.113469002e-3*T + 0.12372674e-6*T2 + 0.12654170e-8*T3 " << endl;
     fin << "             - 0.1371808e-11*T4 - 0.320334e-14*T5;" << endl;
     fin << "  let sq = Math.sqrt(1 - P*P - Q*Q);" << endl;
     fin << "  let p11 = 1 - 2*P*P;" << endl;
     fin << "  let p12 = 2*P*Q;" << endl;
     fin << "  let p13 = 2*P*sq;" << endl;
     fin << "  let p21 = p12;" << endl;
     fin << "  let p22 = 1-2*Q*Q;" << endl;
     fin << "  let p23 = -2*Q*sq;" << endl;
     fin << "  let p31 = -p13;" << endl;
     fin << "  let p32 = -p23;" << endl;
     fin << "  let p33 = 1 - 2*P*P - 2*Q*Q;" << endl << endl;
     fin << "  // Finally, components of position vector wrt J2000.0 mean ecliptic and equinox" << endl;
     fin << "  let X = p11*x0 + p12*y0 + p13*z0;" << endl;
     fin << "  let Y = p21*x0 + p22*y0 + p23*z0;" << endl;
     fin << "  let Z = p31*x0 + p32*y0 + p33*z0;" << endl << endl;
     fin << "  return {X:X, Y:Y, Z:Z, rGeo:r};" << endl;
     fin << "}" << endl << endl;
}

void generate_javascript_code_getX2000_Xdot2000(ostream &fin, const char* funSuffix) {
     fin << "// Calculate the Moon's geocentric X,Y,Z coordinates and velocity with respect to " << endl;
     fin << "// J2000.0 mean ecliptic and equinox." << endl;
     fin << "function getX2000_Xdot2000" << funSuffix << "(T) {" << endl;
     fin << "  let T2 = T*T;" << endl;
     fin << "  let T3 = T*T2;" << endl;
     fin << "  let T4 = T2*T2;" << endl;
     fin << "  let T5 = T2*T3;" << endl;
     fin << "  let fac = " << setprecision(17) << (1.0/36525) << "; // 1/36525" << endl;
     fin << "  let args = compute_Elp_arguments" << funSuffix << "(T);" << endl;
     fin << "  let args_dot = compute_Elp_arguments_dot" << funSuffix << "(T);" << endl << endl;

     fin << "  // Moon's longitude and time derivative" << endl;
     fin << "  let main_long = Elp_main_long_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_longT0 = Elp_pert_longT0_and_derv" << funSuffix << "(args, args_dot);" << endl; 
     fin << "  let pert_longT1 = Elp_pert_longT1_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_longT2 = Elp_pert_longT2_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_longT3 = Elp_pert_longT3_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let longM = args.W1 + main_long.sum + pert_longT0.sum + ";
     fin << "mod2pi" << funSuffix << "(T*pert_longT1.sum) + ";
     fin << "mod2pi" << funSuffix << "(T2*pert_longT2.sum) + ";
     fin << "mod2pi" << funSuffix << "(T3*pert_longT3.sum);" << endl;
     fin << "  let longM_dot = args_dot.W1 + main_long.sum_dot + pert_longT0.sum_dot + T*pert_longT1.sum_dot + T2*pert_longT2.sum_dot + T3*pert_longT3.sum_dot + ";
     fin << "fac*(pert_longT1.sum + 2*T*pert_longT2.sum + 3*T2*pert_longT3.sum);" << endl << endl;

     fin << "  // Moon's latitude and time derivative" << endl;
     fin << "  let main_lat = Elp_main_lat_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_latT0 = Elp_pert_latT0_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_latT1 = Elp_pert_latT1_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_latT2 = Elp_pert_latT2_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let latM = main_lat.sum + pert_latT0.sum + ";
     fin << "mod2pi" << funSuffix << "(T*pert_latT1.sum) + ";
     fin << "mod2pi" << funSuffix << "(T2*pert_latT2.sum);" << endl;
     fin << "  let latM_dot = main_lat.sum_dot + pert_latT0.sum_dot + T*pert_latT1.sum_dot + T2*pert_latT2.sum_dot + ";
     fin << "fac*(pert_latT1.sum + 2*T*pert_latT2.sum);" << endl << endl;

     fin << "  // Moon's distance and time derivative" << endl;
     fin << "  let main_dist = Elp_main_dist_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_distT0 = Elp_pert_distT0_and_derv" << funSuffix << "(args, args_dot);" << endl; 
     fin << "  let pert_distT1 = Elp_pert_distT1_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_distT2 = Elp_pert_distT2_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let pert_distT3 = Elp_pert_distT3_and_derv" << funSuffix << "(args, args_dot);" << endl;
     fin << "  let ra0 = " << setprecision(17) << 384747.961370173/384747.980674318 << ";" << endl;
     fin << "  let r = ra0*(main_dist.sum + pert_distT0.sum + T*pert_distT1.sum + T2*pert_distT2.sum + T3*pert_distT3.sum);" << endl;
     fin << "  let r_dot = ra0*(main_dist.sum_dot + pert_distT0.sum_dot + T*pert_distT1.sum_dot + T2*pert_distT2.sum_dot + T3*pert_distT3.sum_dot + ";
     fin << "fac*(pert_distT1.sum + 2*T*pert_distT2.sum + 3*T2*pert_distT3.sum));" << endl << endl;

     fin << "  // Moon's ecliptic rectangular coordinates and velocity with respect to mean  ecliptic of date" << endl;
     fin << "  let cV = Math.cos(longM), sV = Math.sin(longM), cU = Math.cos(latM), sU = Math.sin(latM);" << endl;
     fin << "  let x0 = r*cV*cU, y0 = r*sV*cU, z0 = r*sU;" << endl;
     fin << "  let x0_dot = r_dot*cV*cU - r*sV*cU*longM_dot - r*cV*sU*latM_dot;" << endl;
     fin << "  let y0_dot = r_dot*sV*cU + r*cV*cU*longM_dot - r*sV*sU*latM_dot;" << endl;
     fin << "  let z0_dot = r_dot*sU + r*cU*latM_dot;" << endl << endl;
     fin << "  // Precession matrix and time derivative" << endl;
     fin << "  let P = 0.10180391e-4*T + 0.47020439e-6*T2 - 0.5417367e-9*T3 - 0.2507948e-11*T4 + 0.463486e-14*T5;" << endl;
     fin << "  let Q = -0.113469002e-3*T + 0.12372674e-6*T2 + 0.12654170e-8*T3 - 0.1371808e-11*T4 - 0.320334e-14*T5;" << endl;
     fin << "  let sq = Math.sqrt(1 - P*P - Q*Q);" << endl;
     fin << "  let p11 = 1 - 2*P*P;" << endl;
     fin << "  let p12 = 2*P*Q;" << endl;
     fin << "  let p13 = 2*P*sq;" << endl;
     fin << "  let p21 = p12;" << endl;
     fin << "  let p22 = 1-2*Q*Q;" << endl;
     fin << "  let p23 = -2*Q*sq;" << endl;
     fin << "  let p31 = -p13;" << endl;
     fin << "  let p32 = -p23;" << endl;
     fin << "  let p33 = 1 - 2*P*P - 2*Q*Q;" << endl << endl;
     fin << "  let P_dot = fac*(0.10180391e-4 + 0.94040878e-6*T - 1.6252101e-9*T2 - 1.0031792e-11*T3 + 2.31743e-14*T4);" << endl;
     fin << "  let Q_dot = fac*(-0.113469002e-3 + 0.24745348e-6*T + 0.3796251e-8*T2 - 0.5487232e-11*T3 - 1.60167e-14*T4);" << endl;
     fin << "  let sq_dot = -(P*P_dot + Q*Q_dot)/sq;" << endl;
     fin << "  let p11_dot = -4*P*P_dot;" << endl;
     fin << "  let p12_dot = 2*(P_dot*Q + P*Q_dot);" << endl;
     fin << "  let p13_dot = 2*(P_dot*sq + P*sq_dot);" << endl;
     fin << "  let p21_dot = p12_dot;" << endl;
     fin << "  let p22_dot = -4*Q*Q_dot;" << endl;
     fin << "  let p23_dot = -2*(Q_dot*sq + Q*sq_dot);" << endl;
     fin << "  let p31_dot = -p13_dot;" << endl;
     fin << "  let p32_dot = -p23_dot;" << endl;
     fin << "  let p33_dot = p11_dot + p22_dot;" << endl << endl;

     fin << "  // Finally, components of position and velocity vector wrt J2000.0 mean ecliptic and equinox" << endl;
     fin << "  let X = p11*x0 + p12*y0 + p13*z0;" << endl;
     fin << "  let Y = p21*x0 + p22*y0 + p23*z0;" << endl;
     fin << "  let Z = p31*x0 + p32*y0 + p33*z0;" << endl;
     fin << "  let Xdot = p11*x0_dot + p12*y0_dot + p13*z0_dot + p11_dot*x0 + p12_dot*y0 + p13_dot*z0;" << endl;
     fin << "  let Ydot = p21*x0_dot + p22*y0_dot + p23*z0_dot + p21_dot*x0 + p22_dot*y0 + p23_dot*z0;" << endl;
     fin << "  let Zdot = p31*x0_dot + p32*y0_dot + p33*z0_dot + p31_dot*x0 + p32_dot*y0 + p33_dot*z0;" << endl;
     fin << "  return {X:X, Y:Y, Z:Z, rGeo:r, Xdot:Xdot, Ydot:Ydot, Zdot:Zdot};" << endl;
     fin << "}" << endl << endl;
}

void generate_javascript_code_compute_Elp_arguments(ostream &fin, const char* funSuffix, 
           Elp_paras &paras) {
     fin << "function compute_Elp_arguments" << funSuffix << "(T) {" << endl;
     fin << "  let T2 = T*T;" << endl;
     fin << "  let T3 = T*T2;" << endl;
     fin << "  let T4 = T2*T2;" << endl << endl;
     const double deg = PI/180.0; // degrees -> radians
     const double sec = PI/648000.0; // arcsecs -> radians
     
     double w10 = (-142.0 + 18.0/60.0 +(59.95571 + paras.Dw1_0)/3600.0)*deg;
     double w11 = (1732559343.73604 + paras.Dw1_1)*sec;
     double w12 = (-6.8084 + paras.Dw1_2)*sec;
     double w13 = (0.006604 + paras.Dw1_3)*sec;
     double w14 = (-3.169e-5 + paras.Dw1_4)*sec;
     fin << "  let W1 = " << setprecision(17) << w10;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w11 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w12 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w13 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w14 << "*T4);" << endl;

     double w20 = (83.0 + 21.0/60.0 + (11.67475 + paras.Dw2_0)/3600.0)*deg;
     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*sec;
     double w22 = (-38.2631 + paras.Dw2_2)*sec;
     double w23 = (-0.045047+ paras.Dw2_3)*sec;
     double w24 = 0.00021301*sec;
     fin << "  let W2 = " << setprecision(17) << w20;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w21 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w22 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w23 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w24 << "*T4);" << endl;

     double w30 = (125.0 + 2.0/60.0 + (40.39816 + paras.Dw3_0)/3600.0)*deg;
     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*sec;
     double w32 = (6.359 + paras.Dw3_2)*sec;
     double w33 = (0.007625 + paras.Dw3_3)*sec;
     double w34 = -3.586e-5*sec;
     fin << "  let W3 = " << setprecision(17) << w30;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w31 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w32 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w33 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w34 << "*T4);" << endl;

     double Ea0 = (100.0 + 27.0/60.0 + (59.13885 + paras.Deart_0)/3600.0)*deg;
     double Ea1 = (129597742.293 + paras.Deart_1)*sec;
     double Ea2 = -0.0202*sec;
     double Ea3 = 9e-6*sec;
     double Ea4 = 1.5e-7*sec;
     fin << "  let Ea = " << setprecision(17) << Ea0;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << Ea1 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << Ea2 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << Ea3 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << Ea4 << "*T4);" << endl;

     double p0 = (102.0 + 56.0/60.0 + (14.45766 + paras.Dperi)/3600.0)*deg;
     double p1 = 1161.24342*sec;
     double p2 = 0.529265*sec;
     double p3 = -1.1814e-4*sec;
     double p4 = 1.1379e-5*sec;
     fin << "  let pomp = " << setprecision(17) << p0;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << p1 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << p2 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << p3 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << p4 << "*T4);" << endl;

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
     fin << "  let Me = " << setprecision(17) << Me;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Me1 << "*T);" << endl;
     fin << "  let Ve = " << setprecision(17) << Ve;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ve1 << "*T);" << endl;
     fin << "  let EM = " << setprecision(17) << EM;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << EM1 << "*T);" << endl;
     fin << "  let Ma = " << setprecision(17) << Ma;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ma1 << "*T);" << endl;
     fin << "  let Ju = " << setprecision(17) << Ju;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ju1 << "*T);" << endl;
     fin << "  let Sa = " << setprecision(17) << Sa;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Sa1 << "*T);" << endl;
     fin << "  let Ur = " << setprecision(17) << Ur;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ur1 << "*T);" << endl;
     fin << "  let Ne = " << setprecision(17) << Ne;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ne1 << "*T);" << endl << endl;
     fin << "  let args = {};" << endl;
     fin << "  // Mean longitude of the Moon" << endl;
     fin << "  args.W1 = mod2pi" << funSuffix << "(W1);" << endl;
     fin << "  // Arguments of Delaunay" << endl;
     fin << "  args.D = mod2pi" << funSuffix << "(W1-Ea + Math.PI);" << endl;
     fin << "  args.F = mod2pi" << funSuffix << "(W1-W3);" << endl;
     fin << "  args.L = mod2pi" << funSuffix << "(W1-W2);" << endl;
     fin << "  args.Lp = mod2pi" << funSuffix << "(Ea-pomp);" << endl;
     fin << endl;
     fin << "  //zeta" << endl;
     fin << "  args.zeta = mod2pi" << funSuffix << "(W1 + 0.02438029560881907*T);" << endl;
     fin << endl;
     fin << "  // Planetary arguments (mean longitudes and mean motions)" << endl;
     fin << "  args.Me = mod2pi" << funSuffix << "(Me);" << endl;
     fin << "  args.Ve = mod2pi" << funSuffix << "(Ve);" << endl;
     fin << "  args.EM = mod2pi" << funSuffix << "(EM);" << endl;
     fin << "  args.Ma = mod2pi" << funSuffix << "(Ma);" << endl;
     fin << "  args.Ju = mod2pi" << funSuffix << "(Ju);" << endl;
     fin << "  args.Sa = mod2pi" << funSuffix << "(Sa);" << endl;
     fin << "  args.Ur = mod2pi" << funSuffix << "(Ur);" << endl;
     fin << "  args.Ne = mod2pi" << funSuffix << "(Ne);" << endl << endl;
     fin << "  return args;" << endl;
     fin << "}" << endl << endl;
}

void generate_javascript_code_compute_Elp_arguments_dot(ostream &fin, const char* funSuffix, 
           Elp_paras &paras) {
     fin << "function compute_Elp_arguments_dot" << funSuffix << "(T) {" << endl;
     fin << "  let T2 = T*T;" << endl;
     fin << "  let T3 = T*T2;" << endl << endl;
     const double fac = PI/648000.0/36525; // arcsecs -> radians/cy
     
     double w11 = (1732559343.73604 + paras.Dw1_1)*fac;
     double w12 = (-6.8084 + paras.Dw1_2)*2*fac;
     double w13 = (0.006604 + paras.Dw1_3)*3*fac;
     double w14 = (-3.169e-5 + paras.Dw1_4)*4*fac;
     fin << "  let W1 = " << setprecision(17) << w11 << w12 << "*T + " 
         << w13 << "*T2" << w14 << "*T3;" << endl;

     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*fac;
     double w22 = (-38.2631 + paras.Dw2_2)*2*fac;
     double w23 = (-0.045047+ paras.Dw2_3)*3*fac;
     double w24 = 0.00021301*4*fac;
     fin << "  let W2 = " << setprecision(17) << w21 << w22 << "*T" << w23 << "*T2 + " 
         << w24 << "*T3;" << endl;

     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*fac;
     double w32 = (6.359 + paras.Dw3_2)*2*fac;
     double w33 = (0.007625 + paras.Dw3_3)*3*fac;
     double w34 = -3.586e-5*4*fac;
     fin << "  let W3 = " << setprecision(17) << w31 << " + " << w32 << "*T + " 
         << w33 << "*T2" << w34 << "*T3;" << endl;

     double Ea1 = (129597742.293 + paras.Deart_1)*fac;
     double Ea2 = -0.0202*2*fac;
     double Ea3 = 9e-6*3*fac;
     double Ea4 = 1.5e-7*4*fac;
     fin << "  let Ea = " << setprecision(17) << Ea1 << Ea2 << "*T + " 
         << Ea3 << "*T2 + " << Ea4 << "*T3;" << endl;

     double p1 = 1161.24342*fac;
     double p2 = 0.529265*2*fac;
     double p3 = -1.1814e-4*3*fac;
     double p4 = 1.1379e-5*4*fac;
     fin << "  let pomp = " << setprecision(17) << p1 << " + " << p2 << "*T" 
         << p3 << "*T2 + " << p4 << "*T3;" << endl << endl;

     fin << "  let args_dot = {};" << endl;
     fin << "  // Mean longitude of the Moon" << endl;
     fin << "  args_dot.W1 = W1; " << endl;
     fin << "  // Arguments of Delaunay" << endl;
     fin << "  args_dot.D = W1 - Ea;" << endl;
     fin << "  args_dot.F = W1 - W3;" << endl;
     fin << "  args_dot.L = W1 - W2;" << endl;
     fin << "  args_dot.Lp = Ea - pomp;" << endl << endl;
     fin << "  //zeta" << endl;
     fin << "  args_dot.zeta = W1 + " << setprecision(17) << 5028.79695*fac << ";" << endl;
     fin << endl;
     fin << "  // Planetary arguments (mean longitudes and mean motions)" << endl;
     fin << "  args_dot.Me = " << setprecision(17) << 538101628.66888*fac << ";" << endl;
     fin << "  args_dot.Ve = " << setprecision(17) << 210664136.45777*fac << ";" << endl;
     fin << "  args_dot.EM = " << setprecision(17) << 129597742.293*fac << ";" << endl;
     fin << "  args_dot.Ma = " << setprecision(17) << 68905077.65936*fac << ";" << endl;
     fin << "  args_dot.Ju = " << setprecision(17) << 10925660.57335*fac << ";" << endl;
     fin << "  args_dot.Sa = " << setprecision(17) << 4399609.33632*fac << ";" << endl;
     fin << "  args_dot.Ur = " << setprecision(17) << 1542482.57845*fac << ";" << endl;
     fin << "  args_dot.Ne = " << setprecision(17) << 786547.897*fac << ";" << endl << endl;
     fin << "  return args_dot;" << endl;
     fin << "}" << endl << endl;
}

// Generate JavaScript code: sum the ELP/MPP02 series: main problem 
void generate_javascript_code_main_problem(ostream &fin, const char* funSuffix, 
                Elp_coefs &coefs) {
     fin << "// Sum the ELP/MPP02 series: main problem, longitude (" 
         << coefs.n_main_long;
     if (coefs.n_main_long > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_main_long" << funSuffix << "(args) {" << endl;
     fin << "  let s = 0.0;" << endl;
     write_main_problem_function(fin, coefs.n_main_long, coefs.i_main_long, 
              coefs.A_main_long, 0);

     fin << "// Sum the ELP/MPP02 series: main problem, latitude (" 
         << coefs.n_main_lat;
     if (coefs.n_main_lat > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_main_lat" << funSuffix << "(args) {" << endl;
     fin << "  let s = 0.0;" << endl;
     write_main_problem_function(fin, coefs.n_main_lat, coefs.i_main_lat,
              coefs.A_main_lat, 0);

     fin << "// Sum the ELP/MPP02 series: main problem, distance (" 
         << coefs.n_main_dist;
     if (coefs.n_main_dist > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_main_dist" << funSuffix << "(args) {" << endl;
     fin << "  let s = " << setprecision(17) << coefs.A_main_dist[0] <<";" << endl;
     write_main_problem_function(fin, coefs.n_main_dist, coefs.i_main_dist,
              coefs.A_main_dist, 1);
}

void generate_javascript_code_main_problem_and_derv(ostream &fin, const char* funSuffix, 
                Elp_coefs &coefs) {
     fin << "// Sum the ELP/MPP02 series: main problem, longitude (" 
         << coefs.n_main_long;
     if (coefs.n_main_long > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_main_long_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     fin << "  let s = 0.0, sdot = 0, phase, phase_dot;" << endl;
     write_main_problem_and_derv_function(fin, coefs.n_main_long, coefs.i_main_long, 
              coefs.A_main_long, 0);

     fin << "// Sum the ELP/MPP02 series: main problem, latitude (" 
         << coefs.n_main_lat;
     if (coefs.n_main_lat > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_main_lat_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     fin << "  let s = 0.0, sdot = 0, phase, phase_dot;" << endl;
     write_main_problem_and_derv_function(fin, coefs.n_main_lat, coefs.i_main_lat,
              coefs.A_main_lat, 0);

     fin << "// Sum the ELP/MPP02 series: main problem, distance (" 
         << coefs.n_main_dist;
     if (coefs.n_main_dist > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_main_dist_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     fin << "  let sdot = 0, phase, phase_dot;" << endl;
     fin << "  let s = " << setprecision(17) << coefs.A_main_dist[0] <<";" << endl;
     write_main_problem_and_derv_function(fin, coefs.n_main_dist, coefs.i_main_dist,
              coefs.A_main_dist, 1);
}

// Generate JavaScript code: sum the perturbation series
void generate_javascript_code_perturbation(ostream &fin, const char* funSuffix, 
          Elp_coefs &coefs) {
     fin << "// Sum the ELP/MPP02 series: perturbation, longitude T^0 (" 
         << coefs.n_pert_longT0;
     if (coefs.n_pert_longT0 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_longT0" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_longT0, coefs.i_pert_longT0, 
           coefs.A_pert_longT0, coefs.ph_pert_longT0);

     fin << "// Sum the ELP/MPP02 series: perturbation, longitude T^1 ("
         << coefs.n_pert_longT1;
     if (coefs.n_pert_longT1 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_longT1" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_longT1, coefs.i_pert_longT1,
           coefs.A_pert_longT1, coefs.ph_pert_longT1);

     fin << "// Sum the ELP/MPP02 series: perturbation, longitude T^2 ("
         << coefs.n_pert_longT2;
     if (coefs.n_pert_longT2 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_longT2" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_longT2, coefs.i_pert_longT2,
           coefs.A_pert_longT2, coefs.ph_pert_longT2);

     fin << "// Sum the ELP/MPP02 series: perturbation, longitude T^3 ("
         << coefs.n_pert_longT3;
     if (coefs.n_pert_longT3 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_longT3" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_longT3, coefs.i_pert_longT3,
           coefs.A_pert_longT3, coefs.ph_pert_longT3);

     fin << "// Sum the ELP/MPP02 series: perturbation, latitude T^0 ("
         << coefs.n_pert_latT0;
     if (coefs.n_pert_latT0 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_latT0" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_latT0, coefs.i_pert_latT0,
           coefs.A_pert_latT0, coefs.ph_pert_latT0);

     fin << "// Sum the ELP/MPP02 series: perturbation, latitude T^1 ("
         << coefs.n_pert_latT1;
     if (coefs.n_pert_latT1 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_latT1" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_latT1, coefs.i_pert_latT1,
           coefs.A_pert_latT1, coefs.ph_pert_latT1);

     fin << "// Sum the ELP/MPP02 series: perturbation, latitude T^2 ("
         << coefs.n_pert_latT2;
     if (coefs.n_pert_latT2 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_latT2" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_latT2, coefs.i_pert_latT2,
           coefs.A_pert_latT2, coefs.ph_pert_latT2);

     fin << "// Sum the ELP/MPP02 series: perturbation, distance T^0 ("
         << coefs.n_pert_distT0;
     if (coefs.n_pert_distT0 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_distT0" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_distT0, coefs.i_pert_distT0,
           coefs.A_pert_distT0, coefs.ph_pert_distT0);

     fin << "// Sum the ELP/MPP02 series: perturbation, distance T^1 ("
         << coefs.n_pert_distT1;
     if (coefs.n_pert_distT1 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_distT1" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_distT1, coefs.i_pert_distT1,
           coefs.A_pert_distT1, coefs.ph_pert_distT1);

     fin << "// Sum the ELP/MPP02 series: perturbation, distance T^2 ("
         << coefs.n_pert_distT2;
     if (coefs.n_pert_distT2 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_distT2" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_distT2, coefs.i_pert_distT2,
           coefs.A_pert_distT2, coefs.ph_pert_distT2);

     fin << "// Sum the ELP/MPP02 series: perturbation, distance T^3 ("
         << coefs.n_pert_distT3;
     if (coefs.n_pert_distT3 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_distT3" << funSuffix << "(args) {" << endl;
     write_perturbation_function(fin, coefs.n_pert_distT3, coefs.i_pert_distT3,
           coefs.A_pert_distT3, coefs.ph_pert_distT3);
}

void generate_javascript_code_perturbation_and_derv(ostream &fin, const char* funSuffix, 
          Elp_coefs &coefs) {
     fin << "// Sum the ELP/MPP02 series: perturbation, longitude T^0 (" 
         << coefs.n_pert_longT0;
     if (coefs.n_pert_longT0 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_longT0_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_longT0, coefs.i_pert_longT0, 
           coefs.A_pert_longT0, coefs.ph_pert_longT0);

     fin << "// Sum the ELP/MPP02 series: perturbation, longitude T^1 ("
         << coefs.n_pert_longT1;
     if (coefs.n_pert_longT1 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_longT1_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_longT1, coefs.i_pert_longT1,
           coefs.A_pert_longT1, coefs.ph_pert_longT1);

     fin << "// Sum the ELP/MPP02 series: perturbation, longitude T^2 ("
         << coefs.n_pert_longT2;
     if (coefs.n_pert_longT2 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_longT2_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_longT2, coefs.i_pert_longT2,
           coefs.A_pert_longT2, coefs.ph_pert_longT2);

     fin << "// Sum the ELP/MPP02 series: perturbation, longitude T^3 ("
         << coefs.n_pert_longT3;
     if (coefs.n_pert_longT3 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_longT3_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_longT3, coefs.i_pert_longT3,
           coefs.A_pert_longT3, coefs.ph_pert_longT3);

     fin << "// Sum the ELP/MPP02 series: perturbation, latitude T^0 ("
         << coefs.n_pert_latT0;
     if (coefs.n_pert_latT0 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_latT0_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_latT0, coefs.i_pert_latT0,
           coefs.A_pert_latT0, coefs.ph_pert_latT0);

     fin << "// Sum the ELP/MPP02 series: perturbation, latitude T^1 ("
         << coefs.n_pert_latT1;
     if (coefs.n_pert_latT1 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_latT1_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_latT1, coefs.i_pert_latT1,
           coefs.A_pert_latT1, coefs.ph_pert_latT1);

     fin << "// Sum the ELP/MPP02 series: perturbation, latitude T^2 ("
         << coefs.n_pert_latT2;
     if (coefs.n_pert_latT2 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_latT2_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_latT2, coefs.i_pert_latT2,
           coefs.A_pert_latT2, coefs.ph_pert_latT2);

     fin << "// Sum the ELP/MPP02 series: perturbation, distance T^0 ("
         << coefs.n_pert_distT0;
     if (coefs.n_pert_distT0 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_distT0_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_distT0, coefs.i_pert_distT0,
           coefs.A_pert_distT0, coefs.ph_pert_distT0);

     fin << "// Sum the ELP/MPP02 series: perturbation, distance T^1 ("
         << coefs.n_pert_distT1;
     if (coefs.n_pert_distT1 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_distT1_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_distT1, coefs.i_pert_distT1,
           coefs.A_pert_distT1, coefs.ph_pert_distT1);

     fin << "// Sum the ELP/MPP02 series: perturbation, distance T^2 ("
         << coefs.n_pert_distT2;
     if (coefs.n_pert_distT2 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_distT2_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_distT2, coefs.i_pert_distT2,
           coefs.A_pert_distT2, coefs.ph_pert_distT2);

     fin << "// Sum the ELP/MPP02 series: perturbation, distance T^3 ("
         << coefs.n_pert_distT3;
     if (coefs.n_pert_distT3 > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_pert_distT3_and_derv" << funSuffix << "(args, args_dot) {" << endl;
     write_perturbation_and_derv_function(fin, coefs.n_pert_distT3, coefs.i_pert_distT3,
           coefs.A_pert_distT3, coefs.ph_pert_distT3);
}

// funSuffix will be added to every JavaScript function. This is useful to 
// prevent conflicts with code that implements other versions of the truncated series.
void generate_javascript_code(const char* outfile, int corr, 
                       double AthU, double AthV, double AthR, double tau, 
                       Elp_paras &paras, Elp_coefs &coefs, const char* funSuffix) {
     ofstream fin(outfile);
     generate_javascript_code_header(fin, corr, AthU, AthV, AthR, tau, funSuffix, 0);
     generate_javascript_code_getX2000(fin, funSuffix);
     generate_javascript_code_compute_Elp_arguments(fin, funSuffix, paras);
     generate_javascript_code_main_problem(fin, funSuffix, coefs);
     generate_javascript_code_perturbation(fin, funSuffix, coefs);
     fin.close();
}

void generate_javascript_code_with_velocity(const char* outfile, int corr,
                       double AthU, double AthV, double AthR, double tau,
                       Elp_paras &paras, Elp_coefs &coefs, const char* funSuffix) { 
     ofstream fin(outfile);
     generate_javascript_code_header(fin, corr, AthU, AthV, AthR, tau, funSuffix, 1);
     generate_javascript_code_getX2000_Xdot2000(fin, funSuffix);
     generate_javascript_code_compute_Elp_arguments(fin, funSuffix, paras);
     generate_javascript_code_compute_Elp_arguments_dot(fin, funSuffix, paras);
     generate_javascript_code_main_problem_and_derv(fin, funSuffix, coefs);
     generate_javascript_code_perturbation_and_derv(fin, funSuffix, coefs);
     fin.close();
}

// ---- Create minified version ------------------
// write function for the main problem
void write_main_problem_function_min(ostream &fin, int n, int ** &i_main, double * &A_main, 
              int istart) {
  string strig;
  if (istart==0) {
    // sine series 
    strig = "Math.sin";
  } else {
    // cosine series. The first term is a constant, so starts at i=1
    strig = "Math.cos";
  }
  const char *trig = strig.c_str();
  char args[4][10] = {"a.a", "a.b", "a.c", "a.d"};
  for (int i=istart; i<n; i++) {
     int firstnonzero = 1;
     fin << "s+=" << setprecision(17) << A_main[i] << "*" << trig << "(";
     for (int j=0; j<4; j++) {
        int k = i_main[i][j];
        if (k != 0) {
          if (firstnonzero==1) {
            firstnonzero = 0;
            if (k==1) {
               fin << args[j];
            } else if (k==-1) {
               fin << "-" << args[j];
            } else { 
               fin << k << "*" << args[j];
            }
          } else {
            if (k > 0) { 
              fin << "+";
            } else {
              fin << "-";
            }
            if (k==1 || k==-1) {
               fin << args[j];
            } else {
               fin << abs(k) << "*" << args[j];
            }
          }
        }
     }
     fin << ");";
  }
  fin << "return s;}"; 
}

void write_main_problem_and_derv_function_min(ostream &fin, int n, int ** &i_main, double * &A_main, 
              int istart) {
  string strig, dstrig;
  if (istart==0) {
    // sine series 
    strig = "Math.sin";
    dstrig = "Math.cos";
  } else {
    // cosine series. The first term is a constant, so starts at i=1
    strig = "Math.cos";
    dstrig = "Math.sin";
  }
  string args[4] = {"a.a", "a.b", "a.c", "a.d"};
  string args_dot[4] = {"b.a", "b.b", "b.c", "b.d"};
  for (int i=istart; i<n; i++) {
     for (int derv=0; derv<2; derv++) {
        int firstnonzero = 1;
        if (derv==0) {
           fin << "p=";
        } else {
           fin << "q=";
        }
        for (int j=0; j<4; j++) {
           int k = i_main[i][j];
           if (k != 0) {
             string arg = (derv==0 ? args[j]:args_dot[j]);
             if (firstnonzero==1) {
               firstnonzero = 0;
               if (k==1) {
                  fin << arg;
               } else if (k==-1) {
                  fin << "-" << arg;
               } else { 
                  fin << k << "*" << arg;
               }
             } else {
               if (k > 0) { 
                 fin << "+";
               } else {
                 fin << "-";
               }
               if (k==1 || k==-1) {
                  fin << arg;
               } else {
                  fin << abs(k) << "*" << arg;
               }
             }
           }
        }
        fin << ";";
     }
     fin << "s+=" << setprecision(17) << A_main[i] << "*" << strig << "(p);";
     fin << (istart==0 ? "v+=":"v-=") << setprecision(17) << A_main[i] << "*" << dstrig << "(p)*q;";
  }
  fin << "return {s:s,v:v};}"; 
}

// write function for perturbation 
void write_perturbation_function_min(ostream &fin, int n, int ** &i_pert, double * &A_pert, 
                           double * &ph_pert) {
  fin << "let s=0;";

  char args[13][10] = {"a.a", "a.b", "a.c", "a.d", "a.e", "a.f", 
                       "a.g", "a.h", "a.i", "a.j", "a.k", "a.l", 
                       "a.m"};
  for (int i=0; i<n; i++) {
     fin << "s+=" << setprecision(17) << A_pert[i] << "*Math.sin(";
     int firstnonzero = 1;
     for (int j=0; j<13; j++) {
        int k = i_pert[i][j];
        if (k != 0) {
          if (firstnonzero==1) {
            firstnonzero = 0;
            if (k==1) {
               fin << args[j];
            } else if (k==-1) {
               fin << "-" << args[j];
            } else {
               fin << k << "*" << args[j];
            }
          } else {
            if (k > 0) {
              fin << "+";
            } else {
              fin << "-";
            }
            if (k==1 || k==-1) {
               fin << args[j];
            } else {
               fin << abs(k) << "*" << args[j];
            }
          }
        }
     }
     if (ph_pert[i] != 0 ) {
       if (ph_pert[i] > 0) {
         fin << "+";
       } else {
         fin << "-";
       }
       fin << fabs(ph_pert[i]);
     }
     fin << ");";
  }
  fin << "return s;}";
}

void write_perturbation_and_derv_function_min(ostream &fin, int n, int ** &i_pert, double * &A_pert, 
                           double * &ph_pert) {
  fin << "let s=0,v=0,p,q;";

  string args[13] = {"a.a", "a.b", "a.c", "a.d", "a.e", "a.f", 
                     "a.g", "a.h", "a.i", "a.j", "a.k", "a.l", 
                     "a.m"};
  string args_dot[13] = {"b.a", "b.b", "b.c", "b.d", "b.e", "b.f",
                         "b.g", "b.h", "b.i", "b.j", "b.k", "b.l",
                         "b.m"};
  for (int i=0; i<n; i++) {
     for (int derv=0; derv < 2; derv++) {
        int firstnonzero = 1;
        if (derv==0) {
           fin << "p=";
        } else {
           fin << "q=";
        }
        for (int j=0; j<13; j++) {
           int k = i_pert[i][j];
           if (k != 0) {
             string arg = (derv==0 ? args[j]:args_dot[j]);
             if (firstnonzero==1) {
               firstnonzero = 0;
               if (k==1) {
                  fin << arg;
               } else if (k==-1) {
                  fin << "-" << arg;
               } else {
                  fin << k << "*" << arg;
               }
             } else {
               if (k > 0) {
                 fin << "+";
               } else {
                 fin << "-";
               }
               if (k==1 || k==-1) {
                  fin << arg;
               } else {
                  fin << abs(k) << "*" << arg;
               }
             }
           }
        }
        if (ph_pert[i] != 0 && derv==0) {
          if (ph_pert[i] > 0) {
            fin << "+";
          } else {
            fin << "-";
          }
          fin << fabs(ph_pert[i]);
        }
        if (derv==1 && firstnonzero==1) {
          // only constant term in the phase, set q = 0.
          fin << 0;
        }
        fin << ";";
     }
     fin << "s+=" << setprecision(17) << A_pert[i] << "*Math.sin(p);";
     fin << "v+=" << setprecision(17) << A_pert[i] << "*Math.cos(p)*q;";
  }
  fin << "return {s:s, v:v};}";
}

void generate_javascript_code_min_getX2000(ostream &fin, const char *funSuffix) {
     fin << "function getX2000" << funSuffix << "(T){";
     fin << "let T2=T*T;";
     fin << "let T3=T*T2;";
     fin << "let T4=T2*T2;";
     fin << "let T5=T2*T3;";
     fin << "let a=compute_Elp_arguments" << funSuffix << "(T);";
     fin << "let U=a.W+Elp_main_long" << funSuffix << "(a)+Elp_pert_longT0" << funSuffix << "(a)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_longT1" << funSuffix << "(a)*T)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_longT2" << funSuffix << "(a)*T2)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_longT3" << funSuffix << "(a)*T3);";
     fin << "let V=Elp_main_lat" << funSuffix << "(a)+Elp_pert_latT0" << funSuffix << "(a)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_latT1" << funSuffix << "(a)*T)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_latT2" << funSuffix << "(a)*T2);";
     fin << "let r=" << setprecision(17) << 384747.961370173/384747.980674318 << "*" 
         << "(Elp_main_dist" << funSuffix << "(a)+Elp_pert_distT0" << funSuffix << "(a)+";
     fin << "Elp_pert_distT1" << funSuffix << "(a)*T+";
     fin << "Elp_pert_distT2" << funSuffix << "(a)*T2+";
     fin << "Elp_pert_distT3" << funSuffix << "(a)*T3);";
     fin << "let x0=r*Math.cos(U)*Math.cos(V);";
     fin << "let y0=r*Math.sin(U)*Math.cos(V);";
     fin << "let z0=r*Math.sin(V);";
     fin << "let P=0.10180391e-4*T+0.47020439e-6*T2-0.5417367e-9*T3";
     fin << "-0.2507948e-11*T4+0.463486e-14*T5;";
     fin << "let Q=-0.113469002e-3*T+0.12372674e-6*T2+0.12654170e-8*T3";
     fin << "-0.1371808e-11*T4-0.320334e-14*T5;";
     fin << "let sq=Math.sqrt(1-P*P-Q*Q);";
     fin << "let p11=1-2*P*P;";
     fin << "let p12=2*P*Q;";
     fin << "let p13=2*P*sq;";
     fin << "let p21=2*P*Q;";
     fin << "let p22=1-2*Q*Q;";
     fin << "let p23=-2*Q*sq;";
     fin << "let p31=-2*P*sq;";
     fin << "let p32=2*Q*sq;";
     fin << "let p33=1-2*P*P-2*Q*Q;";
     fin << "let X=p11*x0+p12*y0+p13*z0;";
     fin << "let Y=p21*x0+p22*y0+p23*z0;";
     fin << "let Z=p31*x0+p32*y0+p33*z0;";
     fin << "return {X:X,Y:Y,Z:Z,rGeo:r};}";
}

void generate_javascript_code_min_getX2000_Xdot2000(ostream &fin, const char *funSuffix) {
     fin << "function getX2000_Xdot2000" << funSuffix << "(T){";
     fin << "let T2=T*T,T3=T*T2,T4=T2*T2,T5=T2*T3,fa=" << setprecision(17) << 1.0/36525 << ";";
     fin << "let a=compute_Elp_arguments" << funSuffix << "(T),";
     fin << "b=compute_Elp_arguments_dot" << funSuffix << "(T);";

     fin << "let d=Elp_main_long_and_derv" << funSuffix << "(a,b),"
         << "e=Elp_pert_longT0_and_derv" << funSuffix << "(a,b),"
         << "f=Elp_pert_longT1_and_derv" << funSuffix << "(a,b),"
         << "g=Elp_pert_longT2_and_derv" << funSuffix << "(a,b),"
         << "h=Elp_pert_longT3_and_derv" << funSuffix << "(a,b);";
     fin << "let V=a.W+d.s+e.s+mod2pi" << funSuffix << "(T*f.s)+mod2pi" 
         << funSuffix << "(T2*g.s)+mod2pi" << funSuffix << "(T3*h.s),"
         << "v=b.W+d.v+e.v+T*f.v+T2*g.v+T3*h.v+fa*(f.s+2*T*g.s+3*T2*h.s);";
     
     fin << "d=Elp_main_lat_and_derv" << funSuffix << "(a,b);"
         << "e=Elp_pert_latT0_and_derv" << funSuffix << "(a,b);"
         << "f=Elp_pert_latT1_and_derv" << funSuffix << "(a,b);"
         << "g=Elp_pert_latT2_and_derv" << funSuffix << "(a,b);";
     fin << "let U=d.s+e.s+mod2pi" << funSuffix << "(T*f.s)+mod2pi" << funSuffix << "(T2*g.s),"
         << "u=d.v+e.v+T*f.v+T2*g.v+fa*(f.s+2*T*g.s);";

     fin << "d=Elp_main_dist_and_derv" << funSuffix << "(a,b);"
         << "e=Elp_pert_distT0_and_derv" << funSuffix << "(a,b);"
         << "f=Elp_pert_distT1_and_derv" << funSuffix << "(a,b);"
         << "g=Elp_pert_distT2_and_derv" << funSuffix << "(a,b);"
         << "h=Elp_pert_distT3_and_derv" << funSuffix << "(a,b);";
     fin << "let ra0=" << setprecision(17) << 384747.961370173/384747.980674318 << ",";
     fin << "r=ra0*(d.s+e.s+T*f.s+T2*g.s+T3*h.s),"
         << "vr=ra0*(d.v+e.v+T*f.v+T2*g.v+T3*h.v+fa*(f.s+2*T*g.s+3*T2*h.s));";

     fin << "let cV=Math.cos(V),sV=Math.sin(V),cU=Math.cos(U),sU=Math.sin(U),"
         << "x0=r*cV*cU,y0=r*sV*cU,z0=r*sU,"
         << "vx0=vr*cV*cU-r*sV*cU*v-r*cV*sU*u,"
         << "vy0=vr*sV*cU+r*cV*cU*v-r*sV*sU*u,"
         << "vz0=vr*sU+r*cU*u;";

     fin << "let P=0.10180391e-4*T+0.47020439e-6*T2-0.5417367e-9*T3";
     fin << "-0.2507948e-11*T4+0.463486e-14*T5,";
     fin << "Q=-0.113469002e-3*T+0.12372674e-6*T2+0.12654170e-8*T3";
     fin << "-0.1371808e-11*T4-0.320334e-14*T5,";
     fin << "sq=Math.sqrt(1-P*P-Q*Q),";
     fin << "p11=1-2*P*P,";
     fin << "p12=2*P*Q,";
     fin << "p13=2*P*sq,";
     fin << "p21=p12,";
     fin << "p22=1-2*Q*Q,";
     fin << "p23=-2*Q*sq,";
     fin << "p31=-p13,";
     fin << "p32=-p23,";
     fin << "p33=1-2*P*P-2*Q*Q;";

     fin << "let Pd=fa*(0.10180391e-4+0.94040878e-6*T-1.6252101e-9*T2-1.0031792e-11*T3+2.31743e-14*T4),"
         << "Qd=fa*(-0.113469002e-3+0.24745348e-6*T+0.3796251e-8*T2-0.5487232e-11*T3-1.60167e-14*T4),"
         << "sqd=-(P*Pd+Q*Qd)/sq,"
         << "q11=-4*P*Pd,q12=2*(Pd*Q+P*Qd),q13=2*(Pd*sq+P*sqd),q21=q12,q22=-4*Q*Qd,"
         << "q23=-2*(Qd*sq+Q*sqd),q31=-q13,q32=-q23,q33=q11+q22;";

     fin << "return {X:p11*x0+p12*y0+p13*z0,Y:p21*x0+p22*y0+p23*z0,Z:p31*x0+p32*y0+p33*z0,rGeo:r,"
         << "Xdot:p11*vx0+p12*vy0+p13*vz0+q11*x0+q12*y0+q13*z0,"
         << "Ydot:p21*vx0+p22*vy0+p23*vz0+q21*x0+q22*y0+q23*z0,"
         << "Zdot:p31*vx0+p32*vy0+p33*vz0+q31*x0+q32*y0+q33*z0};}";
}

void generate_javascript_code_min_compute_Elp_arguments(ostream &fin, const char *funSuffix, 
              Elp_paras &paras) {
     fin << "function compute_Elp_arguments" << funSuffix << "(T){";
     fin << "let T2=T*T,T3=T*T2,T4=T2*T2;";
     const double deg = PI/180.0; // degrees -> radians
     const double sec = PI/648000.0; // arcsecs -> radians
     
     double w10 = (-142.0 + 18.0/60.0 +(59.95571 + paras.Dw1_0)/3600.0)*deg;
     double w11 = (1732559343.73604 + paras.Dw1_1)*sec;
     double w12 = (-6.8084 + paras.Dw1_2)*sec;
     double w13 = (0.006604 + paras.Dw1_3)*sec;
     double w14 = (-3.169e-5 + paras.Dw1_4)*sec;
     fin << "let W1=" << setprecision(17) << w10;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w11 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w12 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w13 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w14 << "*T4),";

     double w20 = (83.0 + 21.0/60.0 + (11.67475 + paras.Dw2_0)/3600.0)*deg;
     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*sec;
     double w22 = (-38.2631 + paras.Dw2_2)*sec;
     double w23 = (-0.045047+ paras.Dw2_3)*sec;
     double w24 = 0.00021301*sec;
     fin << "W2=" << setprecision(17) << w20;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w21 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w22 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w23 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w24 << "*T4),";

     double w30 = (125.0 + 2.0/60.0 + (40.39816 + paras.Dw3_0)/3600.0)*deg;
     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*sec;
     double w32 = (6.359 + paras.Dw3_2)*sec;
     double w33 = (0.007625 + paras.Dw3_3)*sec;
     double w34 = -3.586e-5*sec;
     fin << "W3=" << setprecision(17) << w30;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w31 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w32 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w33 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w34 << "*T4),";

     double Ea0 = (100.0 + 27.0/60.0 + (59.13885 + paras.Deart_0)/3600.0)*deg;
     double Ea1 = (129597742.293 + paras.Deart_1)*sec;
     double Ea2 = -0.0202*sec;
     double Ea3 = 9e-6*sec;
     double Ea4 = 1.5e-7*sec;
     fin << "Ea=" << setprecision(17) << Ea0;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ea1 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ea2 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ea3 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ea4 << "*T4),";

     double p0 = (102.0 + 56.0/60.0 + (14.45766 + paras.Dperi)/3600.0)*deg;
     double p1 = 1161.24342*sec;
     double p2 = 0.529265*sec;
     double p3 = -1.1814e-4*sec;
     double p4 = 1.1379e-5*sec;
     fin << "pp=" << setprecision(17) << p0;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << p1 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << p2 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << p3 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << p4 << "*T4),";

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
     fin << "Me=" << setprecision(17) << Me;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Me1 << "*T),";
     fin << "Ve=" << setprecision(17) << Ve;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ve1 << "*T),";
     fin << "EM=" << setprecision(17) << EM;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << EM1 << "*T),";
     fin << "Ma=" << setprecision(17) << Ma; 
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ma1 << "*T),";
     fin << "Ju=" << setprecision(17) << Ju; 
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ju1 << "*T),";
     fin << "Sa=" << setprecision(17) << Sa;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Sa1 << "*T),";
     fin << "Ur=" << setprecision(17) << Ur;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ur1 << "*T),";
     fin << "Ne=" << setprecision(17) << Ne;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ne1 << "*T);";
     fin << "return {";
     fin << "W:mod2pi" << funSuffix << "(W1),";
     fin << "a:mod2pi" << funSuffix << "(W1-Ea+Math.PI),";
     fin << "b:mod2pi" << funSuffix << "(W1-W3),";
     fin << "c:mod2pi" << funSuffix << "(W1-W2),";
     fin << "d:mod2pi" << funSuffix << "(Ea-pp),";
     fin << "m:mod2pi" << funSuffix << "(W1+0.02438029560881907*T),";
     fin << "e:mod2pi" << funSuffix << "(Me),";
     fin << "f:mod2pi" << funSuffix << "(Ve),";
     fin << "g:mod2pi" << funSuffix << "(EM),";
     fin << "h:mod2pi" << funSuffix << "(Ma),";
     fin << "i:mod2pi" << funSuffix << "(Ju),";
     fin << "j:mod2pi" << funSuffix << "(Sa),";
     fin << "k:mod2pi" << funSuffix << "(Ur),";
     fin << "l:mod2pi" << funSuffix << "(Ne)};}";
}

void generate_javascript_code_min_compute_Elp_arguments_dot(ostream &fin, const char *funSuffix, 
              Elp_paras &paras) {
     fin << "function compute_Elp_arguments_dot" << funSuffix << "(T){";
     fin << "let T2=T*T,T3=T*T2,T4=T2*T2;";
     const double fac = PI/648000.0/36525; // arcsecs -> radians/cy

     double w11 = (1732559343.73604 + paras.Dw1_1)*fac;
     double w12 = (-6.8084 + paras.Dw1_2)*2*fac;
     double w13 = (0.006604 + paras.Dw1_3)*3*fac;
     double w14 = (-3.169e-5 + paras.Dw1_4)*4*fac;
     fin << "let W1=" << setprecision(17) << w11 << w12 << "*T+"
         << w13 << "*T2" << w14 << "*T3,";

     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*fac;
     double w22 = (-38.2631 + paras.Dw2_2)*2*fac;
     double w23 = (-0.045047+ paras.Dw2_3)*3*fac;
     double w24 = 0.00021301*4*fac;
     fin << "W2=" << setprecision(17) << w21 << w22 << "*T" << w23 << "*T2+"
         << w24 << "*T3,";

     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*fac;
     double w32 = (6.359 + paras.Dw3_2)*2*fac;
     double w33 = (0.007625 + paras.Dw3_3)*3*fac;
     double w34 = -3.586e-5*4*fac;
     fin << "W3=" << setprecision(17) << w31 << "+" << w32 << "*T+"
         << w33 << "*T2" << w34 << "*T3,";

     double Ea1 = (129597742.293 + paras.Deart_1)*fac;
     double Ea2 = -0.0202*2*fac;
     double Ea3 = 9e-6*3*fac;
     double Ea4 = 1.5e-7*4*fac;
     fin << "Ea=" << setprecision(17) << Ea1 << Ea2 << "*T+"
         << Ea3 << "*T2+" << Ea4 << "*T3,";

     double p1 = 1161.24342*fac;
     double p2 = 0.529265*2*fac;
     double p3 = -1.1814e-4*3*fac;
     double p4 = 1.1379e-5*4*fac;
     fin << "pp=" << setprecision(17) << p1 << "+" << p2 << "*T"
         << p3 << "*T2+" << p4 << "*T3;";

     fin << "return {W:W1,a:W1-Ea,b:W1-W3,c:W1-W2,d:Ea-pp," 
         << "m:W1+" << setprecision(17) << 5028.79695*fac << ","
         << "e:" << 538101628.66888*fac << ","
         << "f:" << 210664136.45777*fac << ","
         << "g:" << 129597742.293*fac << ","
         << "h:" << 68905077.65936*fac << "," 
         << "i:" << 10925660.57335*fac << "," 
         << "j:" << 4399609.33632*fac << "," 
         << "k:" << 1542482.57845*fac << "," 
         << "l:" << 786547.897*fac << "};}";
}

void generate_javascript_code_min_main_problem(ostream &fin, const char *funSuffix, 
        Elp_coefs &coefs) {
     fin << "function Elp_main_long" << funSuffix << "(a){";
     fin << "let s=0;";
     write_main_problem_function_min(fin, coefs.n_main_long, coefs.i_main_long, 
              coefs.A_main_long, 0);

     fin << "function Elp_main_lat" << funSuffix << "(a){";
     fin << "let s=0;";
     write_main_problem_function_min(fin, coefs.n_main_lat, coefs.i_main_lat,
              coefs.A_main_lat, 0);

     fin << "function Elp_main_dist" << funSuffix << "(a){";
     fin << "let s=" << setprecision(17) << coefs.A_main_dist[0] <<";";
     write_main_problem_function_min(fin, coefs.n_main_dist, coefs.i_main_dist,
              coefs.A_main_dist, 1);
}

void generate_javascript_code_min_main_problem_and_derv(ostream &fin, const char *funSuffix, 
        Elp_coefs &coefs) {
     fin << "function Elp_main_long_and_derv" << funSuffix << "(a,b){";
     fin << "let s=0,v=0,p,q;";
     write_main_problem_and_derv_function_min(fin, coefs.n_main_long, coefs.i_main_long, 
              coefs.A_main_long, 0);

     fin << "function Elp_main_lat_and_derv" << funSuffix << "(a,b){";
     fin << "let s=0,v=0,p,q;";
     write_main_problem_and_derv_function_min(fin, coefs.n_main_lat, coefs.i_main_lat,
              coefs.A_main_lat, 0);

     fin << "function Elp_main_dist_and_derv" << funSuffix << "(a,b){";
     fin << "let s=" << setprecision(17) << coefs.A_main_dist[0] <<",v=0,p,q;";
     write_main_problem_and_derv_function_min(fin, coefs.n_main_dist, coefs.i_main_dist,
              coefs.A_main_dist, 1);
}

void generate_javascript_code_min_perturbation(ostream &fin, const char *funSuffix,
        Elp_coefs &coefs) {
     fin << "function Elp_pert_longT0" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_longT0, coefs.i_pert_longT0, 
           coefs.A_pert_longT0, coefs.ph_pert_longT0);

     fin << "function Elp_pert_longT1" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_longT1, coefs.i_pert_longT1,
           coefs.A_pert_longT1, coefs.ph_pert_longT1);

     fin << "function Elp_pert_longT2" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_longT2, coefs.i_pert_longT2,
           coefs.A_pert_longT2, coefs.ph_pert_longT2);

     fin << "function Elp_pert_longT3" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_longT3, coefs.i_pert_longT3,
           coefs.A_pert_longT3, coefs.ph_pert_longT3);

     fin << "function Elp_pert_latT0" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_latT0, coefs.i_pert_latT0,
           coefs.A_pert_latT0, coefs.ph_pert_latT0);

     fin << "function Elp_pert_latT1" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_latT1, coefs.i_pert_latT1,
           coefs.A_pert_latT1, coefs.ph_pert_latT1);

     fin << "function Elp_pert_latT2" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_latT2, coefs.i_pert_latT2,
           coefs.A_pert_latT2, coefs.ph_pert_latT2);

     fin << "function Elp_pert_distT0" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_distT0, coefs.i_pert_distT0,
           coefs.A_pert_distT0, coefs.ph_pert_distT0);

     fin << "function Elp_pert_distT1" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_distT1, coefs.i_pert_distT1,
           coefs.A_pert_distT1, coefs.ph_pert_distT1);

     fin << "function Elp_pert_distT2" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_distT2, coefs.i_pert_distT2,
           coefs.A_pert_distT2, coefs.ph_pert_distT2);

     fin << "function Elp_pert_distT3" << funSuffix << "(a){";
     write_perturbation_function_min(fin, coefs.n_pert_distT3, coefs.i_pert_distT3,
           coefs.A_pert_distT3, coefs.ph_pert_distT3);
}

void generate_javascript_code_min_perturbation_and_derv(ostream &fin, const char *funSuffix,
        Elp_coefs &coefs) {
     fin << "function Elp_pert_longT0_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_longT0, coefs.i_pert_longT0, 
           coefs.A_pert_longT0, coefs.ph_pert_longT0);

     fin << "function Elp_pert_longT1_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_longT1, coefs.i_pert_longT1,
           coefs.A_pert_longT1, coefs.ph_pert_longT1);

     fin << "function Elp_pert_longT2_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_longT2, coefs.i_pert_longT2,
           coefs.A_pert_longT2, coefs.ph_pert_longT2);

     fin << "function Elp_pert_longT3_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_longT3, coefs.i_pert_longT3,
           coefs.A_pert_longT3, coefs.ph_pert_longT3);

     fin << "function Elp_pert_latT0_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_latT0, coefs.i_pert_latT0,
           coefs.A_pert_latT0, coefs.ph_pert_latT0);

     fin << "function Elp_pert_latT1_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_latT1, coefs.i_pert_latT1,
           coefs.A_pert_latT1, coefs.ph_pert_latT1);

     fin << "function Elp_pert_latT2_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_latT2, coefs.i_pert_latT2,
           coefs.A_pert_latT2, coefs.ph_pert_latT2);

     fin << "function Elp_pert_distT0_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_distT0, coefs.i_pert_distT0,
           coefs.A_pert_distT0, coefs.ph_pert_distT0);

     fin << "function Elp_pert_distT1_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_distT1, coefs.i_pert_distT1,
           coefs.A_pert_distT1, coefs.ph_pert_distT1);

     fin << "function Elp_pert_distT2_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_distT2, coefs.i_pert_distT2,
           coefs.A_pert_distT2, coefs.ph_pert_distT2);

     fin << "function Elp_pert_distT3_and_derv" << funSuffix << "(a,b){";
     write_perturbation_and_derv_function_min(fin, coefs.n_pert_distT3, coefs.i_pert_distT3,
           coefs.A_pert_distT3, coefs.ph_pert_distT3);
}

void generate_javascript_code_min(const char* outfile, Elp_paras &paras, Elp_coefs &coefs, 
                                  const char *funSuffix) {
     ofstream fin(outfile);
     fin << "\"use strict\";";
     fin << "let mod2pi" << funSuffix << "=function(x){";
     fin << "return x-6.283185307179586*Math.floor(0.5*(x*0.3183098861837907+1));};";
     generate_javascript_code_min_getX2000(fin, funSuffix);
     generate_javascript_code_min_compute_Elp_arguments(fin, funSuffix, paras);
     generate_javascript_code_min_main_problem(fin, funSuffix, coefs);
     generate_javascript_code_min_perturbation(fin, funSuffix, coefs);
     fin.close();
}

void generate_javascript_code_with_velocity_min(const char* outfile, Elp_paras &paras, Elp_coefs &coefs,
                                  const char *funSuffix) {
     ofstream fin(outfile);
     fin << "\"use strict\";";
     fin << "let mod2pi" << funSuffix << "=function(x){";
     fin << "return x-6.283185307179586*Math.floor(0.5*(x*0.3183098861837907+1));};";
     generate_javascript_code_min_getX2000_Xdot2000(fin, funSuffix);
     generate_javascript_code_min_compute_Elp_arguments(fin, funSuffix, paras);
     generate_javascript_code_min_compute_Elp_arguments_dot(fin, funSuffix, paras);
     generate_javascript_code_min_main_problem_and_derv(fin, funSuffix, coefs);
     generate_javascript_code_min_perturbation_and_derv(fin, funSuffix, coefs);
     fin.close();
}
