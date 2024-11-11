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

// write function for perturbation 
void write_perturbation_function(ostream &fin, int n, int ** &i_pert, double * &A_pert, 
                           double * &ph_pert) {
  fin << "  var s = 0.0;" << endl;

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
         fin << " + ";
       } else {
         fin << " - ";
       }
       fin << fabs(ph_pert[i]);
     }
     fin << ");" << endl;
  }
  fin << "  return s;" << endl;
  fin << "}" << endl << endl;
}

// funSuffix will be added to every JavaScript function. This is useful to 
// prevent conflicts with code that implements other versions of the truncated series.
void generate_javascript_code(const char* outfile, int corr, 
                       double AthU, double AthV, double AthR, double tau, 
                       Elp_paras &paras, Elp_coefs &coefs, const char* funSuffix) {
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
     fin << "//  This code is generated by \"ElpMpp_JavaScript.cpp\" using the following parameters: " << endl;
     fin << "//  corr = " << corr << "," << endl;
     fin << "//  AthU = " << setprecision(17) << AthU << "," << endl;
     fin << "//  AthV = " << setprecision(17) << AthV << "," << endl;
     fin << "//  AthR = " << setprecision(17) << AthR << "," << endl;
     fin << "//  tau = " << setprecision(17) << tau << "," << endl;
     fin << "//" << endl;
     fin << "//  Usage: Simply call the function getX2000" << funSuffix << "() to" << endl;
     fin << "//         compute the rectangular geocentric coordinates of the Moon" << endl;
     fin << "//         with respect to the mean ecliptic and equinox of J2000.0." << endl;
     fin << "//" << endl;
     fin << "// **Note: You should use the minified version of the JavaScript code" << endl;
     fin << "//         instead of this file to optimize performance.**" << endl;
     fin << "//" << endl;
     fin << "// ---------------------------------------------------------------- " << endl;
     fin << endl;
     fin << "\"use strict\";" << endl << endl;
     fin << "// restrict x to [-pi,pi) " << endl;
     fin << "var mod2pi" << funSuffix << " = function(x) {" << endl;
     fin << "  return x - 6.283185307179586*Math.floor(0.5*(x*0.3183098861837907 + 1));" << endl;
     fin << "};" << endl << endl;
     fin << "// Calculate the Moon's geocentric X,Y,Z coordinates with respect to " << endl;
     fin << "// J2000.0 mean ecliptic and equinox." << endl;
     fin << "function getX2000" << funSuffix << "(T) {" << endl;
     fin << "  var T2 = T*T;" << endl;
     fin << "  var T3 = T*T2;" << endl;
     fin << "  var T4 = T2*T2;" << endl;
     fin << "  var T5 = T2*T3;" << endl << endl;
     fin << "  // Moon's longitude, latitude and distance" << endl;
     fin << "  var args = compute_Elp_arguments" << funSuffix << "(T);" << endl;
     fin << "  var longM = args.W1 + Elp_main_long" << funSuffix << "(args) + Elp_pert_longT0" << funSuffix << "(args) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_longT1" << funSuffix << "(args)*T) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_longT2" << funSuffix << "(args)*T2) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_longT3" << funSuffix << "(args)*T3);" << endl;
     fin << "  var latM =  Elp_main_lat" << funSuffix << "(args) + Elp_pert_latT0" << funSuffix << "(args) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_latT1" << funSuffix << "(args)*T) +" << endl;
     fin << "              mod2pi" << funSuffix << "(Elp_pert_latT2" << funSuffix << "(args)*T2);" << endl;
     fin << "  var r = " << setprecision(17) << 384747.961370173/384747.980674318 << "*" 
         << "(Elp_main_dist" << funSuffix << "(args) + Elp_pert_distT0" << funSuffix << "(args) +" << endl;
     fin << "              Elp_pert_distT1" << funSuffix << "(args)*T +" << endl;
     fin << "              Elp_pert_distT2" << funSuffix << "(args)*T2 +" << endl;
     fin << "              Elp_pert_distT3" << funSuffix << "(args)*T3);" << endl << endl;
     fin << "  var x0 = r*Math.cos(longM)*Math.cos(latM);" << endl;
     fin << "  var y0 = r*Math.sin(longM)*Math.cos(latM);" << endl;
     fin << "  var z0 = r*Math.sin(latM);" << endl;
     fin << "" << endl;
     fin << "  // Precession matrix" << endl;
     fin << "  var P = 0.10180391e-4*T + 0.47020439e-6*T2 - 0.5417367e-9*T3 " << endl;
     fin << "             - 0.2507948e-11*T4 + 0.463486e-14*T5;" << endl;
     fin << "  var Q = -0.113469002e-3*T + 0.12372674e-6*T2 + 0.12654170e-8*T3 " << endl;
     fin << "             - 0.1371808e-11*T4 - 0.320334e-14*T5;" << endl;
     fin << "  var sq = Math.sqrt(1 - P*P - Q*Q);" << endl;
     fin << "  var p11 = 1 - 2*P*P;" << endl;
     fin << "  var p12 = 2*P*Q;" << endl;
     fin << "  var p13 = 2*P*sq;" << endl;
     fin << "  var p21 = 2*P*Q;" << endl;
     fin << "  var p22 = 1-2*Q*Q;" << endl;
     fin << "  var p23 = -2*Q*sq;" << endl;
     fin << "  var p31 = -2*P*sq;" << endl;
     fin << "  var p32 = 2*Q*sq;" << endl;
     fin << "  var p33 = 1 - 2*P*P - 2*Q*Q;" << endl << endl;
     fin << "  // Finally, components of position vector wrt J2000.0 mean ecliptic and equinox" << endl;
     fin << "  var X = p11*x0 + p12*y0 + p13*z0;" << endl;
     fin << "  var Y = p21*x0 + p22*y0 + p23*z0;" << endl;
     fin << "  var Z = p31*x0 + p32*y0 + p33*z0;" << endl << endl;
     fin << "  return {X:X, Y:Y, Z:Z, rGeo:r};" << endl;
     fin << "}" << endl << endl;

     fin << "function compute_Elp_arguments" << funSuffix << "(T) {" << endl;
     fin << "  var T2 = T*T;" << endl;
     fin << "  var T3 = T*T2;" << endl;
     fin << "  var T4 = T2*T2;" << endl << endl;
     const double deg = PI/180.0; // degrees -> radians
     const double sec = PI/648000.0; // arcsecs -> radians
     
     double w10 = (-142.0 + 18.0/60.0 +(59.95571 + paras.Dw1_0)/3600.0)*deg;
     double w11 = (1732559343.73604 + paras.Dw1_1)*sec;
     double w12 = (-6.8084 + paras.Dw1_2)*sec;
     double w13 = (0.006604 + paras.Dw1_3)*sec;
     double w14 = (-3.169e-5 + paras.Dw1_4)*sec;
     fin << "  var W1 = " << setprecision(17) << w10;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w11 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w12 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w13 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w14 << "*T4);" << endl;

     double w20 = (83.0 + 21.0/60.0 + (11.67475 + paras.Dw2_0)/3600.0)*deg;
     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*sec;
     double w22 = (-38.2631 + paras.Dw2_2)*sec;
     double w23 = (-0.045047+ paras.Dw2_3)*sec;
     double w24 = 0.00021301*sec;
     fin << "  var W2 = " << setprecision(17) << w20;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w21 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w22 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w23 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w24 << "*T4);" << endl;

     double w30 = (125.0 + 2.0/60.0 + (40.39816 + paras.Dw3_0)/3600.0)*deg;
     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*sec;
     double w32 = (6.359 + paras.Dw3_2)*sec;
     double w33 = (0.007625 + paras.Dw3_3)*sec;
     double w34 = -3.586e-5*sec;
     fin << "  var W3 = " << setprecision(17) << w30;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w31 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w32 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w33 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << w34 << "*T4);" << endl;

     double Ea0 = (100.0 + 27.0/60.0 + (59.13885 + paras.Deart_0)/3600.0)*deg;
     double Ea1 = (129597742.293 + paras.Deart_1)*sec;
     double Ea2 = -0.0202*sec;
     double Ea3 = 9e-6*sec;
     double Ea4 = 1.5e-7*sec;
     fin << "  var Ea = " << setprecision(17) << Ea0;
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << Ea1 << "*T)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << Ea2 << "*T2)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << Ea3 << "*T3)";
     fin << "  + mod2pi" << funSuffix << "(" << setprecision(17) << Ea4 << "*T4);" << endl;

     double p0 = (102.0 + 56.0/60.0 + (14.45766 + paras.Dperi)/3600.0)*deg;
     double p1 = 1161.24342*sec;
     double p2 = 0.529265*sec;
     double p3 = -1.1814e-4*sec;
     double p4 = 1.1379e-5*sec;
     fin << "  var pomp = " << setprecision(17) << p0;
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
     fin << "  var Me = " << setprecision(17) << Me;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Me1 << "*T);" << endl;
     fin << "  var Ve = " << setprecision(17) << Ve;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ve1 << "*T);" << endl;
     fin << "  var EM = " << setprecision(17) << EM;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << EM1 << "*T);" << endl;
     fin << "  var Ma = " << setprecision(17) << Ma;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ma1 << "*T);" << endl;
     fin << "  var Ju = " << setprecision(17) << Ju;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ju1 << "*T);" << endl;
     fin << "  var Sa = " << setprecision(17) << Sa;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Sa1 << "*T);" << endl;
     fin << "  var Ur = " << setprecision(17) << Ur;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ur1 << "*T);" << endl;
     fin << "  var Ne = " << setprecision(17) << Ne;
     fin << " + mod2pi" << funSuffix << "(" << setprecision(17) << Ne1 << "*T);" << endl << endl;
     fin << "  var args = {};" << endl;
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

     fin << "// Sum the ELP/MPP02 series: main problem, longitude (" 
         << coefs.n_main_long;
     if (coefs.n_main_long > 1) {
       fin << " terms)" << endl;
     } else {
       fin << " term)" << endl;
     }
     fin << "function Elp_main_long" << funSuffix << "(args) {" << endl;
     fin << "  var s = 0.0;" << endl;
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
     fin << "  var s = 0.0;" << endl;
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
     fin << "  var s = " << setprecision(17) << coefs.A_main_dist[0] <<";" << endl;
     write_main_problem_function(fin, coefs.n_main_dist, coefs.i_main_dist,
              coefs.A_main_dist, 1);

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

// write function for perturbation 
void write_perturbation_function_min(ostream &fin, int n, int ** &i_pert, double * &A_pert, 
                           double * &ph_pert) {
  fin << "var s=0;";

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

void generate_javascript_code_min(const char* outfile, int corr, 
                       double AthU, double AthV, double AthR, double tau, 
                       Elp_paras &paras, Elp_coefs &coefs, const char *funSuffix) {
     ofstream fin(outfile);
     fin << "\"use strict\";";
     fin << "var mod2pi" << funSuffix << "=function(x){";
     fin << "return x-6.283185307179586*Math.floor(0.5*(x*0.3183098861837907+1));};";
     fin << "function getX2000" << funSuffix << "(T){";
     fin << "var T2=T*T;";
     fin << "var T3=T*T2;";
     fin << "var T4=T2*T2;";
     fin << "var T5=T2*T3;";
     fin << "var a=compute_Elp_arguments" << funSuffix << "(T);";
     fin << "var U=a.W+Elp_main_long" << funSuffix << "(a)+Elp_pert_longT0" << funSuffix << "(a)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_longT1" << funSuffix << "(a)*T)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_longT2" << funSuffix << "(a)*T2)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_longT3" << funSuffix << "(a)*T3);";
     fin << "var V=Elp_main_lat" << funSuffix << "(a)+Elp_pert_latT0" << funSuffix << "(a)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_latT1" << funSuffix << "(a)*T)+";
     fin << "mod2pi" << funSuffix << "(Elp_pert_latT2" << funSuffix << "(a)*T2);";
     fin << "var r=" << setprecision(17) << 384747.961370173/384747.980674318 << "*" 
         << "(Elp_main_dist" << funSuffix << "(a)+Elp_pert_distT0" << funSuffix << "(a)+";
     fin << "Elp_pert_distT1" << funSuffix << "(a)*T+";
     fin << "Elp_pert_distT2" << funSuffix << "(a)*T2+";
     fin << "Elp_pert_distT3" << funSuffix << "(a)*T3);";
     fin << "var x0=r*Math.cos(U)*Math.cos(V);";
     fin << "var y0=r*Math.sin(U)*Math.cos(V);";
     fin << "var z0=r*Math.sin(V);";
     fin << "var P=0.10180391e-4*T+0.47020439e-6*T2-0.5417367e-9*T3";
     fin << "-0.2507948e-11*T4+0.463486e-14*T5;";
     fin << "var Q=-0.113469002e-3*T+0.12372674e-6*T2+0.12654170e-8*T3";
     fin << "-0.1371808e-11*T4-0.320334e-14*T5;";
     fin << "var sq=Math.sqrt(1-P*P-Q*Q);";
     fin << "var p11=1-2*P*P;";
     fin << "var p12=2*P*Q;";
     fin << "var p13=2*P*sq;";
     fin << "var p21=2*P*Q;";
     fin << "var p22=1-2*Q*Q;";
     fin << "var p23=-2*Q*sq;";
     fin << "var p31=-2*P*sq;";
     fin << "var p32=2*Q*sq;";
     fin << "var p33=1-2*P*P-2*Q*Q;";
     fin << "var X=p11*x0+p12*y0+p13*z0;";
     fin << "var Y=p21*x0+p22*y0+p23*z0;";
     fin << "var Z=p31*x0+p32*y0+p33*z0;";
     fin << "return {X:X,Y:Y,Z:Z,rGeo:r};}";

     fin << "function compute_Elp_arguments" << funSuffix << "(T){";
     fin << "var T2=T*T;";
     fin << "var T3=T*T2;";
     fin << "var T4=T2*T2;";
     const double deg = PI/180.0; // degrees -> radians
     const double sec = PI/648000.0; // arcsecs -> radians
     
     double w10 = (-142.0 + 18.0/60.0 +(59.95571 + paras.Dw1_0)/3600.0)*deg;
     double w11 = (1732559343.73604 + paras.Dw1_1)*sec;
     double w12 = (-6.8084 + paras.Dw1_2)*sec;
     double w13 = (0.006604 + paras.Dw1_3)*sec;
     double w14 = (-3.169e-5 + paras.Dw1_4)*sec;
     fin << "var W1=" << setprecision(17) << w10;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w11 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w12 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w13 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w14 << "*T4);";

     double w20 = (83.0 + 21.0/60.0 + (11.67475 + paras.Dw2_0)/3600.0)*deg;
     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*sec;
     double w22 = (-38.2631 + paras.Dw2_2)*sec;
     double w23 = (-0.045047+ paras.Dw2_3)*sec;
     double w24 = 0.00021301*sec;
     fin << "var W2=" << setprecision(17) << w20;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w21 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w22 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w23 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w24 << "*T4);";

     double w30 = (125.0 + 2.0/60.0 + (40.39816 + paras.Dw3_0)/3600.0)*deg;
     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*sec;
     double w32 = (6.359 + paras.Dw3_2)*sec;
     double w33 = (0.007625 + paras.Dw3_3)*sec;
     double w34 = -3.586e-5*sec;
     fin << "var W3=" << setprecision(17) << w30;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w31 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w32 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w33 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << w34 << "*T4);";

     double Ea0 = (100.0 + 27.0/60.0 + (59.13885 + paras.Deart_0)/3600.0)*deg;
     double Ea1 = (129597742.293 + paras.Deart_1)*sec;
     double Ea2 = -0.0202*sec;
     double Ea3 = 9e-6*sec;
     double Ea4 = 1.5e-7*sec;
     fin << "var Ea=" << setprecision(17) << Ea0;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ea1 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ea2 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ea3 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ea4 << "*T4);";

     double p0 = (102.0 + 56.0/60.0 + (14.45766 + paras.Dperi)/3600.0)*deg;
     double p1 = 1161.24342*sec;
     double p2 = 0.529265*sec;
     double p3 = -1.1814e-4*sec;
     double p4 = 1.1379e-5*sec;
     fin << "var pp = " << setprecision(17) << p0;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << p1 << "*T)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << p2 << "*T2)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << p3 << "*T3)";
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << p4 << "*T4);";

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
     fin << "var Me=" << setprecision(17) << Me;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Me1 << "*T);";
     fin << "var Ve=" << setprecision(17) << Ve;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ve1 << "*T);";
     fin << "var EM=" << setprecision(17) << EM;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << EM1 << "*T);";
     fin << "var Ma=" << setprecision(17) << Ma; 
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ma1 << "*T);";
     fin << "var Ju=" << setprecision(17) << Ju; 
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ju1 << "*T);";
     fin << "var Sa=" << setprecision(17) << Sa;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Sa1 << "*T);";
     fin << "var Ur=" << setprecision(17) << Ur;
     fin << "+mod2pi" << funSuffix << "(" << setprecision(17) << Ur1 << "*T);";
     fin << "var Ne=" << setprecision(17) << Ne;
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

     fin << "function Elp_main_long" << funSuffix << "(a){";
     fin << "var s=0;";
     write_main_problem_function_min(fin, coefs.n_main_long, coefs.i_main_long, 
              coefs.A_main_long, 0);

     fin << "function Elp_main_lat" << funSuffix << "(a){";
     fin << "var s=0;";
     write_main_problem_function_min(fin, coefs.n_main_lat, coefs.i_main_lat,
              coefs.A_main_lat, 0);

     fin << "function Elp_main_dist" << funSuffix << "(a){";
     fin << "var s=" << setprecision(17) << coefs.A_main_dist[0] <<";";
     write_main_problem_function_min(fin, coefs.n_main_dist, coefs.i_main_dist,
              coefs.A_main_dist, 1);

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

     fin.close();
}
