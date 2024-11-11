// ----------------------------------------------------------------
//  This code calculates the ELP/MPP02 series.
//
//  ELP/MPP02 is a semi-analytic solution for the lunar motion developed by
//  J. Chapront and G. Francou in 2002. It is an improvement of the ELP2000-82B
//  lunar theory.
//
//  Source paper:
//    The lunar theory ELP revisited. Introduction of new planetary perturbations
//    by J. Chapront and G. Francou, Astronomy and Astrophysics, v.404, p.735-742 (2003)
//    http://adsabs.harvard.edu/abs/2003A%26A...404..735C
//
//  This code and data files are based on the authors' FORTRAN code and data files on
//  ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/
//
//  The following 14 data files are required:
//    elp_main.long, elp_main.lat, elp_main.dist,
//    elp_pert.longT0, elp_pert.longT1, elp_pert.longT2, elp_pert.longT3,
//    elp_pert.latT0, elp_pert.latT1, elp_pert.latT2,
//    elp_pert.distT0, elp_pert.distT1, elp_pert.distT2, elp_pert.distT3
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
//       the Elp/MPP02 series for the main problem, stored in the struct 
//       Elp_facs.
//    3. Call the function setup_Elp_coefs() to set up the coefficients 
//       for the ELP/MPP02 series. The coefficients are stored in the struct 
//       Elp_coefs.
//    4. Call getX2000() to compute the rectangular geocentric coordinates 
//       of the Moon's position with respect to the mean ecliptic and 
//       equinox of J2000.0 or getX2000_Xdot2000() if time derivatives 
//       are required.
// 
//  See example.cpp for an example of using this code.
//---------------------------------------------------------------- 
#include "ElpMpp02.h"

// restrict x to [-pi, pi)
double mod2pi(double x) {
  const double tpi = 2.0*PI;
  return x - tpi*floor((x + PI)/tpi);
}

// Set up adjustable parameters 
// corr=0: fit to LLR data, corr=1: fit to DE405
void setup_parameters(int corr, Elp_paras &paras, Elp_facs &facs) {
   // PARAMETERS adjusted to fit data
   switch (corr) {
       case 0: 
         paras.Dw1_0   = -0.10525;
         paras.Dw2_0   =  0.16826;
         paras.Dw3_0   = -0.10760;
         paras.Deart_0 = -0.04012;
         paras.Dperi   = -0.04854;
         paras.Dw1_1   = -0.32311;
         paras.Dgam    =  0.00069;
         paras.De      =  0.00005;
         paras.Deart_1 =  0.01442;
         paras.Dep     =  0.00226;
         paras.Dw2_1   =  0.08017;
         paras.Dw3_1   = -0.04317;
         paras.Dw1_2   = -0.03794;
         paras.Dw1_3   =  0.0;
         paras.Dw1_4   =  0.0;
         paras.Dw2_2   =  0.0;
         paras.Dw2_3   =  0.0;
         paras.Dw3_2   =  0.0;
         paras.Dw3_3   =  0.0;
         break;
      case 1: 
         paras.Dw1_0   = -0.07008;
         paras.Dw2_0   =  0.20794;
         paras.Dw3_0   = -0.07215;
         paras.Deart_0 = -0.00033;
         paras.Dperi   = -0.00749;
         paras.Dw1_1   = -0.35106;
         paras.Dgam    =  0.00085;
         paras.De      = -0.00006;
         paras.Deart_1 =  0.00732;
         paras.Dep     =  0.00224;
         paras.Dw2_1   =  0.08017;
         paras.Dw3_1   = -0.04317;
         paras.Dw1_2   = -0.03743;
         paras.Dw1_3   = -0.00018865;
         paras.Dw1_4   = -0.00001024;
         paras.Dw2_2   =  0.00470602;
         paras.Dw2_3   = -0.00025213;
         paras.Dw3_2   = -0.00261070;
         paras.Dw3_3   = -0.00010712;
         break;
   }
   // derived parameters
   const double am = 0.074801329;
   const double alpha = 0.002571881;
   const double dtsm = 2.0*alpha/(3.0*am);
   const double xa = 2.0*alpha/3.0;
   const double sec = PI/648000.0; // arcsecs -> radians
   double bp[5][2] = {{0.311079095, -0.103837907},
                      {-0.004482398, 0.000668287},
                      {-0.001102485, -0.001298072},
                      {0.001056062, -0.000178028},
                      {0.000050928, -0.000037342}};
   double w11 = (1732559343.73604 + paras.Dw1_1)*sec;
   double w21 = (14643420.3171 + paras.Dw2_1)*sec;
   double w31 = (-6967919.5383 + paras.Dw3_1)*sec;
   double x2 = w21/w11;
   double x3 = w31/w11;
   double y2 = am*bp[0][0] + xa*bp[4][0];
   double y3 = am*bp[0][1] + xa*bp[4][1];
   double d21 = x2-y2;
   double d22 = w11*bp[1][0];
   double d23 = w11*bp[2][0];
   double d24 = w11*bp[3][0];
   double d25 = y2/am;
   double d31 = x3-y3;
   double d32 = w11*bp[1][1];
   double d33 = w11*bp[2][1];
   double d34 = w11*bp[3][1];
   double d35 = y3/am;
   paras.Cw2_1 = d21*paras.Dw1_1 + d25*paras.Deart_1 + d22*paras.Dgam +
                  d23*paras.De + d24*paras.Dep;
   paras.Cw3_1 = d31*paras.Dw1_1 + d35*paras.Deart_1 + d32*paras.Dgam +
                  d33*paras.De + d34*paras.Dep;

   // factors multipled by B1-B5 for longitude and latitude
   double delnu_nu = (0.55604 + paras.Dw1_1)*sec/w11;
   double dele = (0.01789 + paras.De)*sec;
   double delg = (-0.08066 + paras.Dgam)*sec;
   double delnp_nu = (-0.06424 + paras.Deart_1)*sec/w11;
   double delep = (-0.12879 + paras.Dep)*sec;
   // factors multipled by B1-B5 for longitude and latitude
   facs.fB1 = -am*delnu_nu + delnp_nu;
   facs.fB2 = delg;
   facs.fB3 = dele;
   facs.fB4 = delep;
   facs.fB5 = -xa*delnu_nu + dtsm*delnp_nu;
   // factor multiplie A_i for distance
   facs.fA = 1.0 - 2.0/3.0*delnu_nu;
}

// Read main problem file
// n is the number of terms in the series, which is stored in the 
//   first line of the data file
void read_main_problem_file(const char *infile, int &n, int ** &i_main, double * &A_main, 
                            double fA, Elp_facs facs) {
  double A,B1,B2,B3,B4,B5,B6;
  int i;
  ifstream file(infile, ios::in);
  if (!file) {
    cerr << "Error in opening " << infile << endl;
    exit(1);
  }
  file >> n;
  i_main = new int *[n];
  A_main = new double[n];
  for (i=0; i<n; i++) i_main[i] = new int[4];
  for (i=0; i<n; i++) {
     if (file.eof()) {
       cerr << "Reached the end of the file " << infile 
            << " before reading all data!" << endl;
       exit(1);
     }
     file >> i_main[i][0] >> i_main[i][1] >> i_main[i][2] >> i_main[i][3]
          >> A >> B1 >> B2 >> B3 >> B4 >> B5 >> B6;
     A_main[i] = fA*A + facs.fB1*B1 + facs.fB2*B2 + facs.fB3*B3 +
                                facs.fB4*B4 + facs.fB5*B5;
  }
  file.close();
}

// Read perturbation file
// n is the number of terms in the series, which is stored in the
//   first line of the data file
void read_perturbation_file(const char *infile, int &n, int ** &i_pert, double * &A_pert,
                            double * &phase) {
  int i;
  ifstream file(infile);
  if (!file) {
    cerr << "Error in opening " << infile << endl;
    exit(1);
  }
  file >> n;
  i_pert = new int *[n];
  A_pert = new double[n];
  phase = new double[n];
  for (i=0; i<n; i++) i_pert[i] = new int[13];
  for (i=0; i<n; i++) {
     if (file.eof()) {
       cerr << "Reached the end of the file " << infile
            << " before reading all data!" << endl;
       exit(1);
     }
     file >> i_pert[i][0] >> i_pert[i][1] >> i_pert[i][2] >> i_pert[i][3]
          >> i_pert[i][4] >> i_pert[i][5] >> i_pert[i][6] >> i_pert[i][7]
          >> i_pert[i][8] >> i_pert[i][9] >> i_pert[i][10] >> i_pert[i][11]
          >> i_pert[i][12] >> A_pert[i] >> phase[i];
  }
  file.close();
}

// set up coefficients for the ELP/MPP02 series
void setup_Elp_coefs(Elp_coefs &coefs, Elp_facs facs) {
  string infile;

  // Main problem
  infile = "elp_main.long";
  read_main_problem_file(infile.c_str(), coefs.n_main_long, coefs.i_main_long,
                         coefs.A_main_long, 1.0, facs);
  infile = "elp_main.lat";
  read_main_problem_file(infile.c_str(), coefs.n_main_lat, coefs.i_main_lat,
                         coefs.A_main_lat, 1.0, facs);
  infile = "elp_main.dist";
  read_main_problem_file(infile.c_str(), coefs.n_main_dist, coefs.i_main_dist, 
                         coefs.A_main_dist, facs.fA, facs);

  // perturbation, longitude
  infile = "elp_pert.longT0";
  read_perturbation_file(infile.c_str(), coefs.n_pert_longT0, coefs.i_pert_longT0, 
                         coefs.A_pert_longT0, coefs.ph_pert_longT0);
  infile = "elp_pert.longT1";
  read_perturbation_file(infile.c_str(), coefs.n_pert_longT1, coefs.i_pert_longT1,
                         coefs.A_pert_longT1, coefs.ph_pert_longT1);
  infile = "elp_pert.longT2";
  read_perturbation_file(infile.c_str(), coefs.n_pert_longT2, coefs.i_pert_longT2,
                         coefs.A_pert_longT2, coefs.ph_pert_longT2);
  infile = "elp_pert.longT3";
  read_perturbation_file(infile.c_str(), coefs.n_pert_longT3, coefs.i_pert_longT3,
                         coefs.A_pert_longT3, coefs.ph_pert_longT3);

  // perturbation, latitude
  infile = "elp_pert.latT0";
  read_perturbation_file(infile.c_str(), coefs.n_pert_latT0, coefs.i_pert_latT0,
                         coefs.A_pert_latT0, coefs.ph_pert_latT0);
  infile = "elp_pert.latT1";
  read_perturbation_file(infile.c_str(), coefs.n_pert_latT1, coefs.i_pert_latT1,
                         coefs.A_pert_latT1, coefs.ph_pert_latT1);
  infile = "elp_pert.latT2";
  read_perturbation_file(infile.c_str(), coefs.n_pert_latT2, coefs.i_pert_latT2,
                         coefs.A_pert_latT2, coefs.ph_pert_latT2);

  // perturbation, distance
  infile = "elp_pert.distT0";
  read_perturbation_file(infile.c_str(), coefs.n_pert_distT0, coefs.i_pert_distT0,
                         coefs.A_pert_distT0, coefs.ph_pert_distT0);
  infile = "elp_pert.distT1";
  read_perturbation_file(infile.c_str(), coefs.n_pert_distT1, coefs.i_pert_distT1,
                         coefs.A_pert_distT1, coefs.ph_pert_distT1);
  infile = "elp_pert.distT2";
  read_perturbation_file(infile.c_str(), coefs.n_pert_distT2, coefs.i_pert_distT2,
                         coefs.A_pert_distT2, coefs.ph_pert_distT2);
  infile = "elp_pert.distT3";
  read_perturbation_file(infile.c_str(), coefs.n_pert_distT3, coefs.i_pert_distT3,
                         coefs.A_pert_distT3, coefs.ph_pert_distT3);
}

// Compute the lunar and planetary arguments used in the ELP/MPP02 series
void compute_Elp_arguments(double T, Elp_paras paras, Elp_args &args) {
     const double deg = PI/180.0; // degrees -> radians
     const double sec = PI/648000.0; // arcsecs -> radians
     double T2 = T*T;
     double T3 = T*T2;
     double T4 = T2*T2;
     double w10 = (-142.0 + 18.0/60.0 +(59.95571 + paras.Dw1_0)/3600.0)*deg;
     double w11 = mod2pi((1732559343.73604 + paras.Dw1_1)*T*sec);
     double w12 = mod2pi((-6.8084 + paras.Dw1_2)*T2*sec);
     double w13 = mod2pi((0.006604 + paras.Dw1_3)*T3*sec);
     double w14 = mod2pi((-3.169e-5 + paras.Dw1_4)*T4*sec);
     double w20 = (83.0 + 21.0/60.0 + (11.67475 + paras.Dw2_0)/3600.0)*deg;
     double w21 = mod2pi((14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*T*sec);
     double w22 = mod2pi((-38.2631 + paras.Dw2_2)*T2*sec);
     double w23 = mod2pi((-0.045047+ paras.Dw2_3)*T3*sec);
     double w24 = mod2pi(0.00021301*T4*sec);
     double w30 = (125.0 + 2.0/60.0 + (40.39816 + paras.Dw3_0)/3600.0)*deg;
     double w31 = mod2pi((-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*T*sec);
     double w32 = mod2pi((6.359 + paras.Dw3_2)*T2*sec);
     double w33 = mod2pi((0.007625 + paras.Dw3_3)*T3*sec);
     double w34 = mod2pi(-3.586e-5*T4*sec);
     double Ea0 = (100.0 + 27.0/60.0 + (59.13885 + paras.Deart_0)/3600.0)*deg;
     double Ea1 = mod2pi((129597742.293 + paras.Deart_1)*T*sec);
     double Ea2 = mod2pi(-0.0202*T2*sec);
     double Ea3 = mod2pi(9e-6*T3*sec);
     double Ea4 = mod2pi(1.5e-7*T4*sec);
     double p0 = (102.0 + 56.0/60.0 + (14.45766 + paras.Dperi)/3600.0)*deg;
     double p1 = mod2pi(1161.24342*T*sec);
     double p2 = mod2pi(0.529265*T2*sec);
     double p3 = mod2pi(-1.1814e-4*T3*sec);
     double p4 = mod2pi(1.1379e-5*T4*sec);
     
     double Me = (-108.0 + 15.0/60.0 + 3.216919/3600.0)*deg;
     Me += mod2pi(538101628.66888*T*sec);
     double Ve = (-179.0 + 58.0/60.0 + 44.758419/3600.0)*deg;
     Ve += mod2pi(210664136.45777*T*sec);
     double EM = (100.0 + 27.0/60.0 + 59.13885/3600.0)*deg;
     EM += mod2pi(129597742.293*T*sec);
     double Ma = (-5.0 + 26.0/60.0 + 3.642778/3600.0)*deg;
     Ma += mod2pi(68905077.65936*T*sec);
     double Ju = (34.0 + 21.0/60.0 + 5.379392/3600.0)*deg;
     Ju += mod2pi(10925660.57335*T*sec);
     double Sa = (50.0 + 4.0/60.0 + 38.902495/3600.0)*deg;
     Sa += mod2pi(4399609.33632*T*sec);
     double Ur = (-46.0 + 3.0/60.0 + 4.354234/3600.0)*deg;
     Ur += mod2pi(1542482.57845*T*sec);
     double Ne = (-56.0 + 20.0/60.0 + 56.808371/3600.0)*deg;
     Ne += mod2pi(786547.897*T*sec);
     
     double W1 = w10+w11+w12+w13+w14; 
     double W2 = w20+w21+w22+w23+w24;
     double W3 = w30+w31+w32+w33+w34;
     double Ea = Ea0+Ea1+Ea2+Ea3+Ea4;
     double pomp = p0+p1+p2+p3+p4;
     
     // Mean longitude of the Moon 
     args.W1 = mod2pi(W1);
     // Arguments of Delaunay
     args.D = mod2pi(W1-Ea + PI);
     args.F = mod2pi(W1-W3);
     args.L = mod2pi(W1-W2);
     args.Lp = mod2pi(Ea-pomp);
     
     // zeta
     args.zeta = mod2pi(W1 + 0.02438029560881907*T);
     
     // Planetary arguments (mean longitudes and mean motions)
     args.Me = mod2pi(Me);
     args.Ve = mod2pi(Ve); 
     args.EM = mod2pi(EM);
     args.Ma = mod2pi(Ma);
     args.Ju = mod2pi(Ju);
     args.Sa = mod2pi(Sa);
     args.Ur = mod2pi(Ur);
     args.Ne = mod2pi(Ne);
}

// Compute the time derivatives of the lunar and planetary arguments used in the ELP/MPP02 series
// return the time derivatives in units of rad/day
void compute_Elp_arguments_dot(double T, Elp_paras paras, Elp_args &args_dot) {
     const double fac = PI/648000.0/36525; // arcsecs -> radians/cy
     double T2 = T*T;
     double T3 = T*T2;
     double w11 = (1732559343.73604 + paras.Dw1_1)*fac;
     double w12 = (-6.8084 + paras.Dw1_2)*2*T*fac;
     double w13 = (0.006604 + paras.Dw1_3)*3*T2*fac;
     double w14 = (-3.169e-5 + paras.Dw1_4)*4*T3*fac;
     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*fac;
     double w22 = (-38.2631 + paras.Dw2_2)*2*T*fac;
     double w23 = (-0.045047+ paras.Dw2_3)*3*T2*fac;
     double w24 = 0.00021301*4*T3*fac;
     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*fac;
     double w32 = (6.359 + paras.Dw3_2)*2*T*fac;
     double w33 = (0.007625 + paras.Dw3_3)*3*T2*fac;
     double w34 = -3.586e-5*4*T3*fac;
     double Ea1 = (129597742.293 + paras.Deart_1)*fac;
     double Ea2 = -0.0202*2*T*fac;
     double Ea3 = 9e-6*3*T2*fac;
     double Ea4 = 1.5e-7*4*T3*fac;
     double p1 = 1161.24342*fac;
     double p2 = 0.529265*2*T*fac;
     double p3 = -1.1814e-4*3*T2*fac;
     double p4 = 1.1379e-5*4*T3*fac;
     
     double W1 = w11+w12+w13+w14; 
     double W2 = w21+w22+w23+w24;
     double W3 = w31+w32+w33+w34;
     double Ea = Ea1+Ea2+Ea3+Ea4;
     double pomp = p1+p2+p3+p4;
     
     // Mean longitude of the Moon 
     args_dot.W1 = W1;
     // Arguments of Delaunay
     args_dot.D = W1-Ea;
     args_dot.F = W1-W3;
     args_dot.L = W1-W2;
     args_dot.Lp = Ea-pomp;
     
     // zeta
     args_dot.zeta = W1 + 5028.79695*fac;
     
     // Planetary arguments (mean longitudes and mean motions)
     args_dot.Me = 538101628.66888*fac;
     args_dot.Ve = 210664136.45777*fac;
     args_dot.EM = 129597742.293*fac;
     args_dot.Ma = 68905077.65936*fac;
     args_dot.Ju = 10925660.57335*fac;
     args_dot.Sa = 4399609.33632*fac;
     args_dot.Ur = 1542482.57845*fac;
     args_dot.Ne = 786547.897*fac;
}

// Sum the ELP/MPP02 series for the main problem
// dist = 0: sine series; dist != 0: cosine series
double Elp_main_sum(int n, int ** &i_main, double * &A_main, Elp_args &args, int dist) {
    int i;
    double sum = 0.0;
    double phase;
    if (dist==0) {
       // sine series 
       for (i=0; i<n; i++) {
          phase = i_main[i][0]*args.D + i_main[i][1]*args.F + i_main[i][2]*args.L + 
                  i_main[i][3]*args.Lp;
          sum += A_main[i]*sin(phase);
       }
    } else {
       // cosine series
       for (i=0; i<n; i++) {
          phase = i_main[i][0]*args.D + i_main[i][1]*args.F + i_main[i][2]*args.L +
                  i_main[i][3]*args.Lp;
          sum += A_main[i]*cos(phase);
       }
    }
    return sum;
}

// Sum the ELP/MPP02 series for perturbations 
double Elp_perturbation_sum(int n, int ** &i_pert, double * &A_pert, double * &ph_pert,
                            Elp_args &args) {
    int i;
    double sum = 0.0;
    double phase;
    for (i=0; i<n; i++) {
       phase = ph_pert[i] + i_pert[i][0]*args.D + i_pert[i][1]*args.F + 
               i_pert[i][2]*args.L + i_pert[i][3]*args.Lp + i_pert[i][4]*args.Me + 
               i_pert[i][5]*args.Ve + i_pert[i][6]*args.EM + i_pert[i][7]*args.Ma + 
               i_pert[i][8]*args.Ju + i_pert[i][9]*args.Sa + i_pert[i][10]*args.Ur +
               i_pert[i][11]*args.Ne + i_pert[i][12]*args.zeta;
       sum += A_pert[i]*sin(phase);
    }
    return sum;
}

// Sum the ELP/MPP02 series and its time derivative for the main problem
// dist = 0: sine series; dist != 0: cosine series
void Elp_main_sum_and_derv(int n, int ** &i_main, double * &A_main, Elp_args &args, 
               Elp_args &args_dot, int dist, double &sum, double &sum_dot) {
    int i;
    sum = 0.0;
    sum_dot = 0.0;
    double phase, phase_dot;
    if (dist==0) {
       // sine series 
       for (i=0; i<n; i++) {
          phase = i_main[i][0]*args.D + i_main[i][1]*args.F + i_main[i][2]*args.L + 
                  i_main[i][3]*args.Lp;
          phase_dot = i_main[i][0]*args_dot.D + i_main[i][1]*args_dot.F + 
                      i_main[i][2]*args_dot.L + i_main[i][3]*args_dot.Lp;
          sum += A_main[i]*sin(phase);
          sum_dot += A_main[i]*cos(phase)*phase_dot;
       }
    } else {
       // cosine series
       for (i=0; i<n; i++) {
          phase = i_main[i][0]*args.D + i_main[i][1]*args.F + i_main[i][2]*args.L +
                  i_main[i][3]*args.Lp;
          phase_dot = i_main[i][0]*args_dot.D + i_main[i][1]*args_dot.F +
                      i_main[i][2]*args_dot.L + i_main[i][3]*args_dot.Lp;
          sum += A_main[i]*cos(phase);
          sum_dot -= A_main[i]*sin(phase)*phase_dot;
       }
    }
}

// Sum the ELP/MPP02 series and its time derivative for perturbations 
void Elp_perturbation_sum_and_derv(int n, int ** &i_pert, double * &A_pert, double * &ph_pert,
                            Elp_args &args, Elp_args &args_dot, double &sum, double &sum_dot) {
    int i;
    sum = 0.0;
    sum_dot = 0.0;
    double phase, phase_dot;
    for (i=0; i<n; i++) {
       phase = ph_pert[i] + i_pert[i][0]*args.D + i_pert[i][1]*args.F + 
               i_pert[i][2]*args.L + i_pert[i][3]*args.Lp + i_pert[i][4]*args.Me + 
               i_pert[i][5]*args.Ve + i_pert[i][6]*args.EM + i_pert[i][7]*args.Ma + 
               i_pert[i][8]*args.Ju + i_pert[i][9]*args.Sa + i_pert[i][10]*args.Ur +
               i_pert[i][11]*args.Ne + i_pert[i][12]*args.zeta;
       phase_dot = i_pert[i][0]*args_dot.D + i_pert[i][1]*args_dot.F +
               i_pert[i][2]*args_dot.L + i_pert[i][3]*args_dot.Lp + i_pert[i][4]*args_dot.Me +
               i_pert[i][5]*args_dot.Ve + i_pert[i][6]*args_dot.EM + i_pert[i][7]*args_dot.Ma +
               i_pert[i][8]*args_dot.Ju + i_pert[i][9]*args_dot.Sa + i_pert[i][10]*args_dot.Ur +
               i_pert[i][11]*args_dot.Ne + i_pert[i][12]*args_dot.zeta;
       sum += A_pert[i]*sin(phase);
       sum_dot += A_pert[i]*cos(phase)*phase_dot;
    }
}

// Calculate the Moon's geocentric X,Y,Z coordinates with respect to 
// J2000.0 mean ecliptic and equinox.
// T is the TDB Julian century from J2000.0 = (TBD JD - 2451545)/36525
void getX2000(double T,  Elp_paras &paras, Elp_coefs &coefs, 
              double &X, double &Y, double &Z) {
  double T2 = T*T;
  double T3 = T*T2;
  double T4 = T2*T2;
  double T5 = T2*T3;
  Elp_args args;
  compute_Elp_arguments(T, paras, args);

  // Sum the ELP/MPP02 series
  // main problem series
  double main_long = Elp_main_sum(coefs.n_main_long, coefs.i_main_long, 
                                  coefs.A_main_long, args, 0);
  double main_lat = Elp_main_sum(coefs.n_main_lat, coefs.i_main_lat,
                                  coefs.A_main_lat, args, 0);
  double main_dist = Elp_main_sum(coefs.n_main_dist, coefs.i_main_dist,
                                  coefs.A_main_dist, args, 1);
  // perturbation, longitude
  double pert_longT0 = Elp_perturbation_sum(coefs.n_pert_longT0, coefs.i_pert_longT0, 
                                            coefs.A_pert_longT0, coefs.ph_pert_longT0, args);
  double pert_longT1 = Elp_perturbation_sum(coefs.n_pert_longT1, coefs.i_pert_longT1,
                                            coefs.A_pert_longT1, coefs.ph_pert_longT1, args);
  double pert_longT2 = Elp_perturbation_sum(coefs.n_pert_longT2, coefs.i_pert_longT2,
                                            coefs.A_pert_longT2, coefs.ph_pert_longT2, args);
  double pert_longT3 = Elp_perturbation_sum(coefs.n_pert_longT3, coefs.i_pert_longT3,
                                            coefs.A_pert_longT3, coefs.ph_pert_longT3, args);
  // perturbation, latitude
  double pert_latT0 = Elp_perturbation_sum(coefs.n_pert_latT0, coefs.i_pert_latT0,
                                            coefs.A_pert_latT0, coefs.ph_pert_latT0, args);
  double pert_latT1 = Elp_perturbation_sum(coefs.n_pert_latT1, coefs.i_pert_latT1,
                                            coefs.A_pert_latT1, coefs.ph_pert_latT1, args);
  double pert_latT2 = Elp_perturbation_sum(coefs.n_pert_latT2, coefs.i_pert_latT2,
                                            coefs.A_pert_latT2, coefs.ph_pert_latT2, args);
  // perturbation, distance
  double pert_distT0 = Elp_perturbation_sum(coefs.n_pert_distT0, coefs.i_pert_distT0,
                                            coefs.A_pert_distT0, coefs.ph_pert_distT0, args);
  double pert_distT1 = Elp_perturbation_sum(coefs.n_pert_distT1, coefs.i_pert_distT1,
                                            coefs.A_pert_distT1, coefs.ph_pert_distT1, args);
  double pert_distT2 = Elp_perturbation_sum(coefs.n_pert_distT2, coefs.i_pert_distT2,
                                            coefs.A_pert_distT2, coefs.ph_pert_distT2, args);
  double pert_distT3 = Elp_perturbation_sum(coefs.n_pert_distT3, coefs.i_pert_distT3,
                                            coefs.A_pert_distT3, coefs.ph_pert_distT3, args);

  // Moon's longitude, latitude and distance
  double longM = args.W1 + main_long + pert_longT0 + mod2pi(pert_longT1*T) + 
                 mod2pi(pert_longT2*T2) + mod2pi(pert_longT3*T3);
  double latM  = main_lat + pert_latT0 + mod2pi(pert_latT1*T) + mod2pi(pert_latT2*T2);
  const double ra0 = 384747.961370173/384747.980674318;
  double r = ra0*(main_dist +  pert_distT0 + pert_distT1*T + pert_distT2*T2 + pert_distT3*T3);
  double x0 = r*cos(longM)*cos(latM);
  double y0 = r*sin(longM)*cos(latM);
  double z0 = r*sin(latM);

  // Precession matrix
  double P = 0.10180391e-4*T + 0.47020439e-6*T2 - 0.5417367e-9*T3 
             - 0.2507948e-11*T4 + 0.463486e-14*T5;
  double Q = -0.113469002e-3*T + 0.12372674e-6*T2 + 0.12654170e-8*T3 
             - 0.1371808e-11*T4 - 0.320334e-14*T5;
  double sq = sqrt(1 - P*P - Q*Q);
  double p11 = 1 - 2*P*P;
  double p12 = 2*P*Q;
  double p13 = 2*P*sq;
  double p21 = 2*P*Q;
  double p22 = 1-2*Q*Q;
  double p23 = -2*Q*sq;
  double p31 = -2*P*sq;
  double p32 = 2*Q*sq;
  double p33 = 1 - 2*P*P - 2*Q*Q;

  // Finally, components of position vector wrt J2000.0 mean ecliptic and equinox
  X = p11*x0 + p12*y0 + p13*z0;
  Y = p21*x0 + p22*y0 + p23*z0;
  Z = p31*x0 + p32*y0 + p33*z0;
}

// Calculate the Moon's geocentric X,Y,Z coordinates and their time derivatives 
// with respect to J2000.0 mean ecliptic and equinox.
// T is the TDB Julian century from J2000.0 = (TBD JD - 2451545)/36525
// Xvec[0] = X, Xvec[1] = Y, Xvec[2] = Z, Xvec[3] = Xdot, Xvec[4] = Ydot, Xvec[5] = Zdot
void getX2000_Xdot2000(double T,  Elp_paras &paras, Elp_coefs &coefs, double (&Xvec)[6]) {
  double T2 = T*T;
  double T3 = T*T2;
  double T4 = T2*T2;
  double T5 = T2*T3;
  double fac = 1.0/36525;
  Elp_args args, args_dot;
  compute_Elp_arguments(T, paras, args);
  compute_Elp_arguments_dot(T, paras, args_dot);

  // Sum the ELP/MPP02 series
  // main problem series
  double main_long, main_long_dot, main_lat, main_lat_dot, main_dist, main_dist_dot;
  Elp_main_sum_and_derv(coefs.n_main_long, coefs.i_main_long,
                        coefs.A_main_long, args, args_dot, 0, main_long, main_long_dot);
  Elp_main_sum_and_derv(coefs.n_main_lat, coefs.i_main_lat, 
                        coefs.A_main_lat, args, args_dot, 0, main_lat, main_lat_dot);
  Elp_main_sum_and_derv(coefs.n_main_dist, coefs.i_main_dist,
                        coefs.A_main_dist, args, args_dot, 1, main_dist, main_dist_dot);
  // perturbation, longitude
  double pert_longT0, pert_longT0_dot, pert_longT1, pert_longT1_dot, 
         pert_longT2, pert_longT2_dot, pert_longT3, pert_longT3_dot;
  Elp_perturbation_sum_and_derv(coefs.n_pert_longT0, coefs.i_pert_longT0, 
                                coefs.A_pert_longT0, coefs.ph_pert_longT0, 
                                args, args_dot, pert_longT0, pert_longT0_dot);
  Elp_perturbation_sum_and_derv(coefs.n_pert_longT1, coefs.i_pert_longT1,
                                coefs.A_pert_longT1, coefs.ph_pert_longT1, 
                                args, args_dot, pert_longT1, pert_longT1_dot);
  Elp_perturbation_sum_and_derv(coefs.n_pert_longT2, coefs.i_pert_longT2,
                                coefs.A_pert_longT2, coefs.ph_pert_longT2, 
                                args, args_dot, pert_longT2, pert_longT2_dot);
  Elp_perturbation_sum_and_derv(coefs.n_pert_longT3, coefs.i_pert_longT3,
                                coefs.A_pert_longT3, coefs.ph_pert_longT3,
                                args, args_dot, pert_longT3, pert_longT3_dot);
  // perturbation, latitude
  double pert_latT0, pert_latT0_dot, pert_latT1, pert_latT1_dot, pert_latT2, pert_latT2_dot;
  Elp_perturbation_sum_and_derv(coefs.n_pert_latT0, coefs.i_pert_latT0,
                                coefs.A_pert_latT0, coefs.ph_pert_latT0, 
                                args, args_dot, pert_latT0, pert_latT0_dot);
  Elp_perturbation_sum_and_derv(coefs.n_pert_latT1, coefs.i_pert_latT1,
                                coefs.A_pert_latT1, coefs.ph_pert_latT1,
                                args, args_dot, pert_latT1, pert_latT1_dot);
  Elp_perturbation_sum_and_derv(coefs.n_pert_latT2, coefs.i_pert_latT2,
                                coefs.A_pert_latT2, coefs.ph_pert_latT2,
                                args, args_dot, pert_latT2, pert_latT2_dot);
  // perturbation, distance
  double pert_distT0, pert_distT0_dot, pert_distT1, pert_distT1_dot, 
         pert_distT2, pert_distT2_dot, pert_distT3, pert_distT3_dot;
  Elp_perturbation_sum_and_derv(coefs.n_pert_distT0, coefs.i_pert_distT0,
                                coefs.A_pert_distT0, coefs.ph_pert_distT0,
                                args, args_dot, pert_distT0, pert_distT0_dot);
  Elp_perturbation_sum_and_derv(coefs.n_pert_distT1, coefs.i_pert_distT1,
                                coefs.A_pert_distT1, coefs.ph_pert_distT1,
                                args, args_dot, pert_distT1, pert_distT1_dot);
  Elp_perturbation_sum_and_derv(coefs.n_pert_distT2, coefs.i_pert_distT2,
                                coefs.A_pert_distT2, coefs.ph_pert_distT2,
                                args, args_dot, pert_distT2, pert_distT2_dot);
  Elp_perturbation_sum_and_derv(coefs.n_pert_distT3, coefs.i_pert_distT3,
                                coefs.A_pert_distT3, coefs.ph_pert_distT3,
                                args, args_dot, pert_distT3, pert_distT3_dot);

  // Moon's longitude, latitude, distance and their time derivatives
  double longM = args.W1 + main_long + pert_longT0 + mod2pi(pert_longT1*T) + 
                 mod2pi(pert_longT2*T2) + mod2pi(pert_longT3*T3);
  double latM  = main_lat + pert_latT0 + mod2pi(pert_latT1*T) + mod2pi(pert_latT2*T2);
  const double ra0 = 384747.961370173/384747.980674318;
  double r = ra0*(main_dist +  pert_distT0 + pert_distT1*T + pert_distT2*T2 + pert_distT3*T3);
  double longM_dot = args_dot.W1 + main_long_dot + pert_longT0_dot + pert_longT1_dot*T +
                     pert_longT2_dot*T2 + pert_longT3_dot*T3 + 
                     fac*(pert_longT1 + 2*T*pert_longT2 + 3*T2*pert_longT3);
  double latM_dot = main_lat_dot + pert_latT0_dot + pert_latT1_dot*T + pert_latT2_dot*T2 + 
                    fac*(pert_latT1 + 2*T*pert_latT2);
  double r_dot = ra0*(main_dist_dot +  pert_distT0_dot + pert_distT1_dot*T + 
                      pert_distT2_dot*T2 + pert_distT3_dot*T3 + fac*(pert_distT1 + 
                      2*T*pert_distT2 + 3*T2*pert_distT3) );
  double cV = cos(longM), sV = sin(longM), cU = cos(latM), sU = sin(latM);
  double x0 = r*cV*cU, y0 = r*sV*cU, z0 = r*sU;
  double x0_dot = r_dot*cV*cU - r*sV*cU*longM_dot - r*cV*sU*latM_dot;
  double y0_dot = r_dot*sV*cU + r*cV*cU*longM_dot - r*sV*sU*latM_dot;
  double z0_dot = r_dot*sU + r*cU*latM_dot;

  // Precession matrix and its time derivatives
  double P = 0.10180391e-4*T + 0.47020439e-6*T2 - 0.5417367e-9*T3 
             - 0.2507948e-11*T4 + 0.463486e-14*T5;
  double Q = -0.113469002e-3*T + 0.12372674e-6*T2 + 0.12654170e-8*T3 
             - 0.1371808e-11*T4 - 0.320334e-14*T5;
  double sq = sqrt(1 - P*P - Q*Q);
  double p11 = 1 - 2*P*P;
  double p12 = 2*P*Q;
  double p13 = 2*P*sq;
  double p21 = p12;
  double p22 = 1-2*Q*Q;
  double p23 = -2*Q*sq;
  double p31 = -p13;
  double p32 = -p23;
  double p33 = 1 - 2*P*P - 2*Q*Q;
  double P_dot = fac*(0.10180391e-4 + 0.94040878e-6*T - 1.6252101e-9*T2
             - 1.0031792e-11*T3 + 2.31743e-14*T4);
  double Q_dot = fac*(-0.113469002e-3 + 0.24745348e-6*T + 0.3796251e-8*T2 
                - 0.5487232e-11*T3 - 1.60167e-14*T4);
  double sq_dot = -(P*P_dot + Q*Q_dot)/sq;
  double p11_dot = -4*P*P_dot;
  double p12_dot = 2*(P_dot*Q + P*Q_dot);
  double p13_dot = 2*(P_dot*sq + P*sq_dot);
  double p21_dot = p12_dot;
  double p22_dot = -4*Q*Q_dot;
  double p23_dot = -2*(Q_dot*sq + Q*sq_dot);
  double p31_dot = -p13_dot;
  double p32_dot = -p23_dot;
  double p33_dot = p11_dot + p22_dot;

  // Finally, components of position and vector vector wrt J2000.0 mean ecliptic and equinox
  Xvec[0] = p11*x0 + p12*y0 + p13*z0;
  Xvec[1] = p21*x0 + p22*y0 + p23*z0;
  Xvec[2] = p31*x0 + p32*y0 + p33*z0;
  Xvec[3] = p11*x0_dot + p12*y0_dot + p13*z0_dot + p11_dot*x0 + p12_dot*y0 + p13_dot*z0;
  Xvec[4] = p21*x0_dot + p22*y0_dot + p23*z0_dot + p21_dot*x0 + p22_dot*y0 + p23_dot*z0;
  Xvec[5] = p31*x0_dot + p32*y0_dot + p33*z0_dot + p31_dot*x0 + p32_dot*y0 + p33_dot*z0;
}
