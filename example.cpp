#include "ElpMpp02.h"

int main() {
   Elp_facs facs; 
   Elp_paras paras;
   Elp_coefs coefs;

   int corr,i,n;
   double X,Y,Z, JD0,JD,dJD, T0,T,dT;

//-----------------------------------------------------------------------
//     Compute Lunar coordinates (J2000.0 mean ecliptic and equinox):
//     5 dates from JD 2444239.5 to 2452239.5,
//     Time interval = 2000 days,
//     Parameters fitted to LLR observations.
//-----------------------------------------------------------------------
   corr = 0; // use parameters fitted to the LLR data
   setup_parameters(corr, paras, facs);
   setup_Elp_coefs(coefs, facs);
   JD0 = 2444239.5;
   dJD = 2000.0;
   n = 5;
   T0 = (JD0 - 2451545.0)/36525.0;
   dT = dJD/36525.0;
   printf("-------------------------------------------------------\n");
   printf("Lunar coordinates (J2000.0 mean ecliptic and equinox):\n");
   printf("5 dates from JD 2444239.5 to 2452239.5,\n");
   printf("Time interval = 2000 days,\n");
   printf("Parameters fitted to LLR observations.\n");
   printf("-------------------------------------------------------\n");
   for (i=0; i<n; i++) {
      JD = JD0 + i*dJD;
      T = T0 + i*dT;
      getX2000(T, paras, coefs, X,Y,Z);
      printf("JD: %f\n", JD);
      printf("X, Y, Z (km): %f   %f   %f\n\n\n", X,Y,Z);
   }

//-----------------------------------------------------------------------
//     Compute Lunar coordinates (J2000.0 mean ecliptic and equinox):
//     5 dates from JD 2500000.5 to 1700000.5,
//     Time interval = -200000 days,
//     Parameters fitted to JPL ephemeris DE405/DE406.
//-----------------------------------------------------------------------
   corr = 1; // use parameters fitted to DE405
   setup_parameters(corr, paras, facs);
   setup_Elp_coefs(coefs, facs);
   JD0 = 2500000.5;
   dJD = -2e5;
   n = 5;
   T0 = (JD0 - 2451545.0)/36525.0;
   dT = dJD/36525.0;
   printf("-------------------------------------------------------\n");
   printf("Lunar coordinates (J2000.0 mean ecliptic and equinox):\n");
   printf("5 dates from JD 2500000.5 to 1700000.5,\n");
   printf("Time interval = -200000 days,\n");
   printf("Parameters fitted to JPL ephemeris DE405/DE406.\n");
   printf("-------------------------------------------------------\n");
   for (i=0; i<n; i++) {
      JD = JD0 + i*dJD;
      T = T0 + i*dT;
      getX2000(T, paras, coefs, X,Y,Z);
      printf("JD: %f\n", JD);
      printf("X, Y, Z (km): %f   %f   %f\n\n", X,Y,Z);
   }

  return 0;
}
