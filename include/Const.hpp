/////////////////////////////////////////////////////
//                                                 //
// Program Created by Arc-Pintade (Aurélien CARLE) //
// Thanks to Plut0n (Xavier Valcarce) for his help //
//                                                 //
/////////////////////////////////////////////////////

#ifndef Const_h
#define Const_h
#define Const_cxx

#include <cmath>


//____________ MMSM TOUT EN GeV ______________//

// Strong Interaction Coupling
   const double gS2 = 4*M_PI*0.1182;
   const double gS4 = gS2 * gS2;
// Weak Interaction Coupling
   const double gW2 = 4*M_PI*10e-5;
   const double gW4 = gW2 * gW2;
// Top Mass
   const double mt = 173.21;         //Incertitude ~0.71
   const double mt2 = mt * mt;
// Top Width
   const double gammat = 1.41;       //Incertitude ~ +0.19-0.15
   const double gammat2 = gammat * gammat;
// W Mass
   const double mW = 80.385;         //Incertitude ~0.015  
   const double mW2 = mW * mW;
// W Width
   const double gammaW = 2.085;      //Incertitude ~0.042
   const double gammaW2 = gammaW * gammaW;

//______________ Rotation ______________//

// Latitude
   const double latit = (46.3099/180) * M_PI;    //Incertitude ~0.0001
// Longitude
   const double longit = (6.0766/180) * M_PI;    //Incertitude ~0.0001
// Azimuth
   const double azim = (101.279/180) * M_PI;     //Incertitude ~0.003
// Tilt
   const double tilt = (0.20632/180) * M_PI;     //Incertitude ~0.00003
// Universal Frequence
   const double omega = 7.29211515e-5;           // en rad/s

// Rotation Matrix Values
   const double s1 = sin(latit);
   const double c1 = cos(latit);
   const double s2 = sin(azim);
   const double c2 = cos(azim);
   const double s3 = sin(tilt);
   const double c3 = cos(tilt);

//______________ Rotation D0______________//

// Latitude D0
   const double colatitD0 = (49.8255/180) * M_PI;    //Incertitude ~0.0001
// Azimuth D0 (Orientation of proton beam)
   const double azimD0 = (42.192/180) * M_PI;     //Incertitude ~0.003

// Rotation Matrix Values
   const double s1D0 = sin(colatitD0);
   const double c1D0 = cos(colatitD0);
   const double s2D0 = sin(azimD0);
   const double c2D0 = cos(azimD0);

// Dephasage D0/CMS
   const double deph = 6;

#endif
