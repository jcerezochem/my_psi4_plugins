#include <cmath>

// CONSTANTS (NIST)
const double  PI      = 3.1415926535897932; // !
const double  NAv     = 6.02214129e23;      // ! Avogadro number
const double  SL      = 2.99792458e8;       // ! Speed of light
const double  plank   = 6.62606957e-34;     // ! Planck constant
const double  plankbar= 1.054571726e-34;    // ! Planck constant over 2PI
const double  kboltz  = 1.3806488e-23;      // ! Boltzman constant
const double  atmass  = 1.660538921e-27;    // ! Atomic mass
const double  cvelau  = 137.0369;           //
                   
// CONVERSION FACTORS
const double  BOHRtoANGS= 5.2917721092e-1;  // !(NIST)
const double  AMUtoKG   = 1.66053873e-27;   //
const double  AMUtoAU   = 1.82288839e3;     //
const double  AUtoKG    = 9.10938291e-31;   // !(NIST)
const double  BOHRtoM   = 5.2917721092e-11; // !(NIST)
const double  BOHRtoNM  = 5.2917721092e-2;  // !(NIST)
const double  AMStoM    = 1.e-10;           // ! exact by definition
const double  ANGStoM   = 1.e-10;           // ! exact by definition
const double  HARTtoJ   = 4.35974434e-18;   // !(NIST)
const double  HtoKCALM  = 627.5095;       //
const double  CALtoJ    = 4.184;            //
const double  HtoeV     = 27.21138505;      // !(NIST)
// FCclasses:                  
const double  autoev    = 27.2113961;     //
const double  autown    = 2.1947463068e5;   // ! Freq from AU to cm-1
const double  evtown    = 8065.5446811132;  // ! Energy from eV to cm-1
const double  autofs    = 2.4189e-2;        //
const double  fstoev    = 4.135667516;      //
const double  autoang   = 0.5291771;      //
const double  pmass     = 1.007825;       //
const double  peratio   = 1836.1515;      //
const double  autoamu   = 1.e0/(peratio/pmass); //
const double  facabs    = 703.300;       //
const double  facecd    = 20.5288;       //
const double  facemi    = 4./(3.*pow(cvelau,3)*pow(autoev,3)*autofs*1.e-6);