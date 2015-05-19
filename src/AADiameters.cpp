/*
Module : AADIAMETERS.CPP
Purpose: Implementation for the algorithms for the semi diameters of the Sun, Moon, Planets and Asteroids
Created: PJN / 15-01-2004
History: None

Copyright (c) 2004 - 2013 by PJ Naughter (Web: www.naughter.com, Email: pjna@naughter.com)

All rights reserved.

Copyright / Usage Details:

You are allowed to include the source code in any product (commercial, shareware, freeware or otherwise) 
when your product is released in binary form. You are allowed to modify the source code in any way you want 
except you cannot modify the copyright details at the top of each module. If you want to distribute source 
code with your application, then you are only allowed to distribute versions released by the author. This is 
to maintain a single distribution point for the source code. 

*/


//////////////////// Includes /////////////////////////////////////////////////

#include "stdafx.h"
#include "AADiameters.h"
#include "AACoordinateTransformation.h"
#include "AAGlobe.h"
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


//////////////////// Implementation ///////////////////////////////////////////

double CAADiameters::SunSemidiameterA(double Delta)
{
  return 959.63/Delta;
}

double CAADiameters::MercurySemidiameterA(double Delta)
{
  return 3.34/Delta;
}

double CAADiameters::VenusSemidiameterA(double Delta)
{
  return 8.41/Delta;
}

double CAADiameters::MarsSemidiameterA(double Delta)
{
  return 4.68/Delta;
}

double CAADiameters::JupiterEquatorialSemidiameterA(double Delta)
{
  return 98.47/Delta;
}

double CAADiameters::JupiterPolarSemidiameterA(double Delta)
{
  return 91.91/Delta;
}

double CAADiameters::SaturnEquatorialSemidiameterA(double Delta)
{
  return 83.33/Delta;
}

double CAADiameters::SaturnPolarSemidiameterA(double Delta)
{
  return 74.57/Delta;
}

double CAADiameters::ApparentSaturnPolarSemidiameterA(double Delta, double B)
{
  double cosB = cos(CAACoordinateTransformation::DegreesToRadians(B));
  return SaturnPolarSemidiameterA(Delta)*sqrt(1 - 0.199197*cosB*cosB);
}

double CAADiameters::UranusSemidiameterA(double Delta)
{
  return 34.28/Delta;
}

double CAADiameters::NeptuneSemidiameterA(double Delta)
{
  return 36.56/Delta;
}

double CAADiameters::MercurySemidiameterB(double Delta)
{
  return 3.36/Delta;
}

double CAADiameters::VenusSemidiameterB(double Delta)
{
  return 8.34/Delta;
}

double CAADiameters::MarsSemidiameterB(double Delta)
{
  return 4.68/Delta;
}

double CAADiameters::JupiterEquatorialSemidiameterB(double Delta)
{
  return 98.44/Delta;
}

double CAADiameters::JupiterPolarSemidiameterB(double Delta)
{
  return 92.06/Delta;
}

double CAADiameters::SaturnEquatorialSemidiameterB(double Delta)
{
  return 82.73/Delta;
}

double CAADiameters::SaturnPolarSemidiameterB(double Delta)
{
  return 73.82/Delta;
}

double CAADiameters::ApparentSaturnPolarSemidiameterB(double Delta, double B)
{
  double cosB = cos(CAACoordinateTransformation::DegreesToRadians(B));
  return SaturnPolarSemidiameterB(Delta)*sqrt(1 - 0.203800*cosB*cosB);
}

double CAADiameters::UranusSemidiameterB(double Delta)
{
  return 35.02/Delta;
}

double CAADiameters::NeptuneSemidiameterB(double Delta)
{
  return 33.50/Delta;
}

double CAADiameters::PlutoSemidiameterB(double Delta)
{
  return 2.07/Delta;
}

double CAADiameters::GeocentricMoonSemidiameter(double Delta)
{
  return CAACoordinateTransformation::RadiansToDegrees(0.272481*6378.14/Delta)*3600;
}

double CAADiameters::TopocentricMoonSemidiameter(double DistanceDelta, double Delta, double H, double Latitude, double Height)
{
  //Convert to radians
  H = CAACoordinateTransformation::HoursToRadians(H);
  Delta = CAACoordinateTransformation::DegreesToRadians(Delta);

  double pi = asin(6378.14/DistanceDelta);
  double A = cos(Delta)*sin(H);
  double B = cos(Delta)*cos(H) - CAAGlobe::RhoCosThetaPrime(Latitude, Height)*sin(pi);
  double C = sin(Delta) - CAAGlobe::RhoSinThetaPrime(Latitude, Height)*sin(pi);
  double q = sqrt(A*A + B*B + C*C);

  double s = CAACoordinateTransformation::DegreesToRadians(GeocentricMoonSemidiameter(DistanceDelta)/3600);
  return CAACoordinateTransformation::RadiansToDegrees(asin(sin(s)/q))*3600;
}

double CAADiameters::AsteroidDiameter(double H, double A)
{
  double x = 3.12 - H/5 - 0.217147*log(A);
  return pow(10.0, x);
}

double CAADiameters::ApparentAsteroidDiameter(double Delta, double d)
{
  return 0.0013788*d/Delta;
}

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
/*#########################################################*/

RcppExport SEXP CAADiameters_SunSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(959.63/Delta));
}

RcppExport SEXP CAADiameters_MercurySemidiameterA(SEXP  Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(3.34/Delta));
}

RcppExport SEXP CAADiameters_VenusSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(8.41/Delta));
}

RcppExport SEXP CAADiameters_MarsSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(4.68/Delta));
}

RcppExport SEXP CAADiameters_JupiterEquatorialSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(98.47/Delta));
}

RcppExport SEXP CAADiameters_JupiterPolarSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(91.91/Delta));
}

RcppExport SEXP CAADiameters_SaturnEquatorialSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(83.33/Delta));
}

RcppExport SEXP CAADiameters_SaturnPolarSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(74.57/Delta));
}

RcppExport SEXP CAADiameters_ApparentSaturnPolarSemidiameterA(SEXP Delta_, SEXP B_)
{
  double Delta  = as<double>(Delta_);
  double B  = as<double>(B_);
  double cosB = cos(CAACoordinateTransformation::DegreesToRadians(B));
  return ( wrap(CAADiameters::SaturnPolarSemidiameterA(Delta)*sqrt(1 - 0.199197*cosB*cosB)));
}

RcppExport SEXP CAADiameters_UranusSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(34.28/Delta));
}

RcppExport SEXP CAADiameters_NeptuneSemidiameterA(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(36.56/Delta));
}

RcppExport SEXP CAADiameters_MercurySemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(3.36/Delta));
}

RcppExport SEXP CAADiameters_VenusSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(8.34/Delta));
}

RcppExport SEXP CAADiameters_MarsSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(4.68/Delta));
}

RcppExport SEXP CAADiameters_JupiterEquatorialSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(98.44/Delta));
}

RcppExport SEXP CAADiameters_JupiterPolarSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(92.06/Delta));
}

RcppExport SEXP CAADiameters_SaturnEquatorialSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(82.73/Delta));
}

RcppExport SEXP CAADiameters_SaturnPolarSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(73.82/Delta));
}

RcppExport SEXP CAADiameters_ApparentSaturnPolarSemidiameterB(SEXP Delta_, SEXP B_)
{
  double Delta  = as<double>(Delta_);
  double B  = as<double>(B_);
  double cosB = cos(CAACoordinateTransformation::DegreesToRadians(B));
  return ( wrap(CAADiameters::SaturnPolarSemidiameterB(Delta)*sqrt(1 - 0.203800*cosB*cosB)));
}

RcppExport SEXP CAADiameters_UranusSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(35.02/Delta));
}

RcppExport SEXP CAADiameters_NeptuneSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(33.50/Delta));
}

RcppExport SEXP CAADiameters_PlutoSemidiameterB(SEXP Delta_)
{
  double Delta  = as<double>(Delta_);
  return ( wrap(2.07/Delta));
}

RcppExport SEXP CAADiameters_GeocentricMoonSemidiameter(SEXP Delta_)
{ 
  double Delta  = as<double>(Delta_);
  return ( wrap(CAACoordinateTransformation::RadiansToDegrees(0.272481*6378.14/Delta)*3600));
}

RcppExport SEXP CAADiameters_TopocentricMoonSemidiameter(SEXP DistanceDelta_, SEXP Delta_, SEXP H_, SEXP Latitude_, SEXP Height_)
{ 
 double DistanceDelta = as<double>(DistanceDelta_); 
 double Delta         = as<double>(Delta_        );
 double H             = as<double>(H_            );
 double Latitude      = as<double>(Latitude_     );
 double Height        = as<double>(Height_       );
 
  //Convert to radians
  H = CAACoordinateTransformation::HoursToRadians(H);
  Delta = CAACoordinateTransformation::DegreesToRadians(Delta);

  double pi = asin(6378.14/DistanceDelta);
  double A = cos(Delta)*sin(H);
  double B = cos(Delta)*cos(H) - CAAGlobe::RhoCosThetaPrime(Latitude, Height)*sin(pi);
  double C = sin(Delta) - CAAGlobe::RhoSinThetaPrime(Latitude, Height)*sin(pi);
  double q = sqrt(A*A + B*B + C*C);

  double s = CAACoordinateTransformation::DegreesToRadians(CAADiameters::GeocentricMoonSemidiameter(DistanceDelta)/3600);
  return ( wrap(CAACoordinateTransformation::RadiansToDegrees(asin(sin(s)/q))*3600));
}

RcppExport SEXP CAADiameters_AsteroidDiameter(SEXP H_, SEXP A_)
{ 
  double H = as<double>(H_);
  double A = as<double>(A_);
  double x = (3.12 - H/5 - 0.217147*log(A));
  return ( wrap(pow(10.0, x)));
}

RcppExport SEXP CAADiameters_ApparentAsteroidDiameter(SEXP Delta_, SEXP d_)
{ 
  double Delta = as<double>(Delta_);
  double d = as<double>(d_);
  return ( wrap(0.0013788*d/Delta));
}
