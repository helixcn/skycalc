/*
Module : AAFK5.CPP
Purpose: Implementation for the algorithms to convert to the FK5 standard reference frame
Created: PJN / 29-12-2003
History: None

Copyright (c) 2003 - 2013 by PJ Naughter (Web: www.naughter.com, Email: pjna@naughter.com)

All rights reserved.

Copyright / Usage Details:

You are allowed to include the source code in any product (commercial, shareware, freeware or otherwise) 
when your product is released in binary form. You are allowed to modify the source code in any way you want 
except you cannot modify the copyright details at the top of each module. If you want to distribute source 
code with your application, then you are only allowed to distribute versions released by the author. This is 
to maintain a single distribution point for the source code. 

*/


/////////////////////// Includes //////////////////////////////////////////////

#include "stdafx.h"
#include "AAFK5.h"
#include "AACoordinateTransformation.h"
#include "AAEarth.h"
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


/////////////////////// Implementation ////////////////////////////////////////

double CAAFK5::CorrectionInLongitude(double Longitude, double Latitude, double JD)
{
  double T = (JD - 2451545) / 36525;
  double Ldash = Longitude - 1.397*T - 0.00031*T*T;

  //Convert to radians
  Ldash = CAACoordinateTransformation::DegreesToRadians(Ldash);
  Latitude = CAACoordinateTransformation::DegreesToRadians(Latitude);

  double value = -0.09033 + 0.03916*(std::cos(Ldash) + std::sin(Ldash))*tan(Latitude);
  return CAACoordinateTransformation::DMSToDegrees(0, 0, value);
}

double CAAFK5::CorrectionInLatitude(double Longitude, double JD)
{
  double T = (JD - 2451545) / 36525;
  double Ldash = Longitude - 1.397*T - 0.00031*T*T;

  //Convert to radians
  Ldash = CAACoordinateTransformation::DegreesToRadians(Ldash);

  double value = 0.03916*(std::cos(Ldash) - std::sin(Ldash));
  return CAACoordinateTransformation::DMSToDegrees(0, 0, value);
}

CAA3DCoordinate CAAFK5::ConvertVSOPToFK5J2000(const CAA3DCoordinate& value)
{
  CAA3DCoordinate result;
  result.X = value.X + 0.000000440360 * value.Y - 0.000000190919 * value.Z;
  result.Y = -0.000000479966 * value.X + 0.917482137087 * value.Y - 0.397776982902 * value.Z;
  result.Z = 0.397776982902 * value.Y + 0.917482137087 * value.Z;

  return result;
}

CAA3DCoordinate CAAFK5::ConvertVSOPToFK5B1950(const CAA3DCoordinate& value)
{
  CAA3DCoordinate result;
  result.X = 0.999925702634 * value.X + 0.012189716217 * value.Y + 0.000011134016 * value.Z;
  result.Y = -0.011179418036 * value.X + 0.917413998946 * value.Y - 0.397777041885 * value.Z;
  result.Z = -0.004859003787 * value.X + 0.397747363646 * value.Y + 0.917482111428 * value.Z;

  return result;
}

CAA3DCoordinate CAAFK5::ConvertVSOPToFK5AnyEquinox(const CAA3DCoordinate& value, double JDEquinox)
{
  double t = (JDEquinox - 2451545.0) / 36525;
  double tsquared = t*t;
  double tcubed  = tsquared * t;

  double sigma = 2306.2181*t + 0.30188*tsquared + 0.017988*tcubed;
  sigma = CAACoordinateTransformation::DegreesToRadians(CAACoordinateTransformation::DMSToDegrees(0, 0, sigma));

  double zeta = 2306.2181*t + 1.09468*tsquared + 0.018203*tcubed;
  zeta = CAACoordinateTransformation::DegreesToRadians(CAACoordinateTransformation::DMSToDegrees(0, 0, zeta));

  double phi = 2004.3109*t - 0.42665*tsquared - 0.041833*tcubed;
  phi = CAACoordinateTransformation::DegreesToRadians(CAACoordinateTransformation::DMSToDegrees(0, 0, phi));

  double cossigma = std::cos(sigma);
  double coszeta = std::cos(zeta);
  double cosphi = std::cos(phi);
  double sinsigma = std::sin(sigma);
  double sinzeta = std::sin(zeta);
  double sinphi = std::sin(phi);

  double xx = cossigma * coszeta * cosphi -sinsigma*sinzeta;
  double xy = sinsigma * coszeta + cossigma * sinzeta * cosphi;
  double xz = cossigma * sinphi;
  double yx = -cossigma * sinzeta - sinsigma * coszeta * cosphi;
  double yy = cossigma * coszeta - sinsigma * sinzeta * cosphi;
  double yz = -sinsigma * sinphi;
  double zx = -coszeta * sinphi;
  double zy = -sinzeta * sinphi;
  double zz = cosphi;

  CAA3DCoordinate result;
  result.X = xx * value.X + yx * value.Y + zx * value.Z;
  result.Y = xy * value.X + yy * value.Y + zy * value.Z;
  result.Z = xz * value.X + yz * value.Y + zz * value.Z;

  return result;
}

////########################################################################
////########################################################################
////########################################################################
////########################################################################

RcppExport SEXP CAAFK5_CorrectionInLongitude(SEXP Longitude_, SEXP Latitude_, SEXP JD_)
{
  double Longitude =as<double>(Longitude_);
  double Latitude =as<double>(Latitude_);
  double JD = as<double>(JD_);

  double T = (JD - 2451545) / 36525;
  double Ldash = Longitude - 1.397*T - 0.00031*T*T;

  //Convert to radians
  Ldash = CAACoordinateTransformation::DegreesToRadians(Ldash);
  Latitude = CAACoordinateTransformation::DegreesToRadians(Latitude);

  double value = -0.09033 + 0.03916*(std::cos(Ldash) + std::sin(Ldash))*tan(Latitude);
  return wrap(CAACoordinateTransformation::DMSToDegrees(0, 0, value));
}

RcppExport SEXP  CAAFK5_CorrectionInLatitude(SEXP Latitude_, SEXP JD_)
{ 
  double Longitude =as<double>(Latitude_);
  double JD = as<double>(JD_);
  
  double T = (JD - 2451545) / 36525;
  double Ldash = Longitude - 1.397*T - 0.00031*T*T;

  //Convert to radians
  Ldash = CAACoordinateTransformation::DegreesToRadians(Ldash);

  double value = 0.03916*(std::cos(Ldash) - std::sin(Ldash));
  return wrap(CAACoordinateTransformation::DMSToDegrees(0, 0, value));
}

RcppExport SEXP CAAFK5_ConvertVSOPToFK5J2000(SEXP value_X_, SEXP value_Y_, SEXP value_Z_)
{ 
  double value_X   = as<double>(value_X_);
  double value_Y   = as<double>(value_Y_);
  double value_Z   = as<double>(value_Z_);
  CAA3DCoordinate result;
  result.X = value_X + 0.000000440360 * value_Y - 0.000000190919 * value_Z;
  result.Y = -0.000000479966 * value_X + 0.917482137087 * value_Y - 0.397776982902 * value_Z;
  result.Z = 0.397776982902 * value_Y + 0.917482137087 * value_Z;

  return List::create(Named("X") = result.X, Named("Y") = result.Y, Named("Z") = result.Z);
}

RcppExport SEXP  CAAFK5_ConvertVSOPToFK5B1950(SEXP value_X_, SEXP value_Y_, SEXP value_Z_)
{ 
  double value_X   = as<double>(value_X_);
  double value_Y   = as<double>(value_Y_);
  double value_Z   = as<double>(value_Z_);
  CAA3DCoordinate result;
  result.X = 0.999925702634 * value_X + 0.012189716217 * value_Y + 0.000011134016 * value_Z;
  result.Y = -0.011179418036 * value_X + 0.917413998946 * value_Y - 0.397777041885 * value_Z;
  result.Z = -0.004859003787 * value_X + 0.397747363646 * value_Y + 0.917482111428 * value_Z;

  return List::create(Named("X") = result.X, Named("Y") = result.Y, Named("Z") = result.Z);
}

RcppExport SEXP  CAAFK5_ConvertVSOPToFK5AnyEquinox(SEXP value_X_, SEXP value_Y_, SEXP value_Z_, SEXP JDEquinox_)
{ 
  double value_X   = as<double>(value_X_);
  double value_Y   = as<double>(value_Y_);
  double value_Z   = as<double>(value_Z_);
  double JDEquinox  = as<double>(JDEquinox_);
  
  double t = (JDEquinox - 2451545.0) / 36525;
  double tsquared = t*t;
  double tcubed  = tsquared * t;

  double sigma = 2306.2181*t + 0.30188*tsquared + 0.017988*tcubed;
  sigma = CAACoordinateTransformation::DegreesToRadians(CAACoordinateTransformation::DMSToDegrees(0, 0, sigma));

  double zeta = 2306.2181*t + 1.09468*tsquared + 0.018203*tcubed;
  zeta = CAACoordinateTransformation::DegreesToRadians(CAACoordinateTransformation::DMSToDegrees(0, 0, zeta));

  double phi = 2004.3109*t - 0.42665*tsquared - 0.041833*tcubed;
  phi = CAACoordinateTransformation::DegreesToRadians(CAACoordinateTransformation::DMSToDegrees(0, 0, phi));

  double cossigma = std::cos(sigma);
  double coszeta = std::cos(zeta);
  double cosphi = std::cos(phi);
  double sinsigma = std::sin(sigma);
  double sinzeta = std::sin(zeta);
  double sinphi = std::sin(phi);

  double xx = cossigma * coszeta * cosphi -sinsigma*sinzeta;
  double xy = sinsigma * coszeta + cossigma * sinzeta * cosphi;
  double xz = cossigma * sinphi;
  double yx = -cossigma * sinzeta - sinsigma * coszeta * cosphi;
  double yy = cossigma * coszeta - sinsigma * sinzeta * cosphi;
  double yz = -sinsigma * sinphi;
  double zx = -coszeta * sinphi;
  double zy = -sinzeta * sinphi;
  double zz = cosphi;

  CAA3DCoordinate result;
  result.X = xx * value_X + yx * value_Y + zx * value_Z;
  result.Y = xy * value_X + yy * value_Y + zy * value_Z;
  result.Z = xz * value_X + yz * value_Y + zz * value_Z;

  return List::create(Named("X") = result.X, Named("Y") = result.Y, Named("Z") = result.Z);
}
