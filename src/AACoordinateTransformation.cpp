/*
Module : AACOORDINATETRANSFORMATION.CPP
Purpose: Implementation for the algorithms which convert between the various celestial coordinate systems
Created: PJN / 29-12-2003
History: PJN / 14-02-2004 1. Fixed a "minus zero" bug in the function CAACoordinateTransformation::DMSToDegrees.
                          The sign of the value is now taken explicitly from the new bPositive boolean
                          parameter. Thanks to Patrick Wallace for reporting this problem.
         PJN / 02-06-2005 1. Most of the angular conversion functions have now been reimplemented as simply
                          numeric constants. All of the AA+ code has also been updated to use these new constants.
         PJN / 25-01-2007 1. Fixed a minor compliance issue with GCC in the AACoordinateTransformation.h to do
                          with the declaration of various methods. Thanks to Mathieu Peyréga for reporting this
                          issue.

Copyright (c) 2003 - 2013 by PJ Naughter (Web: www.naughter.com, Email: pjna@naughter.com)

All rights reserved.

Copyright / Usage Details:

You are allowed to include the source code in any product (commercial, shareware, freeware or otherwise) 
when your product is released in binary form. You are allowed to modify the source code in any way you want 
except you cannot modify the copyright details at the top of each module. If you want to distribute source 
code with your application, then you are only allowed to distribute versions released by the author. This is 
to maintain a single distribution point for the source code. 

*/


//////////////////////// Includes /////////////////////////////////////////////

#include "stdafx.h"
#include "AACoordinateTransformation.h"
#include <cmath>
#include <cassert>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


/////////////////////// Implementation ////////////////////////////////////////

CAA2DCoordinate CAACoordinateTransformation::Equatorial2Ecliptic(double Alpha, double Delta, double Epsilon)
{
  Alpha = HoursToRadians(Alpha);
  Delta = DegreesToRadians(Delta);
  Epsilon = DegreesToRadians(Epsilon);

  CAA2DCoordinate Ecliptic;
  Ecliptic.X = RadiansToDegrees(atan2(sin(Alpha)*cos(Epsilon) + tan(Delta)*sin(Epsilon), cos(Alpha)));
  if (Ecliptic.X < 0)
    Ecliptic.X += 360;
  Ecliptic.Y = RadiansToDegrees(asin(sin(Delta)*cos(Epsilon) - cos(Delta)*sin(Epsilon)*sin(Alpha)));

  return Ecliptic;
}

CAA2DCoordinate CAACoordinateTransformation::Ecliptic2Equatorial(double Lambda, double Beta, double Epsilon)
{
  Lambda = DegreesToRadians(Lambda);
  Beta = DegreesToRadians(Beta);
  Epsilon = DegreesToRadians(Epsilon);

  CAA2DCoordinate Equatorial;
  Equatorial.X = RadiansToHours(atan2(sin(Lambda)*cos(Epsilon) - tan(Beta)*sin(Epsilon), cos(Lambda)));
  if (Equatorial.X < 0)
    Equatorial.X += 24;
  Equatorial.Y = RadiansToDegrees(asin(sin(Beta)*cos(Epsilon) + cos(Beta)*sin(Epsilon)*sin(Lambda)));
  
  return Equatorial;
}

CAA2DCoordinate CAACoordinateTransformation::Equatorial2Horizontal(double LocalHourAngle, double Delta, double Latitude)
{
  LocalHourAngle = HoursToRadians(LocalHourAngle);
  Delta = DegreesToRadians(Delta);
  Latitude = DegreesToRadians(Latitude);

  CAA2DCoordinate Horizontal;
  Horizontal.X = RadiansToDegrees(atan2(sin(LocalHourAngle), cos(LocalHourAngle)*sin(Latitude) - tan(Delta)*cos(Latitude)));
  if (Horizontal.X < 0)
    Horizontal.X += 360;
  Horizontal.Y = RadiansToDegrees(asin(sin(Latitude)*sin(Delta) + cos(Latitude)*cos(Delta)*cos(LocalHourAngle)));
    
  return Horizontal;
}

CAA2DCoordinate CAACoordinateTransformation::Horizontal2Equatorial(double Azimuth, double Altitude, double Latitude)
{
  //Convert from degress to radians
  Azimuth = DegreesToRadians(Azimuth);
  Altitude = DegreesToRadians(Altitude);
  Latitude = DegreesToRadians(Latitude);

  CAA2DCoordinate Equatorial;
  Equatorial.X = RadiansToHours(atan2(sin(Azimuth), cos(Azimuth)*sin(Latitude) + tan(Altitude)*cos(Latitude)));
  if (Equatorial.X < 0)
    Equatorial.X += 24;
  Equatorial.Y = RadiansToDegrees(asin(sin(Latitude)*sin(Altitude) - cos(Latitude)*cos(Altitude)*cos(Azimuth)));
  
  return Equatorial;
}

CAA2DCoordinate CAACoordinateTransformation::Equatorial2Galactic(double Alpha, double Delta)
{
  Alpha = 192.25 - HoursToDegrees(Alpha);
  Alpha = DegreesToRadians(Alpha);
  Delta = DegreesToRadians(Delta);

  CAA2DCoordinate Galactic;
  Galactic.X = RadiansToDegrees(atan2(sin(Alpha), cos(Alpha)*sin(DegreesToRadians(27.4)) - tan(Delta)*cos(DegreesToRadians(27.4))));
  Galactic.X = 303 - Galactic.X;
  if (Galactic.X >= 360)
    Galactic.X -= 360;
  Galactic.Y = RadiansToDegrees(asin(sin(Delta)*sin(DegreesToRadians(27.4)) + cos(Delta)*cos(DegreesToRadians(27.4))*cos(Alpha)));
    
  return Galactic;
}

CAA2DCoordinate CAACoordinateTransformation::Galactic2Equatorial(double l, double b)
{
  l -= 123;
  l = DegreesToRadians(l);
  b = DegreesToRadians(b);

  CAA2DCoordinate Equatorial;
  Equatorial.X = RadiansToDegrees(atan2(sin(l), cos(l)*sin(DegreesToRadians(27.4)) - tan(b)*cos(DegreesToRadians(27.4))));
  Equatorial.X += 12.25;
  if (Equatorial.X < 0)
    Equatorial.X += 360;
  Equatorial.X = DegreesToHours(Equatorial.X);
  Equatorial.Y = RadiansToDegrees(asin(sin(b)*sin(DegreesToRadians(27.4)) + cos(b)*cos(DegreesToRadians(27.4))*cos(l)));
    
  return Equatorial;
}

double CAACoordinateTransformation::DMSToDegrees(double Degrees, double Minutes, double Seconds, bool bPositive)
{
  //validate our parameters
  if (!bPositive)
  {
    assert(Degrees >= 0);  //All parameters should be non negative if the "bPositive" parameter is false
    assert(Minutes >= 0);
    assert(Seconds >= 0);
  }

  if (bPositive)
    return Degrees + Minutes/60 + Seconds/3600;
  else
    return -Degrees - Minutes/60 - Seconds/3600;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//////#########################################################################################//////
//////#########################################################################################//////
//////#########################################################################################//////
//////#########################################################################################//////
//////#########################################################################################//////
/////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP CAACoordinateTransformation_Equatorial2Ecliptic(SEXP Alpha_, SEXP Delta_, SEXP Epsilon_)
{ 
  double Alpha     = as<double>(Alpha_  );
  double Delta     = as<double>(Delta_  );
  double Epsilon   = as<double>(Epsilon_);
  
  Alpha = CAACoordinateTransformation::HoursToRadians(Alpha);
  Delta = CAACoordinateTransformation::DegreesToRadians(Delta);
  Epsilon = CAACoordinateTransformation::DegreesToRadians(Epsilon);

  CAA2DCoordinate Ecliptic;
  Ecliptic.X = CAACoordinateTransformation::RadiansToDegrees(atan2(sin(Alpha)*cos(Epsilon) + tan(Delta)*sin(Epsilon), cos(Alpha)));
  if (Ecliptic.X < 0)
    Ecliptic.X += 360;
  Ecliptic.Y = CAACoordinateTransformation::RadiansToDegrees(asin(sin(Delta)*cos(Epsilon) - cos(Delta)*sin(Epsilon)*sin(Alpha)));

  return List::create(Named("Ecliptic.X") = Ecliptic.X, Named("Ecliptic.Y") = Ecliptic.Y);
}

RcppExport SEXP CAACoordinateTransformation_Ecliptic2Equatorial(SEXP Lambda_, SEXP Beta_, SEXP Epsilon_)
{ 
  double Lambda  = as<double>(Lambda_);
  double Beta    = as<double>(Beta_);
  double Epsilon = as<double>(Epsilon_);
  
  Lambda = CAACoordinateTransformation::DegreesToRadians(Lambda);
  Beta = CAACoordinateTransformation::DegreesToRadians(Beta);
  Epsilon = CAACoordinateTransformation::DegreesToRadians(Epsilon);

  CAA2DCoordinate Equatorial;
  Equatorial.X = CAACoordinateTransformation::RadiansToHours(atan2(sin(Lambda)*cos(Epsilon) - tan(Beta)*sin(Epsilon), cos(Lambda)));
  if (Equatorial.X < 0)
    Equatorial.X += 24;
  Equatorial.Y = CAACoordinateTransformation::RadiansToDegrees(asin(sin(Beta)*cos(Epsilon) + cos(Beta)*sin(Epsilon)*sin(Lambda)));
  
  return List::create(Named("Equatorial.X") = Equatorial.X, Named("Equatorial.Y") = Equatorial.Y);
}

RcppExport SEXP CAACoordinateTransformation_Equatorial2Horizontal(SEXP LocalHourAngle_, SEXP Delta_, SEXP Latitude_)
{ 
  double LocalHourAngle = as<double>(LocalHourAngle_);
  double Delta          = as<double>(Delta_);
  double Latitude       = as<double>(Latitude_);
  
  LocalHourAngle = CAACoordinateTransformation::HoursToRadians(LocalHourAngle);
  Delta = CAACoordinateTransformation::DegreesToRadians(Delta);
  Latitude = CAACoordinateTransformation::DegreesToRadians(Latitude);

  CAA2DCoordinate Horizontal;
  Horizontal.X = CAACoordinateTransformation::RadiansToDegrees(atan2(sin(LocalHourAngle), cos(LocalHourAngle)*sin(Latitude) - tan(Delta)*cos(Latitude)));
  if (Horizontal.X < 0)
    Horizontal.X += 360;
  Horizontal.Y = CAACoordinateTransformation::RadiansToDegrees(asin(sin(Latitude)*sin(Delta) + cos(Latitude)*cos(Delta)*cos(LocalHourAngle)));
    
  return List::create(Named("Horizontal.X") = Horizontal.X, Named("Horizontal.Y") = Horizontal.Y);
}

RcppExport SEXP CAACoordinateTransformation_Horizontal2Equatorial(SEXP Azimuth_, SEXP Altitude_, SEXP Latitude_)
{
double Azimuth  = as<double>(Azimuth_);
double Altitude = as<double>(Altitude_);
double Latitude = as<double>(Latitude_);

  //Convert from degress to radians
  Azimuth = CAACoordinateTransformation::DegreesToRadians(Azimuth);
  Altitude = CAACoordinateTransformation::DegreesToRadians(Altitude);
  Latitude = CAACoordinateTransformation::DegreesToRadians(Latitude);

  CAA2DCoordinate Equatorial;
  Equatorial.X = CAACoordinateTransformation::RadiansToHours(atan2(sin(Azimuth), cos(Azimuth)*sin(Latitude) + tan(Altitude)*cos(Latitude)));
  if (Equatorial.X < 0)
    Equatorial.X += 24;
  Equatorial.Y = CAACoordinateTransformation::RadiansToDegrees(asin(sin(Latitude)*sin(Altitude) - cos(Latitude)*cos(Altitude)*cos(Azimuth)));
  
  return  List::create(Named("Equatorial.X") = Equatorial.X, Named("Equatorial.Y") = Equatorial.Y);
}

RcppExport SEXP CAACoordinateTransformation_Equatorial2Galactic(SEXP Alpha_, SEXP Delta_)
{
  double Alpha = as<double>(Alpha_);
  double Delta = as<double>(Delta_);
  
  Alpha = 192.25 - CAACoordinateTransformation::HoursToDegrees(Alpha);
  Alpha = CAACoordinateTransformation::DegreesToRadians(Alpha);
  Delta = CAACoordinateTransformation::DegreesToRadians(Delta);

  CAA2DCoordinate Galactic;
  Galactic.X = CAACoordinateTransformation::RadiansToDegrees(atan2(sin(Alpha), cos(Alpha)*sin(CAACoordinateTransformation::DegreesToRadians(27.4)) - tan(Delta)*cos(CAACoordinateTransformation::DegreesToRadians(27.4))));
  Galactic.X = 303 - Galactic.X;
  if (Galactic.X >= 360)
    Galactic.X -= 360;
  Galactic.Y = CAACoordinateTransformation::RadiansToDegrees(asin(sin(Delta)*sin(CAACoordinateTransformation::DegreesToRadians(27.4)) + cos(Delta)*cos(CAACoordinateTransformation::DegreesToRadians(27.4))*cos(Alpha)));
    
  return  List::create(Named("Galactic.X") = Galactic.X, Named("Galactic.Y") = Galactic.Y);
}

RcppExport SEXP CAACoordinateTransformation_Galactic2Equatorial(SEXP l_, SEXP b_)
{ 
  double l = as<double>(l_);
  double b = as<double>(b_);

  l -= 123;
  l = CAACoordinateTransformation::DegreesToRadians(l);
  b = CAACoordinateTransformation::DegreesToRadians(b);

  CAA2DCoordinate Equatorial;
  Equatorial.X = CAACoordinateTransformation::RadiansToDegrees(atan2(sin(l), cos(l)*sin(CAACoordinateTransformation::DegreesToRadians(27.4)) - tan(b)*cos(CAACoordinateTransformation::DegreesToRadians(27.4))));
  Equatorial.X += 12.25;
  if (Equatorial.X < 0)
    Equatorial.X += 360;
  Equatorial.X = CAACoordinateTransformation::DegreesToHours(Equatorial.X);
  Equatorial.Y = CAACoordinateTransformation::RadiansToDegrees(asin(sin(b)*sin(CAACoordinateTransformation::DegreesToRadians(27.4)) + cos(b)*cos(CAACoordinateTransformation::DegreesToRadians(27.4))*cos(l)));
    
  return List::create(Named("Equatorial.X") = Equatorial.X, Named("Equatorial.Y") = Equatorial.Y);
}

RcppExport SEXP CAACoordinateTransformation_DMSToDegrees(SEXP Degrees_, SEXP Minutes_, SEXP Seconds_, SEXP bPositive_)
{ 
  double Degrees = as<double>(Degrees_); 
  double Minutes = as<double>(Minutes_); 
  double Seconds = as<double>(Seconds_); 
  bool bPositive = as<bool>(bPositive_);
  
  //validate our parameters
  if (!bPositive)
  {
    assert(Degrees >= 0);  //All parameters should be non negative if the "bPositive" parameter is false
    assert(Minutes >= 0);
    assert(Seconds >= 0);
  }

  if (bPositive)
    return wrap(Degrees + Minutes/60 + Seconds/3600);
  else
    return wrap(-Degrees - Minutes/60 - Seconds/3600);
}
