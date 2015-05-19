/*
Module : AANODES.CPP
Purpose: Implementation for the algorithms which calculate passage thro the nodes
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


///////////////////////////// Includes ////////////////////////////////////////

#include "stdafx.h"
#include "AANodes.h"
#include "AACoordinateTransformation.h"
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


///////////////////////////// Implementation //////////////////////////////////

CAANodeObjectDetails CAANodes::PassageThroAscendingNode(const CAAEllipticalObjectElements& elements)
{
  double v = CAACoordinateTransformation::MapTo0To360Range(-elements.w);
  v = CAACoordinateTransformation::DegreesToRadians(v);
  double E = atan(sqrt((1 - elements.e) / (1 + elements.e)) * tan(v/2))*2;
  double M = E - elements.e*sin(E);
  M = CAACoordinateTransformation::RadiansToDegrees(M);
  double n = CAAElliptical::MeanMotionFromSemiMajorAxis(elements.a);

  CAANodeObjectDetails details;
  details.t = elements.T + M/n;
  details.radius = elements.a*(1 - elements.e*cos(E));

  return details;
}

CAANodeObjectDetails CAANodes::PassageThroDescendingNode(const CAAEllipticalObjectElements& elements)
{
  double v = CAACoordinateTransformation::MapTo0To360Range(180 - elements.w);
  v = CAACoordinateTransformation::DegreesToRadians(v);
  double E = atan(sqrt((1 - elements.e) / (1 + elements.e)) * tan(v/2))*2;
  double M = E - elements.e*sin(E);
  M = CAACoordinateTransformation::RadiansToDegrees(M);
  double n = CAAElliptical::MeanMotionFromSemiMajorAxis(elements.a);

  CAANodeObjectDetails details;
  details.t = elements.T + M/n;
  details.radius = elements.a*(1 - elements.e*cos(E));

  return details;
}

CAANodeObjectDetails CAANodes::PassageThroAscendingNode(const CAAParabolicObjectElements& elements)
{
  double v = CAACoordinateTransformation::MapTo0To360Range(-elements.w);
  v = CAACoordinateTransformation::DegreesToRadians(v);
  double s = tan(v / 2);
  double s2 = s*s;

  CAANodeObjectDetails details;
  details.t = elements.T + 27.403895*(s2*s + 3*s)*elements.q*sqrt(elements.q);
  details.radius = elements.q*(1 + s2);

  return details;
}

CAANodeObjectDetails CAANodes::PassageThroDescendingNode(const CAAParabolicObjectElements& elements)
{
  double v = CAACoordinateTransformation::MapTo0To360Range(180 - elements.w);
  v = CAACoordinateTransformation::DegreesToRadians(v);

  double s = tan(v / 2);
  double s2 = s*s;

  CAANodeObjectDetails details;
  details.t = elements.T + 27.403895*(s2*s + 3*s)*elements.q*sqrt(elements.q);
  details.radius = elements.q*(1 + s2);

  return details;
}

////////// #############################################################
////////// #############################################################
////////// #############################################################
////////// #############################################################

RcppExport SEXP CAANodes_PassageThroAscendingNode_CAAElliptical(SEXP elements_w_,SEXP elements_e_,SEXP elements_a_,SEXP elements_T_)
{ 
  double elements_w = as<double>(elements_w_);
  double elements_e = as<double>(elements_e_);
  double elements_a = as<double>(elements_a_);
  double elements_T = as<double>(elements_T_);
  
  double v = CAACoordinateTransformation::MapTo0To360Range(-elements_w);
  v = CAACoordinateTransformation::DegreesToRadians(v);
  double E = atan(sqrt((1 - elements_e) / (1 + elements_e)) * tan(v/2))*2;
  double M = E - elements_e*sin(E);
  M = CAACoordinateTransformation::RadiansToDegrees(M);
  double n = CAAElliptical::MeanMotionFromSemiMajorAxis(elements_a);

  CAANodeObjectDetails details;
  details.t = elements_T + M/n;
  details.radius = elements_a*(1 - elements_e*cos(E));

  return List::create(Named("details.t") = details.t, Named("details.radius") = details.radius);
}

RcppExport SEXP  CAANodes_PassageThroDescendingNode_CAAElliptical(SEXP elements_w_,SEXP elements_e_,SEXP elements_a_,SEXP elements_T_)
{
  double elements_w = as<double>(elements_w_);
  double elements_e = as<double>(elements_e_);
  double elements_a = as<double>(elements_a_);
  double elements_T = as<double>(elements_T_);
  
  double v = CAACoordinateTransformation::MapTo0To360Range(180 - elements_w);
  v = CAACoordinateTransformation::DegreesToRadians(v);
  double E = atan(sqrt((1 - elements_e) / (1 + elements_e)) * tan(v/2))*2;
  double M = E - elements_e*sin(E);
  M = CAACoordinateTransformation::RadiansToDegrees(M);
  double n = CAAElliptical::MeanMotionFromSemiMajorAxis(elements_a);

  CAANodeObjectDetails details;
  details.t = elements_T + M/n;
  details.radius = elements_a*(1 - elements_e*cos(E));

  return List::create(Named("details.t") = details.t, Named("details.radius") = details.radius);
}

RcppExport SEXP CAANodes_PassageThroAscendingNode_Parabolic(SEXP elements_w_, SEXP elements_T_, SEXP elements_q_)
{ 
  double elements_w = as<double>(elements_w_);
  double elements_T = as<double>(elements_T_);
  double elements_q = as<double>(elements_q_);

  double v = CAACoordinateTransformation::MapTo0To360Range(-elements_w);
  v = CAACoordinateTransformation::DegreesToRadians(v);
  double s = tan(v / 2);
  double s2 = s*s;

  CAANodeObjectDetails details;
  details.t = elements_T + 27.403895*(s2*s + 3*s)*elements_q*sqrt(elements_q);
  details.radius = elements_q*(1 + s2);

  return  List::create(Named("details.t") = details.t, Named("details.radius") = details.radius);
}

RcppExport SEXP  CAANodes_PassageThroDescendingNode_Parabolic(SEXP elements_w_, SEXP elements_T_, SEXP elements_q_)
{ 
  double elements_w = as<double>(elements_w_);
  double elements_T = as<double>(elements_T_);
  double elements_q = as<double>(elements_q_);
  
  double v = CAACoordinateTransformation::MapTo0To360Range(180 - elements_w);
  v = CAACoordinateTransformation::DegreesToRadians(v);

  double s = tan(v / 2);
  double s2 = s*s;

  CAANodeObjectDetails details;
  details.t = elements_T + 27.403895*(s2*s + 3*s)*elements_q*sqrt(elements_q);
  details.radius = elements_q*(1 + s2);

  return  List::create(Named("details.t") = details.t, Named("details.radius") = details.radius);
}
