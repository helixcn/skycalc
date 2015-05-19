/*
Module : AAPLANETARYPHENOMENA.CPP
Purpose: Implementation for the algorithms which obtain the dates of various planetary phenomena
Created: PJN / 29-12-2003
History: PJN / 08-05-2011 1. Fixed a compilation issue on GCC where size_t was undefined in various methods. Thanks to 
                          Carsten A. Arnholm and Andrew Hammond for reporting this bug.
         PJN / 18-03-2012 1. All global "g_*" tables are now const. Thanks to Roger Dahl for reporting this 
                          issue when compiling AA+ on ARM.

Copyright (c) 2003 - 2013 by PJ Naughter (Web: www.naughter.com, Email: pjna@naughter.com)

All rights reserved.

Copyright / Usage Details:

You are allowed to include the source code in any product (commercial, shareware, freeware or otherwise) 
when your product is released in binary form. You are allowed to modify the source code in any way you want 
except you cannot modify the copyright details at the top of each module. If you want to distribute source 
code with your application, then you are only allowed to distribute versions released by the author. This is 
to maintain a single distribution point for the source code. 

*/


//////////////////////////// Includes /////////////////////////////////////////

#include "stdafx.h"
#include "AAPlanetaryPhenomena.h"
#include "AACoordinateTransformation.h"
#include <cmath>
#include <cassert>
#include <cstddef>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


//////////////////////////// Macros / Defines /////////////////////////////////

struct PlanetaryPhenomenaCoefficient1
{
  double A;
  double B;
  double M0;
  double M1;
};

const PlanetaryPhenomenaCoefficient1 g_PlanetaryPhenomenaCoefficient1[] = 
{
  { 2451612.023, 115.8774771, 63.5867,  114.2088742 },
  { 2451554.084, 115.8774771, 6.4822,   114.2088742 },
  { 2451996.706, 583.921361,  82.7311,  215.513058  },
  { 2451704.746, 583.921361,  154.9745, 215.513058  },
  { 2452097.382, 779.936104,  181.9573, 48.705244   },
  { 2451707.414, 779.936104,  157.6047, 48.705244   },
  { 2451870.628, 398.884046,  318.4681, 33.140229   },
  { 2451671.186, 398.884046,  121.8980, 33.140229   },
  { 2451870.170, 378.091904,  318.0172, 12.647487   },
  { 2451681.124, 378.091904,  131.6934, 12.647487   },
  { 2451764.317, 369.656035,  213.6884, 4.333093    },
  { 2451579.489, 369.656035,  31.5219,  4.333093    },
  { 2451753.122, 367.486703,  202.6544, 2.194998    },
  { 2451569.379, 367.486703,  21.5569,  2.194998    },
};


/////////////////////////// Implementation ////////////////////////////////////

double CAAPlanetaryPhenomena::K(double Year, PlanetaryObject object, EventType type)
{
  size_t nCoefficient;
  if (object >= MARS)
  {
    assert(type == OPPOSITION || type == CONJUNCTION);

    if (type == OPPOSITION)
      nCoefficient = object*2;
    else
      nCoefficient = object*2 + 1;
  }
  else
  {
    assert(type == INFERIOR_CONJUNCTION || type == SUPERIOR_CONJUNCTION);

    if (type == INFERIOR_CONJUNCTION)
      nCoefficient = object*2;
    else
      nCoefficient = object*2 + 1;
  }
  assert(nCoefficient >= 0 && nCoefficient < sizeof(g_PlanetaryPhenomenaCoefficient1)/sizeof(PlanetaryPhenomenaCoefficient1));

  double k = (365.2425*Year + 1721060 - g_PlanetaryPhenomenaCoefficient1[nCoefficient].A) / g_PlanetaryPhenomenaCoefficient1[nCoefficient].B;
  return floor(k + 0.5);
}

double CAAPlanetaryPhenomena::Mean(double k, PlanetaryObject object, EventType type)
{
  size_t nCoefficient;
  if (object >= MARS)
  {
    assert(type == OPPOSITION || type == CONJUNCTION);

    if (type == OPPOSITION)
      nCoefficient = object*2;
    else
      nCoefficient = object*2 + 1;
  }
  else
  {
    assert(type == INFERIOR_CONJUNCTION || type == SUPERIOR_CONJUNCTION);

    if (type == INFERIOR_CONJUNCTION)
      nCoefficient = object*2;
    else
      nCoefficient = object*2 + 1;
  }
  assert(nCoefficient >= 0 && nCoefficient < sizeof(g_PlanetaryPhenomenaCoefficient1)/sizeof(PlanetaryPhenomenaCoefficient1));

  return g_PlanetaryPhenomenaCoefficient1[nCoefficient].A + g_PlanetaryPhenomenaCoefficient1[nCoefficient].B*k;
}

double CAAPlanetaryPhenomena::True(double k, PlanetaryObject object, EventType type)
{
  double JDE0;
  
  if (type == WESTERN_ELONGATION || type == EASTERN_ELONGATION || type == STATION1 || type == STATION2)
  {
    if (object >= MARS)
      JDE0 = Mean(k, object, OPPOSITION);
    else
      JDE0 = Mean(k, object, INFERIOR_CONJUNCTION);
  }
  else
    JDE0 = Mean(k, object, type);

  size_t nCoefficient;
  if (object >= MARS)
  {
    assert(type == OPPOSITION || type == CONJUNCTION || type == STATION1 || type == STATION2);

    if (type == OPPOSITION || type == STATION1 || type == STATION2)
      nCoefficient = object*2;
    else
      nCoefficient = object*2 + 1;
  }
  else
  {
    assert(type == INFERIOR_CONJUNCTION || type == SUPERIOR_CONJUNCTION || type == EASTERN_ELONGATION || 
           type == WESTERN_ELONGATION || type == STATION1 || type == STATION2);

    if (type == INFERIOR_CONJUNCTION || type == EASTERN_ELONGATION || type == WESTERN_ELONGATION || type == STATION1 || type == STATION2)
      nCoefficient = object*2;
    else
      nCoefficient = object*2 + 1;
  }
  assert(nCoefficient >= 0 && nCoefficient < sizeof(g_PlanetaryPhenomenaCoefficient1)/sizeof(PlanetaryPhenomenaCoefficient1));

  double M = CAACoordinateTransformation::MapTo0To360Range(g_PlanetaryPhenomenaCoefficient1[nCoefficient].M0 + g_PlanetaryPhenomenaCoefficient1[nCoefficient].M1*k);
  M = CAACoordinateTransformation::DegreesToRadians(M);   //convert M to radians

  double T = (JDE0 - 2451545) / 36525;
  double T2 = T*T;

  double a=0;
  double b=0;
  double c=0;
  double d=0;
  double e=0;
  double f=0;
  double g=0;

  if (object == JUPITER)
  {
    a = CAACoordinateTransformation::MapTo0To360Range(82.74 + 40.76*T);
    a = CAACoordinateTransformation::DegreesToRadians(a);    
  }
  else if (object == SATURN)
  {
    a = CAACoordinateTransformation::MapTo0To360Range(82.74 + 40.76*T);
    a = CAACoordinateTransformation::DegreesToRadians(a);    
    b = CAACoordinateTransformation::MapTo0To360Range(29.86 + 1181.36*T);
    b = CAACoordinateTransformation::DegreesToRadians(b);    
    c = CAACoordinateTransformation::MapTo0To360Range(14.13 + 590.68*T);
    c = CAACoordinateTransformation::DegreesToRadians(c);    
    d = CAACoordinateTransformation::MapTo0To360Range(220.02 + 1262.87*T);
    d = CAACoordinateTransformation::DegreesToRadians(d);    
  }
  else if (object == URANUS)
  {
    e = CAACoordinateTransformation::MapTo0To360Range(207.83 + 8.51*T);
    e = CAACoordinateTransformation::DegreesToRadians(e);    
    f = CAACoordinateTransformation::MapTo0To360Range(108.84 + 419.96*T);
    f = CAACoordinateTransformation::DegreesToRadians(f);    
  }
  else if (object == NEPTUNE)
  {
    e = CAACoordinateTransformation::MapTo0To360Range(207.83 + 8.51*T);
    e = CAACoordinateTransformation::DegreesToRadians(e);    
    g = CAACoordinateTransformation::MapTo0To360Range(276.74 + 209.98*T);
    g = CAACoordinateTransformation::DegreesToRadians(g);    
  }

  double delta;
  if (object == MERCURY)
  {
    if (type == INFERIOR_CONJUNCTION)
    {
      delta = (0.0545 + 0.0002*T) + 
               sin(M) * (-6.2008 + 0.0074*T + 0.00003*T2) +
               cos(M) * (-3.2750 - 0.0197*T + 0.00001*T2) +
               sin(2*M) * (0.4737 - 0.0052*T - 0.00001*T2) +
               cos(2*M) * (0.8111 + 0.0033*T - 0.00002*T2) +
               sin(3*M) * (0.0037 + 0.0018*T) +
               cos(3*M) * (-0.1768 + 0.00001*T2) +
               sin(4*M) * (-0.0211 - 0.0004*T) +
               cos(4*M) * (0.0326 - 0.0003*T) +
               sin(5*M) * (0.0083 + 0.0001*T) +
               cos(5*M) * (-0.0040 + 0.0001*T);
    }
    else if (type == SUPERIOR_CONJUNCTION)
    {
      delta = (-0.0548 - 0.0002*T) +
              sin(M) * (7.3894 - 0.0100*T - 0.00003*T2) +
              cos(M) * (3.2200 + 0.0197*T - 0.00001*T2) +
              sin(2*M) * (0.8383 - 0.0064*T - 0.00001*T2) +
              cos(2*M) * (0.9666 + 0.0039*T - 0.00003*T2) +
              sin(3*M) * (0.0770 - 0.0026*T) +
              cos(3*M) * (0.2758 + 0.0002*T - 0.00002*T2) +
              sin(4*M) * (-0.0128 - 0.0008*T) +
              cos(4*M) * (0.0734 - 0.0004*T - 0.00001*T2) +
              sin(5*M) * (-0.0122 - 0.0002*T) +
              cos(5*M) * (0.0173 - 0.0002*T); 
    }  
    else if (type == EASTERN_ELONGATION)
    {
      delta = (-21.6101 + 0.0002*T) +
              sin(M) * (-1.9803 - 0.0060*T + 0.00001*T2) +
              cos(M) * (1.4151 - 0.0072*T - 0.00001*T2) +
              sin(2*M) * (0.5528 - 0.0005*T - 0.00001*T2) +
              cos(2*M) * (0.2905 + 0.0034*T + 0.00001*T2) +
              sin(3*M) * (-0.1121 - 0.0001*T + 0.00001*T2) +
              cos(3*M) * (-0.0098 - 0.0015*T) +
              sin(4*M) * (0.0192) +
              cos(4*M) * (0.0111 + 0.0004*T) +
              sin(5*M) * (-0.0061) +
              cos(5*M) * (-0.0032 - 0.0001*T2);
    }
    else if (type == WESTERN_ELONGATION)
    {
      delta = (21.6249 - 0.0002*T) +
              sin(M) * (0.1306 + 0.0065*T) +
              cos(M) * (-2.7661 - 0.0011*T + 0.00001*T2) +
              sin(2*M) * (0.2438 - 0.0024*T - 0.00001*T2) +
              cos(2*M) * (0.5767 + 0.0023*T) +
              sin(3*M) * (0.1041) +
              cos(3*M) * (-0.0184 + 0.0007*T) +
              sin(4*M) * (-0.0051 - 0.0001*T) +
              cos(4*M) * (0.0048 + 0.0001*T) +
              sin(5*M) * (0.0026) +
              cos(5*M) * (0.0037);
    }
    else if (type == STATION1)
    {
      delta = (-11.0761 + 0.0003*T) +
              sin(M) * (-4.7321 + 0.0023*T + 0.00002*T2) + 
              cos(M) * (-1.3230 - 0.0156*T) +
              sin(2*M) * (0.2270 - 0.0046*T) +
              cos(2*M) * (0.7184 + 0.0013*T - 0.00002*T2) +
              sin(3*M) * (0.0638 + 0.0016*T) +
              cos(3*M) * (-0.1655 + 0.0007*T) +
              sin(4*M) * (-0.0395 - 0.0003*T) + 
              cos(4*M) * (0.0247 - 0.0006*T) +
              sin(5*M) * (0.0131) +
              cos(5*M) * (0.0008 + 0.0002*T); 
    }
    else
    {
      assert(type == STATION2);

      delta = (11.1343 - 0.0001*T) +
              sin(M) * (-3.9137 + 0.0073*T + 0.00002*T2) + 
              cos(M) * (-3.3861 - 0.0128*T + 0.00001*T2) +
              sin(2*M) * (0.5222 - 0.0040*T - 0.00002*T2) +
              cos(2*M) * (0.5929 + 0.0039*T - 0.00002*T2) +
              sin(3*M) * (-0.0593 + 0.0018*T) +
              cos(3*M) * (-0.1733 - 0.0007*T + 0.00001*T2) +
              sin(4*M) * (-0.0053 - 0.0006*T) + 
              cos(4*M) * (0.0476 - 0.0001*T) +
              sin(5*M) * (0.0070 + 0.0002*T) +
              cos(5*M) * (-0.0115 + 0.0001*T); 
    }
  }
  else if (object == VENUS)
  {
    if (type == INFERIOR_CONJUNCTION)
    {
      delta = (-0.0096 + 0.0002*T - 0.00001*T2) +
              sin(M) * (2.0009 - 0.0033*T - 0.00001*T2) +
              cos(M) * (0.5980 - 0.0104*T + 0.00001*T2) +
              sin(2*M) * (0.0967 - 0.0018*T - 0.00003*T2) +
              cos(2*M) * (0.0913 + 0.0009*T - 0.00002*T2) +
              sin(3*M) * (0.0046 - 0.0002*T) +
              cos(3*M) * (0.0079 + 0.0001*T);
    }
    else if (type == SUPERIOR_CONJUNCTION)
    {
      delta = (0.0099 - 0.0002*T - 0.00001*T2) +
              sin(M) * (4.1991 - 0.0121*T - 0.00003*T2) +
              cos(M) * (-0.6095 + 0.0102*T - 0.00002*T2) +
              sin(2*M) * (0.2500 - 0.0028*T - 0.00003*T2) +
              cos(2*M) * (0.0063 + 0.0025*T - 0.00002*T2) +
              sin(3*M) * (0.0232 - 0.0005*T - 0.00001*T2) + 
              cos(3*M) * (0.0031 + 0.0004*T);
    }  
    else if (type == EASTERN_ELONGATION)
    {
      delta = (-70.7600 + 0.0002*T - 0.00001*T2) +
              sin(M) * (1.0282 - 0.0010*T - 0.00001*T2) +
              cos(M) * (0.2761 - 0.0060*T) +
              sin(2*M) * (-0.0438 - 0.0023*T + 0.00002*T2) +
              cos(2*M) * (0.1660 - 0.0037*T - 0.00004*T2) +
              sin(3*M) * (0.0036 + 0.0001*T) +
              cos(3*M) * (-0.0011 + 0.00001*T2);
    }
    else if (type == WESTERN_ELONGATION)
    {
      delta = (70.7462 - 0.00001*T2) +
              sin(M) * (1.1218 - 0.0025*T - 0.00001*T2) +
              cos(M) * (0.4538 - 0.0066*T) +
              sin(2*M) * (0.1320 + 0.0020*T - 0.00003*T2) +
              cos(2*M) * (-0.0702 + 0.0022*T + 0.00004*T2) +
              sin(3*M) * (0.0062 - 0.0001*T) +
              cos(3*M) * (0.0015 - 0.00001*T2);
    }
    else if (type == STATION1)
    {
      delta = (-21.0672 + 0.0002*T - 0.00001*T2) +
          sin(M) * (1.9396 - 0.0029*T - 0.00001*T2) + 
          cos(M) * (1.0727 - 0.0102*T) +
          sin(2*M) * (0.0404 - 0.0023*T - 0.00001*T2) +
          cos(2*M) * (0.1305 - 0.0004*T - 0.00003*T2) +
          sin(3*M) * (-0.0007 - 0.0002*T) +
          cos(3*M) * (0.0098);
    }
    else
    {
      assert(type == STATION2);

      delta = (21.0623 - 0.00001*T2) +
              sin(M) * (1.9913 - 0.0040*T - 0.00001*T2) + 
              cos(M) * (-0.0407 - 0.0077*T) +
              sin(2*M) * (0.1351 - 0.0009*T - 0.00004*T2) +
              cos(2*M) * (0.0303 + 0.0019*T) +
              sin(3*M) * (0.0089 - 0.0002*T) +
              cos(3*M) * (0.0043 + 0.0001*T);
    }
  }
  else if (object == MARS)
  {
    if (type == OPPOSITION)
    {
      delta = (-0.3088 + 0.00002*T2) +
              sin(M) * (-17.6965 + 0.0363*T + 0.00005*T2) +
              cos(M) * (18.3131 + 0.0467*T - 0.00006*T2) +
              sin(2*M) * (-0.2162 - 0.0198*T - 0.00001*T2) +
              cos(2*M) * (-4.5028 - 0.0019*T + 0.00007*T2) +
              sin(3*M) * (0.8987 + 0.0058*T - 0.00002*T2) +
              cos(3*M) * (0.7666 - 0.0050*T - 0.00003*T2) +
              sin(4*M) * (-0.3636 - 0.0001*T + 0.00002*T2) +
              cos(4*M) * (0.0402 + 0.0032*T) +
              sin(5*M) * (0.0737 - 0.0008*T) +
              cos(5*M) * (-0.0980 - 0.0011*T);
    }
    else if (type == CONJUNCTION)
    {
      delta = (0.3102 - 0.0001*T + 0.00001*T2) +
              sin(M) * (9.7273 - 0.0156*T + 0.00001*T2) +
              cos(M) * (-18.3195 - 0.0467*T + 0.00009*T2) +
              sin(2*M) * (-1.6488 - 0.0133*T + 0.00001*T2) +
              cos(2*M) * (-2.6117 - 0.0020*T + 0.00004*T2) +
              sin(3*M) * (-0.6827 - 0.0026*T + 0.00001*T2) +
              cos(3*M) * (0.0281 + 0.0035*T + 0.00001*T2) +
              sin(4*M) * (-0.0823 + 0.0006*T + 0.00001*T2) +
              cos(4*M) * (0.1584 + 0.0013*T) +
              sin(5*M) * (0.0270 + 0.0005*T) +
              cos(5*M) * (0.0433);
    }  
    else if (type == STATION1)
    {
      delta = (-37.0790 - 0.0009*T + 0.00002*T2) +
              sin(M) * (-20.0651 + 0.0228*T + 0.00004*T2) + 
              cos(M) * (14.5205 + 0.0504 - 0.00001*T2) +
              sin(2*M) * (1.1737 - 0.0169*T) +
              cos(2*M) * (-4.2550 - 0.0075*T + 0.00008*T2) +
              sin(3*M) * (0.4897 + 0.0074*T - 0.00001*T2) +
              cos(3*M) * (1.1151 - 0.0021*T - 0.00005*T2) +
              sin(4*M) * (-0.3636 - 0.0020*T + 0.00001*T2) + 
              cos(4*M) * (-0.1769 + 0.0028*T + 0.00002*T2) +
              sin(5*M) * (0.1437 - 0.0004*T) +
              cos(5*M) * (-0.0383 - 0.0016*T); 
    }
    else
    {
      assert(type == STATION2);

      delta = (36.7191 + 0.0016*T + 0.00003*T2) +
              sin(M) * (-12.6163 + 0.0417*T - 0.00001*T2) + 
              cos(M) * (20.1218 + 0.0379*T - 0.00006*T2) +
              sin(2*M) * (-1.6360 - 0.0190*T) +
              cos(2*M) * (-3.9657 + 0.0045*T + 0.00007*T2) +
              sin(3*M) * (1.1546 + 0.0029*T - 0.00003*T2) +
              cos(3*M) * (0.2888 - 0.0073*T - 0.00002*T2) +
              sin(4*M) * (-0.3128 + 0.0017*T + 0.00002*T2) + 
              cos(4*M) * (0.2513 + 0.0026*T - 0.00002*T2) +
              sin(5*M) * (-0.0021 - 0.0016*T) +
              cos(5*M) * (-0.1497 - 0.0006*T); 
    }
  }
  else if (object == JUPITER)
  {
    if (type == OPPOSITION)
    {
      delta = (-0.1029 - 0.00009*T2) +
              sin(M) * (-1.9658 - 0.0056*T + 0.00007*T2) +
              cos(M) * (6.1537 + 0.0210*T - 0.00006*T2) +
              sin(2*M) * (-0.2081 - 0.0013*T) +
              cos(2*M) * (-0.1116 - 0.0010*T) +
              sin(3*M) * (0.0074 + 0.0001*T) +
              cos(3*M) * (-0.0097 - 0.0001*T) +
              sin(a) * (0.0144*T - 0.00008*T2) +
              cos(a) * (0.3642 - 0.0019*T - 0.00029*T2);
    }
    else if (type == CONJUNCTION)
    {
      delta = (0.1027 + 0.0002*T - 0.00009*T2) +
              sin(M) * (-2.2637 + 0.0163*T - 0.00003*T2) +
              cos(M) * (-6.1540 - 0.0210*T + 0.00008*T2) +
              sin(2*M) * (-0.2021 - 0.0017*T + 0.00001*T2) +
              cos(2*M) * (0.1310 - 0.0008*T) +
              sin(3*M) * (0.0086) +
              cos(3*M) * (0.0087 + 0.0002*T) +
              sin(a) * (0.0144*T - 0.00008*T2) +
              cos(a) * (0.3642 - 0.0019*T - 0.00029*T2);
    }  
    else if (type == STATION1)
    {
      delta = (-60.3670 - 0.0001*T - 0.00009*T2) +
              sin(M) * (-2.3144 - 0.0124*T + 0.00007*T2) + 
              cos(M) * (6.7439 + 0.0166*T - 0.00006*T2) +
              sin(2*M) * (-0.2259 - 0.0010*T) +
              cos(2*M) * (-0.1497 - 0.0014*T) +
              sin(3*M) * (0.0105 + 0.0001*T) +
              cos(3*M) * (-0.0098) +
              sin(a) * (0.0144*T - 0.00008*T2) + 
              cos(a) * (0.3642 - 0.0019*T - 0.00029*T2);
    }
    else
    {
      assert(type == STATION2);

      delta = (60.3023 + 0.0002*T - 0.00009*T2) +
              sin(M) * (0.3506 - 0.0034*T + 0.00004*T2) + 
              cos(M) * (5.3635 + 0.0247*T - 0.00007*T2) +
              sin(2*M) * (-0.1872 - 0.0016*T) +
              cos(2*M) * (-0.0037 - 0.0005*T) +
              sin(3*M) * (0.0012 + 0.0001*T) +
              cos(3*M) * (-0.0096 - 0.0001*T) +
              sin(a) * (0.0144*T - 0.00008*T2) + 
              cos(a) * (0.3642 - 0.0019*T - 0.00029*T2);
    }
  }
  else if (object == SATURN)
  {
    if (type == OPPOSITION)
    {
      delta = (-0.0209 + 0.0006*T + 0.00023*T2) +
              sin(M) * (4.5795 - 0.0312*T - 0.00017*T2) +
              cos(M) * (1.1462 - 0.0351*T + 0.00011*T2) +
              sin(2*M) * (0.0985 - 0.0015*T) +
              cos(2*M) * (0.0733 - 0.0031*T + 0.00001*T2) +
              sin(3*M) * (0.0025 - 0.0001*T) +
              cos(3*M) * (0.0050 - 0.0002*T) +
              sin(a) * (-0.0337*T + 0.00018*T2) +
              cos(a) * (-0.8510 + 0.0044*T + 0.00068*T2) +
              sin(b) * (-0.0064*T + 0.00004*T2) +
              cos(b) * (0.2397 - 0.0012*T - 0.00008*T2) +
              sin(c) * (-0.0010*T) +
              cos(c) * (0.1245 + 0.0006*T) +
              sin(d) * (0.0024*T - 0.00003*T2) +
              cos(d) * (0.0477 - 0.0005*T - 0.00006*T2);
    }
    else if (type == CONJUNCTION)
    {
      delta = (0.0172 - 0.0006*T + 0.00023*T2) +
              sin(M) * (-8.5885 + 0.0411*T + 0.00020*T2) +
              cos(M) * (-1.1470 + 0.0352*T - 0.00011*T2) +
              sin(2*M) * (0.3331 - 0.0034*T - 0.00001*T2) +
              cos(2*M) * (0.1145 - 0.0045*T + 0.00002*T2) +
              sin(3*M) * (-0.0169 + 0.0002*T) +
              cos(3*M) * (-0.0109 + 0.0004*T) +
              sin(a) * (-0.0337*T + 0.00018*T2) +
              cos(a) * (-0.8510 + 0.0044*T + 0.00068*T2) +
              sin(b) * (-0.0064*T + 0.00004*T2) +
              cos(b) * (0.2397 - 0.0012*T - 0.00008*T2) +
              sin(c) * (-0.0010*T) +
              cos(c) * (0.1245 + 0.0006*T) +
              sin(d) * (0.0024*T - 0.00003*T2) +
              cos(d) * (0.0477 - 0.0005*T - 0.00006*T2);
    }  
    else if (type == STATION1)
    {
      delta = (-68.8840 + 0.0009*T + 0.00023*T2) +
              sin(M) * (5.5452 - 0.0279*T - 0.00020*T2) + 
              cos(M) * (3.0727 - 0.0430*T + 0.00007*T2) +
              sin(2*M) * (0.1101 - 0.0006*T - 0.00001*T2) +
              cos(2*M) * (0.1654 - 0.0043*T + 0.00001*T2) +
              sin(3*M) * (0.0010 + 0.0001*T) +
              cos(3*M) * (0.0095 - 0.0003*T) +
              sin(a) * (-0.0337*T + 0.00018*T2) + 
              cos(a) * (-0.8510 + 0.0044*T + 0.00068*T2) +
              sin(b) * (-0.0064*T + 0.00004*T2) + 
              cos(b) * (0.2397 - 0.0012*T - 0.00008*T2) +
              sin(c) * (-0.0010*T) + 
              cos(c) * (0.1245 + 0.0006*T) +
              sin(d) * (0.0024*T - 0.00003*T2) + 
              cos(d) * (0.0477 - 0.0005*T - 0.00006*T2);
    }
    else
    {
      assert(type == STATION2);

      delta = (68.8720 - 0.0007*T + 0.00023*T2) +
              sin(M) * (5.9399 - 0.0400*T - 0.00015*T2) + 
              cos(M) * (-0.7998 - 0.0266*T + 0.00014*T2) +
              sin(2*M) * (0.1738 - 0.0032*T) +
              cos(2*M) * (-0.0039 - 0.0024*T + 0.00001*T2) +
              sin(3*M) * (0.0073 - 0.0002*T) +
              cos(3*M) * (0.0020 - 0.0002*T) +
              sin(a) * (-0.0337*T + 0.00018*T2) + 
              cos(a) * (-0.8510 + 0.0044*T + 0.00068*T2) +
              sin(b) * (-0.0064*T + 0.00004*T2) + 
              cos(b) * (0.2397 - 0.0012*T - 0.00008*T2) +
              sin(c) * (-0.0010*T) + 
              cos(c) * (0.1245 + 0.0006*T) +
              sin(d) * (0.0024*T - 0.00003*T2) + 
              cos(d) * (0.0477 - 0.0005*T - 0.00006*T2);
    }
  }
  else if (object == URANUS)
  {
    if (type == OPPOSITION)
    {
      delta = (0.0844 - 0.0006*T) +
              sin(M) * (-0.1048 + 0.0246*T) +
              cos(M) * (-5.1221 + 0.0104*T + 0.00003*T2) +
              sin(2*M) * (-0.1428 - 0.0005*T) +
              cos(2*M) * (-0.0148 - 0.0013*T) + 
              cos(3*M) * 0.0055 +
              cos(e) * 0.8850 + 
              cos(f) * 0.2153;
    }
    else
    {
      assert(type == CONJUNCTION);

      delta = (-0.0859 + 0.0003*T) +
              sin(M) * (-3.8179 - 0.0148*T + 0.00003*T2) +
              cos(M) * (5.1228 - 0.0105*T - 0.00002*T2) +
              sin(2*M) * (-0.0803 + 0.0011*T) +
              cos(2*M) * (-0.1905 - 0.0006*T) +
              sin(3*M) * (0.0088 + 0.0001*T) +
              cos(e) * 0.8850 +
              cos(f) * 0.2153;
    }
  }
  else
  {
    assert(object == NEPTUNE);
    
    if (type == OPPOSITION)
    {
      delta = (-0.0140 + 0.00001*T2) +
              sin(M) * (-1.3486 + 0.0010*T + 0.00001*T2) +
              cos(M) * (0.8597 + 0.0037*T) +
              sin(2*M) * (-0.0082 - 0.0002*T + 0.00001*T2) +
              cos(2*M) * (0.0037 - 0.0003*T) +
              cos(e) * (-0.5964) +
              cos(g) * (0.0728);
    }
    else
    {
      assert(type == CONJUNCTION);

      delta = (0.0168) +
              sin(M) * (-2.5606 + 0.0088*T + 0.00002*T2) +
              cos(M) * (-0.8611 - 0.0037*T + 0.00002*T2) +
              sin(2*M) * (0.0118 - 0.0004*T + 0.00001*T2) +
              cos(2*M) * (0.0307 - 0.0003*T) +
              cos(e) * (-0.5964) +
              cos(g) * (0.0728);
    }  
  }
  
  return JDE0 + delta;
}

double CAAPlanetaryPhenomena::ElongationValue(double k, PlanetaryObject object, bool bEastern)
{
  double JDE0 = Mean(k, object, INFERIOR_CONJUNCTION);

  assert(object < MARS);

  size_t nCoefficient = object*2;
  assert(nCoefficient >= 0 && nCoefficient < sizeof(g_PlanetaryPhenomenaCoefficient1)/sizeof(PlanetaryPhenomenaCoefficient1));

  double M = CAACoordinateTransformation::MapTo0To360Range(g_PlanetaryPhenomenaCoefficient1[nCoefficient].M0 + g_PlanetaryPhenomenaCoefficient1[nCoefficient].M1*k);
  M = CAACoordinateTransformation::DegreesToRadians(M);   //convert M to radians

  double T = (JDE0 - 2451545) / 36525;
  double T2 = T*T;

  double value = 0;
  if (object == MERCURY)
  {
    if (bEastern)
    {
      value = (22.4697) + 
               sin(M) * (-4.2666 + 0.0054*T + 0.00002*T2) +
               cos(M) * (-1.8537 - 0.0137*T) +
               sin(2*M) * (0.3598 + 0.0008*T - 0.00001*T2) +
               cos(2*M) * (-0.0680 + 0.0026*T) +
               sin(3*M) * (-0.0524 - 0.0003*T) +
               cos(3*M) * (0.0052 - 0.0006*T) +
               sin(4*M) * (0.0107 + 0.0001*T) +
               cos(4*M) * (-0.0013 + 0.0001*T) +
               sin(5*M) * (-0.0021) +
               cos(5*M) * (0.0003);
    }
    else
    {
      value = (22.4143 - 0.0001*T) +
              sin(M) * (4.3651 - 0.0048*T - 0.00002*T2) +
              cos(M) * (2.3787 + 0.0121*T - 0.00001*T2) +
              sin(2*M) * (0.2674 + 0.0022*T) +
              cos(2*M) * (-0.3873 + 0.0008*T + 0.00001*T2) +
              sin(3*M) * (-0.0369 - 0.0001*T) +
              cos(3*M) * (0.0017 - 0.0001*T) +
              sin(4*M) * (0.0059) +
              cos(4*M) * (0.0061 + 0.0001*T) +
              sin(5*M) * (0.0007) +
              cos(5*M) * (-0.0011);
    }  
  }
  else if (object == VENUS)
  {
    if (bEastern)
    {
      value = (46.3173 + 0.0001*T) +
              sin(M) * (0.6916 - 0.0024*T) +
              cos(M) * (0.6676 - 0.0045*T) +
              sin(2*M) * (0.0309 - 0.0002*T) +
              cos(2*M) * (0.0036 - 0.0001*T);
    }
    else
    {
      value = (46.3245) +
              sin(M) * (-0.5366 - 0.0003*T + 0.00001*T2) +
              cos(M) * (0.3097 + 0.0016*T - 0.00001*T2) +
              sin(2*M) * (-0.0163) +
              cos(2*M) * (-0.0075 + 0.0001*T);
    }  
  }
  
  return value;
}

//###############################################################################
//###############################################################################
//###############################################################################
// 
////////////////////////////////////////////////////////////////////////////////
//##############################################################################

RcppExport SEXP MERCURYPlanet (SEXP year_){
  double Year = as <double> (year_); 
  
  //1 MERCURY STATION1
  double Kpp_STATION1          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::STATION1);
  double MERCURY_STATION1_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION1, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::STATION1);
  double MERCURY_STATION1_true = CAAPlanetaryPhenomena::True(Kpp_STATION1, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::STATION1);
  
  //2 MERCURY STATION2
  double Kpp_STATION2          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::STATION2);
  double MERCURY_STATION2_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION2, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::STATION2);
  double MERCURY_STATION2_true = CAAPlanetaryPhenomena::True(Kpp_STATION2, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::STATION2);
    
  //3 MERCURY OPPOSITION
  double Kpp_OPPOSITION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::OPPOSITION);
  double MERCURY_OPPOSITION_mean = CAAPlanetaryPhenomena::Mean(Kpp_OPPOSITION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::OPPOSITION);
  double MERCURY_OPPOSITION_true = CAAPlanetaryPhenomena::True(Kpp_OPPOSITION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::OPPOSITION);
  
  //4 MERCURY CONJUNCTION
  double Kpp_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::CONJUNCTION);
  double MERCURY_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::CONJUNCTION);
  double MERCURY_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::CONJUNCTION);
  
  //5 MERCURY EASTERN_ELONGATION
  double Kpp_EASTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double MERCURY_EASTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double MERCURY_EASTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double MERCURY_EASTERN_ELONGATION_Value =  CAAPlanetaryPhenomena::ElongationValue(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::MERCURY, true);
  
  //6 MERCURY WESTERN_ELONGATION
  double Kpp_WESTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double MERCURY_WESTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double MERCURY_WESTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double MERCURY_WESTERN_ELONGATION_Value =  CAAPlanetaryPhenomena::ElongationValue(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::MERCURY, false);
  
  //7 MERCURY INFERIOR_CONJUNCTION
  double Kpp_INFERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double MERCURY_INFERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double MERCURY_INFERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  
  
  //8 MERCURY_SUPERIOR_CONJUNCTION
  double Kpp_SUPERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double MERCURY_SUPERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double MERCURY_SUPERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::MERCURY, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  
  
  
  List res = List::create( Named("STATION1")            = List::create(Named("Kpp_STATION1")             = Kpp_STATION1,                     Named("STATION1_mean") = MERCURY_STATION1_mean,                                   Named("STATION1_true")= MERCURY_STATION1_true                                 ),
                           Named("STATION2")             = List::create(Named("Kpp_STATION2")             = Kpp_STATION2                     , Named("STATION2_mean") = MERCURY_STATION2_mean ,                          Named("STATION2_true") = MERCURY_STATION2_true                        ),
                           Named("OPPOSITION")           = List::create(Named("Kpp_OPPOSITION")           = Kpp_OPPOSITION                   , Named("OPPOSITION_mean") = MERCURY_OPPOSITION_mean,                       Named("OPPOSITION_true") = MERCURY_OPPOSITION_true                    ),
                           Named("CONJUNCTION")          = List::create(Named("Kpp_CONJUNCTION")          = Kpp_CONJUNCTION                  , Named("CONJUNCTION_mean") = MERCURY_CONJUNCTION_mean,                     Named("CONJUNCTION_true") = MERCURY_CONJUNCTION_true                  ),
                           Named("EASTERN_ELONGATION")   = List::create(Named("Kpp_EASTERN_ELONGATION")   = Kpp_EASTERN_ELONGATION           , Named("EASTERN_ELONGATION_mean") = MERCURY_EASTERN_ELONGATION_mean ,      Named("EASTERN_ELONGATION_true ") = MERCURY_EASTERN_ELONGATION_true ,   Named("MERCURY_EASTERN_ELONGATION_Value") = MERCURY_EASTERN_ELONGATION_Value  ),
                           Named("WESTERN_ELONGATION")   = List::create(Named("Kpp_WESTERN_ELONGATION")   = Kpp_WESTERN_ELONGATION           , Named("WESTERN_ELONGATION_mean") = MERCURY_WESTERN_ELONGATION_mean ,      Named("WESTERN_ELONGATION_true ") = MERCURY_WESTERN_ELONGATION_true  ,  Named("MERCURY_WESTERN_ELONGATION_Value") = MERCURY_WESTERN_ELONGATION_Value ),
                           Named("INFERIOR_CONJUNCTION") = List::create(Named("Kpp_INFERIOR_CONJUNCTION") = Kpp_INFERIOR_CONJUNCTION         , Named("INFERIOR_CONJUNCTION_mean") = MERCURY_INFERIOR_CONJUNCTION_mean,   Named("INFERIOR_CONJUNCTION_true") = MERCURY_INFERIOR_CONJUNCTION_true),
                           Named("SUPERIOR_CONJUNCTION") = List::create(Named("Kpp_SUPERIOR_CONJUNCTION") = Kpp_SUPERIOR_CONJUNCTION         , Named("SUPERIOR_CONJUNCTION_mean") = MERCURY_SUPERIOR_CONJUNCTION_mean,   Named("SUPERIOR_CONJUNCTION_true") = MERCURY_SUPERIOR_CONJUNCTION_true));
  return res;
}
  
  
///// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////// ################################################################################################################################
// Venus 
RcppExport SEXP VENUSPlanet (SEXP year_){
  double Year = as <double> (year_); 
  
  //1 VENUS STATION1
  double Kpp_STATION1          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::STATION1);
  double VENUS_STATION1_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION1, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::STATION1);
  double VENUS_STATION1_true = CAAPlanetaryPhenomena::True(Kpp_STATION1, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::STATION1);
  
  //2 VENUS STATION2
  double Kpp_STATION2          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::STATION2);
  double VENUS_STATION2_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION2, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::STATION2);
  double VENUS_STATION2_true = CAAPlanetaryPhenomena::True(Kpp_STATION2, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::STATION2);
    
  //3 VENUS OPPOSITION
  double Kpp_OPPOSITION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::OPPOSITION);
  double VENUS_OPPOSITION_mean = CAAPlanetaryPhenomena::Mean(Kpp_OPPOSITION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::OPPOSITION);
  double VENUS_OPPOSITION_true = CAAPlanetaryPhenomena::True(Kpp_OPPOSITION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::OPPOSITION);
  
  //4 VENUS CONJUNCTION
  double Kpp_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::CONJUNCTION);
  double VENUS_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::CONJUNCTION);
  double VENUS_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::CONJUNCTION);
  
  //5 VENUS EASTERN_ELONGATION
  double Kpp_EASTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double VENUS_EASTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double VENUS_EASTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
    double VENUS_EASTERN_ELONGATION_Value =  CAAPlanetaryPhenomena::ElongationValue(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::VENUS, true);
    
  //6 VENUS WESTERN_ELONGATION
  double Kpp_WESTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double VENUS_WESTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double VENUS_WESTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double VENUS_WESTERN_ELONGATION_Value =  CAAPlanetaryPhenomena::ElongationValue(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::VENUS, false);
    
  //7 VENUS INFERIOR_CONJUNCTION
  double Kpp_INFERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double VENUS_INFERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double VENUS_INFERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  
  //8 VENUS_SUPERIOR_CONJUNCTION
  double Kpp_SUPERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double VENUS_SUPERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double VENUS_SUPERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::VENUS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  
  List res = List::create( Named("STATION1")             = List::create(Named("Kpp_STATION1")             = Kpp_STATION1,                     Named("STATION1_mean") = VENUS_STATION1_mean,                           Named("STATION1_true")= VENUS_STATION1_true                                 ),
                           Named("STATION2")             = List::create(Named("Kpp_STATION2")             = Kpp_STATION2                     , Named("STATION2_mean") = VENUS_STATION2_mean ,                          Named("STATION2_true") = VENUS_STATION2_true                        ),
                           Named("OPPOSITION")           = List::create(Named("Kpp_OPPOSITION")           = Kpp_OPPOSITION                   , Named("OPPOSITION_mean") = VENUS_OPPOSITION_mean,                       Named("OPPOSITION_true") = VENUS_OPPOSITION_true                    ),
                           Named("CONJUNCTION")          = List::create(Named("Kpp_CONJUNCTION")          = Kpp_CONJUNCTION                  , Named("CONJUNCTION_mean") = VENUS_CONJUNCTION_mean,                     Named("CONJUNCTION_true") = VENUS_CONJUNCTION_true                  ),
                           Named("EASTERN_ELONGATION")   = List::create(Named("Kpp_EASTERN_ELONGATION")   = Kpp_EASTERN_ELONGATION           , Named("EASTERN_ELONGATION_mean") = VENUS_EASTERN_ELONGATION_mean ,      Named("EASTERN_ELONGATION_true ") = VENUS_EASTERN_ELONGATION_true   , Named("VENUS_EASTERN_ELONGATION_Value") = VENUS_EASTERN_ELONGATION_Value ),
                           Named("WESTERN_ELONGATION")   = List::create(Named("Kpp_WESTERN_ELONGATION")   = Kpp_WESTERN_ELONGATION           , Named("WESTERN_ELONGATION_mean") = VENUS_WESTERN_ELONGATION_mean ,      Named("WESTERN_ELONGATION_true ") = VENUS_WESTERN_ELONGATION_true   , Named("VENUS_WESTERN_ELONGATION_Value") = VENUS_WESTERN_ELONGATION_Value ),
                           Named("INFERIOR_CONJUNCTION") = List::create(Named("Kpp_INFERIOR_CONJUNCTION") = Kpp_INFERIOR_CONJUNCTION         , Named("INFERIOR_CONJUNCTION_mean") = VENUS_INFERIOR_CONJUNCTION_mean,   Named("INFERIOR_CONJUNCTION_true") = VENUS_INFERIOR_CONJUNCTION_true),
                           Named("SUPERIOR_CONJUNCTION") = List::create(Named("Kpp_SUPERIOR_CONJUNCTION") = Kpp_SUPERIOR_CONJUNCTION         , Named("SUPERIOR_CONJUNCTION_mean") = VENUS_SUPERIOR_CONJUNCTION_mean,   Named("SUPERIOR_CONJUNCTION_true") = VENUS_SUPERIOR_CONJUNCTION_true));
  return res;
}
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
   
///// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////// ################################################################################################################################
// Mars 
RcppExport SEXP MARSPlanet (SEXP year_){
  double Year = as <double> (year_); 
  
  //1 MARS STATION1
  double Kpp_STATION1          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::STATION1);
  double MARS_STATION1_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION1, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::STATION1);
  double MARS_STATION1_true = CAAPlanetaryPhenomena::True(Kpp_STATION1, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::STATION1);
  
  //2 MARS STATION2
  double Kpp_STATION2          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::STATION2);
  double MARS_STATION2_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION2, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::STATION2);
  double MARS_STATION2_true = CAAPlanetaryPhenomena::True(Kpp_STATION2, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::STATION2);
    
  //3 MARS OPPOSITION
  double Kpp_OPPOSITION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::OPPOSITION);
  double MARS_OPPOSITION_mean = CAAPlanetaryPhenomena::Mean(Kpp_OPPOSITION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::OPPOSITION);
  double MARS_OPPOSITION_true = CAAPlanetaryPhenomena::True(Kpp_OPPOSITION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::OPPOSITION);
  
  //4 MARS CONJUNCTION
  double Kpp_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::CONJUNCTION);
  double MARS_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::CONJUNCTION);
  double MARS_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::CONJUNCTION);
  
  //5 MARS EASTERN_ELONGATION
  double Kpp_EASTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double MARS_EASTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double MARS_EASTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  
  //6 MARS WESTERN_ELONGATION
  double Kpp_WESTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double MARS_WESTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double MARS_WESTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  
  //7 MARS INFERIOR_CONJUNCTION
  double Kpp_INFERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double MARS_INFERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double MARS_INFERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  
  //8 MARS_SUPERIOR_CONJUNCTION
  double Kpp_SUPERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double MARS_SUPERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double MARS_SUPERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::MARS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  
  List res = List::create( Named("STATION1")             = List::create(Named("Kpp_STATION1")             = Kpp_STATION1,                     Named("STATION1_mean") = MARS_STATION1_mean,                            Named("STATION1_true")= MARS_STATION1_true                                 ),
                           Named("STATION2")             = List::create(Named("Kpp_STATION2")             = Kpp_STATION2                     , Named("STATION2_mean") = MARS_STATION2_mean ,                          Named("STATION2_true") = MARS_STATION2_true                        ),
                           Named("OPPOSITION")           = List::create(Named("Kpp_OPPOSITION")           = Kpp_OPPOSITION                   , Named("OPPOSITION_mean") = MARS_OPPOSITION_mean,                       Named("OPPOSITION_true") = MARS_OPPOSITION_true                    ),
                           Named("CONJUNCTION")          = List::create(Named("Kpp_CONJUNCTION")          = Kpp_CONJUNCTION                  , Named("CONJUNCTION_mean") = MARS_CONJUNCTION_mean,                     Named("CONJUNCTION_true") = MARS_CONJUNCTION_true                  ),
                           Named("EASTERN_ELONGATION")   = List::create(Named("Kpp_EASTERN_ELONGATION")   = Kpp_EASTERN_ELONGATION           , Named("EASTERN_ELONGATION_mean") = MARS_EASTERN_ELONGATION_mean ,      Named("EASTERN_ELONGATION_true ") = MARS_EASTERN_ELONGATION_true   ),
                           Named("WESTERN_ELONGATION")   = List::create(Named("Kpp_WESTERN_ELONGATION")   = Kpp_WESTERN_ELONGATION           , Named("WESTERN_ELONGATION_mean") = MARS_WESTERN_ELONGATION_mean ,      Named("WESTERN_ELONGATION_true ") = MARS_WESTERN_ELONGATION_true   ),
                           Named("INFERIOR_CONJUNCTION") = List::create(Named("Kpp_INFERIOR_CONJUNCTION") = Kpp_INFERIOR_CONJUNCTION         , Named("INFERIOR_CONJUNCTION_mean") = MARS_INFERIOR_CONJUNCTION_mean,   Named("INFERIOR_CONJUNCTION_true") = MARS_INFERIOR_CONJUNCTION_true),
                           Named("SUPERIOR_CONJUNCTION") = List::create(Named("Kpp_SUPERIOR_CONJUNCTION") = Kpp_SUPERIOR_CONJUNCTION         , Named("SUPERIOR_CONJUNCTION_mean") = MARS_SUPERIOR_CONJUNCTION_mean,   Named("SUPERIOR_CONJUNCTION_true") = MARS_SUPERIOR_CONJUNCTION_true));
  return res;
}
  
  ///// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////// ################################################################################################################################
// Jupiter

RcppExport SEXP JUPITERPlanet (SEXP year_){
  double Year = as <double> (year_); 
  
  //1 JUPITER STATION1
  double Kpp_STATION1          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::STATION1);
  double JUPITER_STATION1_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION1, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::STATION1);
  double JUPITER_STATION1_true = CAAPlanetaryPhenomena::True(Kpp_STATION1, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::STATION1);
  
  //2 JUPITER STATION2
  double Kpp_STATION2          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::STATION2);
  double JUPITER_STATION2_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION2, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::STATION2);
  double JUPITER_STATION2_true = CAAPlanetaryPhenomena::True(Kpp_STATION2, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::STATION2);
    
  //3 JUPITER OPPOSITION
  double Kpp_OPPOSITION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::OPPOSITION);
  double JUPITER_OPPOSITION_mean = CAAPlanetaryPhenomena::Mean(Kpp_OPPOSITION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::OPPOSITION);
  double JUPITER_OPPOSITION_true = CAAPlanetaryPhenomena::True(Kpp_OPPOSITION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::OPPOSITION);
  
  //4 JUPITER CONJUNCTION
  double Kpp_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::CONJUNCTION);
  double JUPITER_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::CONJUNCTION);
  double JUPITER_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::CONJUNCTION);
  
  //5 JUPITER EASTERN_ELONGATION
  double Kpp_EASTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double JUPITER_EASTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double JUPITER_EASTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  
  //6 JUPITER WESTERN_ELONGATION
  double Kpp_WESTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double JUPITER_WESTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double JUPITER_WESTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  
  //7 JUPITER INFERIOR_CONJUNCTION
  double Kpp_INFERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double JUPITER_INFERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double JUPITER_INFERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  
  //8 JUPITER_SUPERIOR_CONJUNCTION
  double Kpp_SUPERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double JUPITER_SUPERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double JUPITER_SUPERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::JUPITER, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  
  List res = List::create( Named("STATION1")             = List::create(Named("Kpp_STATION1")             = Kpp_STATION1,                     Named("STATION1_mean") = JUPITER_STATION1_mean,                            Named("STATION1_true")= JUPITER_STATION1_true                                 ),
                           Named("STATION2")             = List::create(Named("Kpp_STATION2")             = Kpp_STATION2                     , Named("STATION2_mean") = JUPITER_STATION2_mean ,                          Named("STATION2_true") = JUPITER_STATION2_true                        ),
                           Named("OPPOSITION")           = List::create(Named("Kpp_OPPOSITION")           = Kpp_OPPOSITION                   , Named("OPPOSITION_mean") = JUPITER_OPPOSITION_mean,                       Named("OPPOSITION_true") = JUPITER_OPPOSITION_true                    ),
                           Named("CONJUNCTION")          = List::create(Named("Kpp_CONJUNCTION")          = Kpp_CONJUNCTION                  , Named("CONJUNCTION_mean") = JUPITER_CONJUNCTION_mean,                     Named("CONJUNCTION_true") = JUPITER_CONJUNCTION_true                  ),
                           Named("EASTERN_ELONGATION")   = List::create(Named("Kpp_EASTERN_ELONGATION")   = Kpp_EASTERN_ELONGATION           , Named("EASTERN_ELONGATION_mean") = JUPITER_EASTERN_ELONGATION_mean ,      Named("EASTERN_ELONGATION_true ") = JUPITER_EASTERN_ELONGATION_true   ),
                           Named("WESTERN_ELONGATION")   = List::create(Named("Kpp_WESTERN_ELONGATION")   = Kpp_WESTERN_ELONGATION           , Named("WESTERN_ELONGATION_mean") = JUPITER_WESTERN_ELONGATION_mean ,      Named("WESTERN_ELONGATION_true ") = JUPITER_WESTERN_ELONGATION_true   ),
                           Named("INFERIOR_CONJUNCTION") = List::create(Named("Kpp_INFERIOR_CONJUNCTION") = Kpp_INFERIOR_CONJUNCTION         , Named("INFERIOR_CONJUNCTION_mean") = JUPITER_INFERIOR_CONJUNCTION_mean,   Named("INFERIOR_CONJUNCTION_true") = JUPITER_INFERIOR_CONJUNCTION_true),
                           Named("SUPERIOR_CONJUNCTION") = List::create(Named("Kpp_SUPERIOR_CONJUNCTION") = Kpp_SUPERIOR_CONJUNCTION         , Named("SUPERIOR_CONJUNCTION_mean") = JUPITER_SUPERIOR_CONJUNCTION_mean,   Named("SUPERIOR_CONJUNCTION_true") = JUPITER_SUPERIOR_CONJUNCTION_true));
  return res;
}
  

  ///// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////// ################################################################################################################################
// Saturn

RcppExport SEXP SATURNPlanet (SEXP year_){
  double Year = as <double> (year_); 
  
  //1 SATURN STATION1
  double Kpp_STATION1          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::STATION1);
  double SATURN_STATION1_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION1, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::STATION1);
  double SATURN_STATION1_true = CAAPlanetaryPhenomena::True(Kpp_STATION1, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::STATION1);
  
  //2 SATURN STATION2
  double Kpp_STATION2          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::STATION2);
  double SATURN_STATION2_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION2, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::STATION2);
  double SATURN_STATION2_true = CAAPlanetaryPhenomena::True(Kpp_STATION2, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::STATION2);
    
  //3 SATURN OPPOSITION
  double Kpp_OPPOSITION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::OPPOSITION);
  double SATURN_OPPOSITION_mean = CAAPlanetaryPhenomena::Mean(Kpp_OPPOSITION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::OPPOSITION);
  double SATURN_OPPOSITION_true = CAAPlanetaryPhenomena::True(Kpp_OPPOSITION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::OPPOSITION);
  
  //4 SATURN CONJUNCTION
  double Kpp_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::CONJUNCTION);
  double SATURN_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::CONJUNCTION);
  double SATURN_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::CONJUNCTION);
  
  //5 SATURN EASTERN_ELONGATION
  double Kpp_EASTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double SATURN_EASTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double SATURN_EASTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  
  //6 SATURN WESTERN_ELONGATION
  double Kpp_WESTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double SATURN_WESTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double SATURN_WESTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  
  //7 SATURN INFERIOR_CONJUNCTION
  double Kpp_INFERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double SATURN_INFERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double SATURN_INFERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  
  //8 SATURN_SUPERIOR_CONJUNCTION
  double Kpp_SUPERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double SATURN_SUPERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double SATURN_SUPERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::SATURN, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  
  List res = List::create( Named("STATION1")             = List::create(Named("Kpp_STATION1")             = Kpp_STATION1,                     Named("STATION1_mean") = SATURN_STATION1_mean,                           Named("STATION1_true")= SATURN_STATION1_true                                 ),
                           Named("STATION2")             = List::create(Named("Kpp_STATION2")             = Kpp_STATION2                     , Named("STATION2_mean") = SATURN_STATION2_mean ,                          Named("STATION2_true") = SATURN_STATION2_true                        ),
                           Named("OPPOSITION")           = List::create(Named("Kpp_OPPOSITION")           = Kpp_OPPOSITION                   , Named("OPPOSITION_mean") = SATURN_OPPOSITION_mean,                       Named("OPPOSITION_true") = SATURN_OPPOSITION_true                    ),
                           Named("CONJUNCTION")          = List::create(Named("Kpp_CONJUNCTION")          = Kpp_CONJUNCTION                  , Named("CONJUNCTION_mean") = SATURN_CONJUNCTION_mean,                     Named("CONJUNCTION_true") = SATURN_CONJUNCTION_true                  ),
                           Named("EASTERN_ELONGATION")   = List::create(Named("Kpp_EASTERN_ELONGATION")   = Kpp_EASTERN_ELONGATION           , Named("EASTERN_ELONGATION_mean") = SATURN_EASTERN_ELONGATION_mean ,      Named("EASTERN_ELONGATION_true ") = SATURN_EASTERN_ELONGATION_true   ),
                           Named("WESTERN_ELONGATION")   = List::create(Named("Kpp_WESTERN_ELONGATION")   = Kpp_WESTERN_ELONGATION           , Named("WESTERN_ELONGATION_mean") = SATURN_WESTERN_ELONGATION_mean ,      Named("WESTERN_ELONGATION_true ") = SATURN_WESTERN_ELONGATION_true   ),
                           Named("INFERIOR_CONJUNCTION") = List::create(Named("Kpp_INFERIOR_CONJUNCTION") = Kpp_INFERIOR_CONJUNCTION         , Named("INFERIOR_CONJUNCTION_mean") = SATURN_INFERIOR_CONJUNCTION_mean,   Named("INFERIOR_CONJUNCTION_true") = SATURN_INFERIOR_CONJUNCTION_true),
                           Named("SUPERIOR_CONJUNCTION") = List::create(Named("Kpp_SUPERIOR_CONJUNCTION") = Kpp_SUPERIOR_CONJUNCTION         , Named("SUPERIOR_CONJUNCTION_mean") = SATURN_SUPERIOR_CONJUNCTION_mean,   Named("SUPERIOR_CONJUNCTION_true") = SATURN_SUPERIOR_CONJUNCTION_true));
  return res;
}
  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////// ################################################################################################################################
// Uranus 

RcppExport SEXP URANUSPlanet (SEXP year_){
  double Year = as <double> (year_); 
  
  //1 URANUS STATION1
  double Kpp_STATION1          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::STATION1);
  double URANUS_STATION1_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION1, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::STATION1);
  double URANUS_STATION1_true = CAAPlanetaryPhenomena::True(Kpp_STATION1, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::STATION1);
  
  //2 URANUS STATION2
  double Kpp_STATION2          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::STATION2);
  double URANUS_STATION2_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION2, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::STATION2);
  double URANUS_STATION2_true = CAAPlanetaryPhenomena::True(Kpp_STATION2, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::STATION2);
    
  //3 URANUS OPPOSITION
  double Kpp_OPPOSITION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::OPPOSITION);
  double URANUS_OPPOSITION_mean = CAAPlanetaryPhenomena::Mean(Kpp_OPPOSITION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::OPPOSITION);
  double URANUS_OPPOSITION_true = CAAPlanetaryPhenomena::True(Kpp_OPPOSITION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::OPPOSITION);
  
  //4 URANUS CONJUNCTION
  double Kpp_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::CONJUNCTION);
  double URANUS_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::CONJUNCTION);
  double URANUS_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::CONJUNCTION);
  
  //5 URANUS EASTERN_ELONGATION
  double Kpp_EASTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double URANUS_EASTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double URANUS_EASTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  
  //6 URANUS WESTERN_ELONGATION
  double Kpp_WESTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double URANUS_WESTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double URANUS_WESTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  
  //7 URANUS INFERIOR_CONJUNCTION
  double Kpp_INFERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double URANUS_INFERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double URANUS_INFERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  
  //8 URANUS_SUPERIOR_CONJUNCTION
  double Kpp_SUPERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double URANUS_SUPERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double URANUS_SUPERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::URANUS, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  
  List res = List::create( Named("STATION1")             = List::create(Named("Kpp_STATION1")             = Kpp_STATION1,                     Named("STATION1_mean") = URANUS_STATION1_mean,                           Named("STATION1_true")= URANUS_STATION1_true                                 ),
                           Named("STATION2")             = List::create(Named("Kpp_STATION2")             = Kpp_STATION2                     , Named("STATION2_mean") = URANUS_STATION2_mean ,                          Named("STATION2_true") = URANUS_STATION2_true                        ),
                           Named("OPPOSITION")           = List::create(Named("Kpp_OPPOSITION")           = Kpp_OPPOSITION                   , Named("OPPOSITION_mean") = URANUS_OPPOSITION_mean,                       Named("OPPOSITION_true") = URANUS_OPPOSITION_true                    ),
                           Named("CONJUNCTION")          = List::create(Named("Kpp_CONJUNCTION")          = Kpp_CONJUNCTION                  , Named("CONJUNCTION_mean") = URANUS_CONJUNCTION_mean,                     Named("CONJUNCTION_true") = URANUS_CONJUNCTION_true                  ),
                           Named("EASTERN_ELONGATION")   = List::create(Named("Kpp_EASTERN_ELONGATION")   = Kpp_EASTERN_ELONGATION           , Named("EASTERN_ELONGATION_mean") = URANUS_EASTERN_ELONGATION_mean ,      Named("EASTERN_ELONGATION_true ") = URANUS_EASTERN_ELONGATION_true   ),
                           Named("WESTERN_ELONGATION")   = List::create(Named("Kpp_WESTERN_ELONGATION")   = Kpp_WESTERN_ELONGATION           , Named("WESTERN_ELONGATION_mean") = URANUS_WESTERN_ELONGATION_mean ,      Named("WESTERN_ELONGATION_true ") = URANUS_WESTERN_ELONGATION_true   ),
                           Named("INFERIOR_CONJUNCTION") = List::create(Named("Kpp_INFERIOR_CONJUNCTION") = Kpp_INFERIOR_CONJUNCTION         , Named("INFERIOR_CONJUNCTION_mean") = URANUS_INFERIOR_CONJUNCTION_mean,   Named("INFERIOR_CONJUNCTION_true") = URANUS_INFERIOR_CONJUNCTION_true),
                           Named("SUPERIOR_CONJUNCTION") = List::create(Named("Kpp_SUPERIOR_CONJUNCTION") = Kpp_SUPERIOR_CONJUNCTION         , Named("SUPERIOR_CONJUNCTION_mean") = URANUS_SUPERIOR_CONJUNCTION_mean,   Named("SUPERIOR_CONJUNCTION_true") = URANUS_SUPERIOR_CONJUNCTION_true));
  return res;
}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////################################################################################################################################
// Neptune

RcppExport SEXP NEPTUNEPlanet (SEXP year_){
  double Year = as <double> (year_); 
  
  //1 NEPTUNE STATION1
  double Kpp_STATION1          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::STATION1);
  double NEPTUNE_STATION1_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION1, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::STATION1);
  double NEPTUNE_STATION1_true = CAAPlanetaryPhenomena::True(Kpp_STATION1, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::STATION1);
  
  //2 NEPTUNE STATION2
  double Kpp_STATION2          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::STATION2);
  double NEPTUNE_STATION2_mean = CAAPlanetaryPhenomena::Mean(Kpp_STATION2, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::STATION2);
  double NEPTUNE_STATION2_true = CAAPlanetaryPhenomena::True(Kpp_STATION2, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::STATION2);
    
  //3 NEPTUNE OPPOSITION
  double Kpp_OPPOSITION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::OPPOSITION);
  double NEPTUNE_OPPOSITION_mean = CAAPlanetaryPhenomena::Mean(Kpp_OPPOSITION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::OPPOSITION);
  double NEPTUNE_OPPOSITION_true = CAAPlanetaryPhenomena::True(Kpp_OPPOSITION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::OPPOSITION);
  
  //4 NEPTUNE CONJUNCTION
  double Kpp_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::CONJUNCTION);
  double NEPTUNE_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::CONJUNCTION);
  double NEPTUNE_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_CONJUNCTION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::CONJUNCTION);
  
  //5 NEPTUNE EASTERN_ELONGATION
  double Kpp_EASTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double NEPTUNE_EASTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  double NEPTUNE_EASTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_EASTERN_ELONGATION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::EASTERN_ELONGATION);
  
  //6 NEPTUNE WESTERN_ELONGATION
  double Kpp_WESTERN_ELONGATION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double NEPTUNE_WESTERN_ELONGATION_mean = CAAPlanetaryPhenomena::Mean(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  double NEPTUNE_WESTERN_ELONGATION_true = CAAPlanetaryPhenomena::True(Kpp_WESTERN_ELONGATION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::WESTERN_ELONGATION);
  
  //7 NEPTUNE INFERIOR_CONJUNCTION
  double Kpp_INFERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double NEPTUNE_INFERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  double NEPTUNE_INFERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_INFERIOR_CONJUNCTION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::INFERIOR_CONJUNCTION);
  
  //8 NEPTUNE_SUPERIOR_CONJUNCTION
  double Kpp_SUPERIOR_CONJUNCTION          = CAAPlanetaryPhenomena::K(Year, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double NEPTUNE_SUPERIOR_CONJUNCTION_mean = CAAPlanetaryPhenomena::Mean(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  double NEPTUNE_SUPERIOR_CONJUNCTION_true = CAAPlanetaryPhenomena::True(Kpp_SUPERIOR_CONJUNCTION, CAAPlanetaryPhenomena::NEPTUNE, CAAPlanetaryPhenomena::SUPERIOR_CONJUNCTION);
  
  List res = List::create( Named("STATION1")             = List::create(Named("Kpp_STATION1")             = Kpp_STATION1,                     Named("STATION1_mean") = NEPTUNE_STATION1_mean,                            Named("STATION1_true")= NEPTUNE_STATION1_true                         ),
                           Named("STATION2")             = List::create(Named("Kpp_STATION2")             = Kpp_STATION2                     , Named("STATION2_mean") = NEPTUNE_STATION2_mean ,                          Named("STATION2_true") = NEPTUNE_STATION2_true                        ),
                           Named("OPPOSITION")           = List::create(Named("Kpp_OPPOSITION")           = Kpp_OPPOSITION                   , Named("OPPOSITION_mean") = NEPTUNE_OPPOSITION_mean,                       Named("OPPOSITION_true") = NEPTUNE_OPPOSITION_true                    ),
                           Named("CONJUNCTION")          = List::create(Named("Kpp_CONJUNCTION")          = Kpp_CONJUNCTION                  , Named("CONJUNCTION_mean") = NEPTUNE_CONJUNCTION_mean,                     Named("CONJUNCTION_true") = NEPTUNE_CONJUNCTION_true                  ),
                           Named("EASTERN_ELONGATION")   = List::create(Named("Kpp_EASTERN_ELONGATION")   = Kpp_EASTERN_ELONGATION           , Named("EASTERN_ELONGATION_mean") = NEPTUNE_EASTERN_ELONGATION_mean ,      Named("EASTERN_ELONGATION_true ") = NEPTUNE_EASTERN_ELONGATION_true   ),
                           Named("WESTERN_ELONGATION")   = List::create(Named("Kpp_WESTERN_ELONGATION")   = Kpp_WESTERN_ELONGATION           , Named("WESTERN_ELONGATION_mean") = NEPTUNE_WESTERN_ELONGATION_mean ,      Named("WESTERN_ELONGATION_true ") = NEPTUNE_WESTERN_ELONGATION_true   ),
                           Named("INFERIOR_CONJUNCTION") = List::create(Named("Kpp_INFERIOR_CONJUNCTION") = Kpp_INFERIOR_CONJUNCTION         , Named("INFERIOR_CONJUNCTION_mean") = NEPTUNE_INFERIOR_CONJUNCTION_mean,   Named("INFERIOR_CONJUNCTION_true") = NEPTUNE_INFERIOR_CONJUNCTION_true),
                           Named("SUPERIOR_CONJUNCTION") = List::create(Named("Kpp_SUPERIOR_CONJUNCTION") = Kpp_SUPERIOR_CONJUNCTION         , Named("SUPERIOR_CONJUNCTION_mean") = NEPTUNE_SUPERIOR_CONJUNCTION_mean,   Named("SUPERIOR_CONJUNCTION_true") = NEPTUNE_SUPERIOR_CONJUNCTION_true));
  return res;
}
  
  
