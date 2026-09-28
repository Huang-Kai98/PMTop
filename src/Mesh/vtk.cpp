// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "vtk.h"
#include <cmath>
#include <iostream>
#include <ostream>

bool VTKGeometry::IsLagrange(int vtk_geom)
{
    return vtk_geom >= LAGRANGE_SEGMENT && vtk_geom <= LAGRANGE_PYRAMID;
}

bool VTKGeometry::IsQuadratic(int vtk_geom)
{
    return vtk_geom >= QUADRATIC_SEGMENT
           && vtk_geom <= BIQUADRATIC_QUADRATIC_PRISM;
}

int VTKGeometry::GetOrder(int vtk_geom, int npoints)
{
   if (IsQuadratic(vtk_geom))
   {
      return 2;
   }
   else if (IsLagrange(vtk_geom))
   {
      switch (vtk_geom)
      {
         case LAGRANGE_SEGMENT:
            return npoints - 1;
         case LAGRANGE_TRIANGLE:
            return static_cast<int>(std::sqrt(8*npoints + 1) - 3)/2;
         case LAGRANGE_SQUARE:
            return static_cast<int>(std::round(std::sqrt(npoints))) - 1;
         case LAGRANGE_TETRAHEDRON:
            switch (npoints)
            {
               // Note that for given order, npoints is given by
               // npoints_order = (order + 1)*(order + 2)*(order + 3)/6,
               case 4: return 1;
               case 10: return 2;
               case 20: return 3;
               case 35: return 4;
               case 56: return 5;
               case 84: return 6;
               case 120: return 7;
               case 165: return 8;
               case 220: return 9;
               case 286: return 10;
               default:
               {
                  constexpr int max_order = 20;
                  int order = 11, npoints_order;
                  for (; order<max_order; ++order)
                  {
                     npoints_order = (order + 1)*(order + 2)*(order + 3)/6;
                     if (npoints_order == npoints) { break; }
                  }
                  if (npoints_order != npoints) {
                     std::cerr<<"VTKGeometry::GetOrder"<<std::endl;
                  }
                  return order;
               }
            }
         case LAGRANGE_CUBE:
            return static_cast<int>(std::round(std::cbrt(npoints))) - 1;
         case LAGRANGE_PRISM:
         {
            const double n = npoints;
            static const double third = 1.0/3.0;
            static const double ninth = 1.0/9.0;
            static const double twentyseventh = 1.0/27.0;
            const double term =
               std::cbrt(third*sqrt(third)*sqrt((27.0*n - 2.0)*n) + n
                         - twentyseventh);
            return static_cast<int>(std::round(term + ninth / term - 4*third));
         }
         case LAGRANGE_PYRAMID:
            std::cerr<<"Lagrange pyramids not currently supported in VTK."<<std::endl;
            return 0;
      }
   }
   return 1;
}

int AnasysToVTKFormat(std::string anasys_format) {
   if (anasys_format == "Shell181" || anasys_format== "SHELL181"
      ||anasys_format == "Shell63" || anasys_format == "SHELL63") {
      return VTKGeometry::SQUARE;
      }
   else if (anasys_format == "Beam188" || anasys_format == "BEAM188"
      || anasys_format == "COMBIN14" || anasys_format == "Combin14") {
      return  VTKGeometry::SEGMENT;
   }
   else if (anasys_format =="Solid185"|| anasys_format == "SOLID185") {
    return VTKGeometry::CUBE;
   }
   else {
      std::cerr << "AnasysToVTKFormat: Unknown format " << anasys_format;
   }
   return -1;
}

