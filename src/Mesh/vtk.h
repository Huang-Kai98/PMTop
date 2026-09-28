/*=========================================================================

Program:   Visualization Toolkit
  Module:    vtkCellType.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkCellType
 * @brief   define types of cells
 *
 * vtkCellType defines the allowable cell types in the visualization
 * library (vtk). In vtk, datasets consist of collections of cells.
 * Different datasets consist of different cell types. The cells may be
 * explicitly represented (as in vtkPolyData), or may be implicit to the
 * data type (as in vtkStructuredPoints).
*/

// To add a new cell type, define a new integer type flag here, then
// create a subclass of vtkCell to implement the proper behavior. You
// may have to modify the following methods: vtkDataSet (and subclasses)
// GetCell() and vtkGenericCell::SetCellType(). Also, to do the job right,
// you'll also have to modify some filters (vtkGeometryFilter...) and
// regression tests (example scripts) to reflect the new cell addition.
// Also, make sure to update vtkCellTypesStrings in vtkCellTypes.cxx
// and the vtkCellTypes::IsLinear method in vtkCellTypes.h.

// .SECTION Caveats
// An unstructured grid stores the types of its cells as a
// unsigned char array. Therefore, the maximum encoding number for a cell type
// is 255.

#ifndef VTK_H
#define VTK_H

#include <cstdint>
#include <string>


// Helpers for reading and writing VTK format

/// @brief Helper class for converting between MFEM and VTK geometry types.
///
/// Note: The VTK element types defined are at: https://git.io/JvZLm
struct VTKGeometry
{
   /// @name VTK geometry types
   ///@{
   static const int POINT = 1;

   /// @name Low-order (linear, straight-sided) VTK geometric types
   ///@{
   static const int SEGMENT = 3;
   static const int TRIANGLE = 5;
   static const int SQUARE = 9;
   static const int TETRAHEDRON = 10;
   static const int CUBE = 12;
   static const int PRISM = 13;
   static const int PYRAMID = 14;
   ///@}

   /// @name Legacy quadratic VTK geometric types
   ///@{
   static const int QUADRATIC_SEGMENT = 21;
   static const int QUADRATIC_TRIANGLE = 22;
   static const int BIQUADRATIC_SQUARE = 28;
   static const int QUADRATIC_TETRAHEDRON = 24;
   static const int TRIQUADRATIC_CUBE = 29;
   static const int QUADRATIC_PRISM = 26;
   static const int BIQUADRATIC_QUADRATIC_PRISM = 32;
   static const int QUADRATIC_PYRAMID = 27;
   ///@}

   /// @name Arbitrary-order VTK geometric types
   ///@{
   static const int LAGRANGE_SEGMENT = 68;
   static const int LAGRANGE_TRIANGLE = 69;
   static const int LAGRANGE_SQUARE = 70;
   static const int LAGRANGE_TETRAHEDRON = 71;
   static const int LAGRANGE_CUBE = 72;
   static const int LAGRANGE_PRISM = 73;
   static const int LAGRANGE_PYRAMID = 74;
   ///@}
   ///@}
   /// @brief Does the given VTK geometry type describe an arbitrary-order
   /// Lagrange element?
   static bool IsLagrange(int vtk_geom);
   /// @brief Does the given VTK geometry type describe a legacy quadratic
   /// element?
   static bool IsQuadratic(int vtk_geom);
   /// @brief For the given VTK geometry type and number of points, return the
   /// order of the element.
   static int GetOrder(int vtk_geom, int npoints);
};

/// Data array format for VTK and VTU files.
enum class VTKFormat
{
   /// Data arrays will be written in ASCII format.
   ASCII,
   /// Data arrays will be written in binary format. Floating point numbers will
   /// be output with 64 bits of precision.
   BINARY,
   /// Data arrays will be written in binary format. Floating point numbers will
   /// be output with 32 bits of precision.
   BINARY32
};

int AnasysToVTKFormat(std::string vtk_format);

#endif //VTK_H
