//
// Created by huangkai on 25-1-21.
//

#ifndef BILINEARFORM_H
#define BILINEARFORM_H
#include "AnsysLoad.h"
#include "NumericalAlgebra.hpp"
#include <unordered_set>
#include <functional>

enum class AssemblyLevel
{
    /// In the case of a BilinearForm LEGACY corresponds to a fully assembled
    /// form, i.e. a global sparse matrix in MFEM, Hypre or PETSC format.
    /// In the case of a NonlinearForm LEGACY corresponds to an operator that
    /// is fully evaluated on the fly.
    /// This assembly level is ALWAYS performed on the host.
    LEGACY = 0,
    /// @deprecated Use LEGACY instead.
    LEGACYFULL = 0,
    /// Fully assembled form, i.e. a global sparse matrix in MFEM format. This
     /// assembly is compatible with device execution.
    FULL,
    /// Form assembled at element level, which computes and stores dense element
     /// matrices.
    ELEMENT,
    /// Partially-assembled form, which computes and stores data only at
     /// quadrature points.
    PARTIAL,
    /// "Matrix-free" form that computes all of its action on-the-fly without any
     /// substantial storage.
    NONE,
 };

class BilinearForm {
 public:
   /// Creates bilinear form associated with FE space @a *f.
   /** The pointer @a f is not owned by the newly constructed object. */
   BilinearForm(AnsysLoad *f):fes(f){stif=nullptr;mass=nullptr;damp=nullptr;}

   /// Get the size of the BilinearForm as a square matrix.
   size_t Size() const { return fes->getNumberOfDof();}

   void assemble(const FloatArray& x, int skip_zeros = 1);
   void assembleStif(FloatArray x, int skip_zeros = 1);
   void assembleMass(FloatArray x, int skip_zeros = 1);
   void assemble(FloatArray x, std::function<double(double x)> StifRho, std::function<double(double x)> MassRho, int skip_zeros = 1);
   void assemble(FloatArray x, std::function<double(double x)> StifRho, std::function<double(double x)> MassRho, std::function<double(double x)> DampRho, int skip_zeros = 1);
   void assembleStif(FloatArray x,std::function<double(double x)> StifRho, int skip_zeros = 1);
   void assembleMass(FloatArray x,std::function<double(double x)> MassRho, int skip_zeros = 1);
   void assembleDamp(FloatArray x,std::function<double(double x)> DampRho, int skip_zeros = 1);
   void EliminateEssentialBC(std::unordered_set<int> ConstrainedDofs, int d = 0);
   virtual void Finalize (int skip_zeros = 1);
   void convertToMKL();
   SparseMatrix &SpStif() { return *stif; }
   SparseMatrix &SpMass() { return *mass; }
   SparseMatrix &SpDamp() { return *damp; }
   ~BilinearForm() {
    delete stif;
    delete mass;
    delete damp;
   };
private:
   /// Sparse matrix $ K,M $ to be associated with the form. Owned.
   SparseMatrix *stif;
   SparseMatrix *mass;
   SparseMatrix *damp;


   /// FE space on which the form lives. Not owned.
   AnsysLoad *fes;

   /** @brief The ::AssemblyLevel of the form (AssemblyLevel::LEGACY,
       AssemblyLevel::FULL, AssemblyLevel::ELEMENT, AssemblyLevel::PARTIAL) */
   AssemblyLevel assembly;



};



#endif //BILINEARFORM_H
