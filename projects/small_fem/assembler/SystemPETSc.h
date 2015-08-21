#ifndef _SYSTEMPETSC_H_
#define _SYSTEMPETSC_H_

#include "SystemAbstract.h"
#include "petscmat.h"
#include "petscvec.h"

/**
   @class SystemPETSc
   @brief This class assembles a linear system

   This class assembles a linear system.

   The Solver used is PETSc.
 */

template<typename scalar>
class SystemPETSc: public SystemAbstract<scalar>{
 protected:
  Mat A;
  Vec b;
  Vec xPetsc;
  fullVector<scalar>* x;

 public:
  SystemPETSc(void);
  virtual ~SystemPETSc(void);

  virtual size_t getNComputedSolution(void)                           const;
  virtual void   getSolution(fullVector<scalar>& sol, size_t nSol)    const;
  virtual void   getSolution(std::map<Dof, scalar>& sol, size_t nSol) const;
  virtual void   getSolution(FEMSolution<scalar>& feSol,
                             const FunctionSpace& fs,
                             const GroupOfElement& domain)            const;

  virtual void   getSolution(FEMSolution<scalar>& feSol,
                             const FunctionSpace& fs,
                             const std::vector<const GroupOfElement*>& domain)
                                                                      const;
  virtual void assemble(void);
  virtual void solve(void);

  virtual void writeMatrix(std::string fileName,
                           std::string matrixName) const;

 private:
  PetscScalar* toPetscScalar(scalar* a, size_t size);
  void getSolution(void);
};


/**
   @fn SystemPETSc::SystemPETSc
   Instantiates a new SystemPETSc
   **

   @fn SystemPETSc::~SystemPETSc
   Deletes this SystemPETSc
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SystemPETScInclusion.h"


#endif
