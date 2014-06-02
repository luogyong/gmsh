#ifndef _SYSTEMABSTRACT_H_
#define _SYSTEMABSTRACT_H_

#include <string>
#include <list>

#include "DofManager.h"
#include "Formulation.h"
#include "FormulationBlock.h"
#include "FormulationCoupled.h"
#include "FEMSolution.h"

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "fullMatrix.h"

/**
   @interface SystemAbstract
   @brief Common interface for linear systems assemblers

   This is a common interface for linear systems assemblers.
   All SystemAbstract%s are able to assemble the terms of a Formulation.
 */

template<typename scalar>
class SystemAbstract{
 protected:
  typedef
    scalar (FormulationBlock<scalar>::*formulationPtr)(size_t dofI,
                                                       size_t dofJ,
                                                       size_t elementId) const;
 protected:
  static const scalar minusSign;

 protected:
  bool assembled;
  bool solved;

  std::list<const FormulationBlock<scalar>*> formulation;
  DofManager<scalar> dofM;

 public:
  virtual ~SystemAbstract(void);

  bool   isAssembled(void) const;
  bool   isSolved(void)    const;
  size_t getSize(void)     const;

  void addFormulation(const Formulation<scalar>& formulation);
  void constraint(const std::map<Dof, scalar>& constr);

  virtual void assemble(void) = 0;
  virtual void solve(void)    = 0;

  virtual size_t getNComputedSolution(void)                           const = 0;
  virtual void   getSolution(fullVector<scalar>& sol, size_t nSol)    const = 0;
  virtual void   getSolution(std::map<Dof, scalar>& sol, size_t nSol) const = 0;
  virtual void   getSolution(FEMSolution<scalar>& feSol,
                             const FunctionSpace& fs,
                             const GroupOfElement& domain)            const = 0;

  virtual void writeMatrix(std::string fileName,
                           std::string matrixName) const = 0;

 protected:
  void addFormulationBlock(const FormulationBlock<scalar>& formulation,
                           std::list<const FormulationBlock<scalar>*>& fList);

  void addFormulationCoupled(const FormulationCoupled<scalar>& formulation,
                             std::list<const FormulationBlock<scalar>*>& fList);

  void assemble(SolverMatrix<scalar>& A,
                SolverVector<scalar>& b,
                size_t elementId,
                const std::vector<Dof>& dofField,
                const std::vector<Dof>& dofTest,
                formulationPtr& term,
                const FormulationBlock<scalar>& formulation);

  void assembleRHSOnly(SolverVector<scalar>& b,
                       size_t elementId,
                       const std::vector<Dof>& dofTest,
                       const FormulationBlock<scalar>& formulation);
};


/**
   @fn SystemAbstract::~SystemAbstract
   Deletes this SystemAbstract
   **

   @fn SystemAbstract::isAssembled
   @return Returns:
   @li true, if the system has been assembled
   @li false otherwise
   **

   @fn SystemAbstract::isSolved
   @return Returns:
   @li true, if the system has been solved
   @li false otherwise
   **

   @fn SystemAbstract::getSize
   @return Returns the number of unknowns in this linear system
   **

   @fn SystemAbstract::addFormulation(const Formulation<scalar>& formulation)
   @param formulation A Formulation
   Adds the given Formulation to the Formulation%s that will be assembled
   **

   @fn SystemAbstract::constraint
   @param constr A map of Dof%s and scalar
   Constraints this SystemAbstract with the given Dof%s
   and their associated values
   **

   @fn SystemAbstract::assemble(void)
   Assembles this linear system

   All the given Formulation%s will be processed
   @see SytemAbstract::addFormulation
   **

   @fn SystemAbstract::solve(void)
   Solves this linear system
   **

   @fn SystemAbstract::getNComputedSolution(void)
   @return The number of computed solution by SystemAbstract::solve()
   **

   @fn SystemAbstract::getSolution(fullVector<scalar>&, size_t) const = 0
   @param sol A vector
   @param nSol An integer
   Allocates and populates the given vector with the nSolth solution vector
   computed by SystemAbstract::solve()
   **

   @fn SystemAbstract::getSolution(std::map<Dof, scalar>&, size_t) const = 0
   @param sol A map mapping a Dof to a scalar
   @param nSol An integer
   Takes every Dof in the given map and set its assoicated value to
   the nSolth solution of this SystemAbstract
   **

   @fn SystemAbstract::getSolution(FEMSolution<scalar>&,const FunctionSpace&,const GroupOfElement&) const = 0
   @param feSol A FEMSolution
   @param fs A FunctionSpace
   @param domain A Domain

   Adds to the given FEMSolution the computed finite element solutions,
   for the given FunctionSpace and the given domain.
   If no solution has been computed, and Exception is throw.
   **

   @fn SystemAbstract::writeMatrix
   @param fileName A string
   @param matrixName A string

   Writes this system matrix in Octave/Matlab format,
   with the given name and into the given file
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SystemAbstractInclusion.h"

#endif
