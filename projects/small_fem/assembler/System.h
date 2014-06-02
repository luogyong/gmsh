#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "SystemAbstract.h"
#include "SolverMUMPS.h"

/**
   @class System
   @brief This class assembles a linear system

   This class assembles a linear system.

   The Solver used is <a href="http://graal.ens-lyon.fr/MUMPS/index.php">MUMPS
   </a>.
 */

template<typename scalar>
class System: public SystemAbstract<scalar>{
 protected:
  SolverMatrix<scalar>* A;
  SolverVector<scalar>* b;
  fullVector<scalar>*   x;

  SolverMUMPS<scalar> solver;

 public:
  System(void);
  virtual ~System(void);

  virtual size_t getNComputedSolution(void)                           const;
  virtual void   getSolution(fullVector<scalar>& sol, size_t nSol)    const;
  virtual void   getSolution(std::map<Dof, scalar>& sol, size_t nSol) const;
  virtual void   getSolution(FEMSolution<scalar>& feSol,
                             const FunctionSpace& fs,
                             const GroupOfElement& domain)            const;

  virtual void assemble(void);
  virtual void solve(void);

  void assembleAgainRHS(void);
  void solveAgain(void);

  virtual void writeMatrix(std::string fileName,
                           std::string matrixName) const;
};


/**
   @fn System::System
   Instantiates a new System
   **

   @fn System::~System
   Deletes this System
   **

   @fn System::assembleAgainRHS
   This method reassemble the right hand sides of the Formulation%s given by
   System::addFormulation().

   This method works only if the full system was already assembled by
   System::assemble().
   **

   @fn System::solveAgain
   This method solve the new system given by System::assembleAgainRHS.

   This method works only if the full sysem was already solved by
   System::solve().
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SystemInclusion.h"


#endif
