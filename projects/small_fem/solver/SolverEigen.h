#ifndef _SOLVEREIGEN_H_
#define _SOLVEREIGEN_H_

#include "slepceps.h"
#include "petscmat.h"
#include "fullMatrix.h"
#include "SmallFem.h"

/**
   @class SolverEigen
   @brief This class is an eigenvalue solver

   This class is an eigenvalue solver.

   The eigenvalue problem can be generalized or not:
   @li An eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$
   @li A Generalized Eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$

   In the case of a generalized problem, a second set of Formulation%s
   can be used to assemble matrix @f$\mathbf{B}@f$

   The Solver used is <a href="http://www.grycap.upv.es/slepc/">SLEPc</a>.
 */

class SolverEigen{
 private:
  bool general;

  Mat A;
  Mat B;

  PetscInt       nEigenValues;
  PetscInt       maxIt;
  PetscReal      tol;
  PetscScalar    target;
  EPSProblemType whichProblem;
  EPSWhich       whichEigenpair;

  fullVector<Complex>* eigenValue;
  std::vector<fullVector<Complex> >* eigenVector;

 public:
   SolverEigen(void);
  ~SolverEigen(void);

  bool isGeneral(void)                                   const;
  void getEigenValues(fullVector<Complex>& val)          const;
  void getEigenVector(fullVector<Complex>& vec, int eig) const;
  int  getNComputedSolution(void)                        const;

  void setMatrixA(Mat A);
  void setMatrixB(Mat B);

  void setProblem(std::string type);
  void setWhichEigenpair(std::string type);
  void setTarget(Complex target);
  void setNumberOfEigenValues(int nEigenValues);
  void setMaxIteration(int maxIt);
  void setTolerance(double tol);

  void solve(void);
};


/**
   @fn SolverEigen::SolverEigen
   Instantiates a new SolverEigen
   **

   @fn SolverEigen::~SolverEigen
   Deletes this SolverEigen
   **

   @fn SolverEigen::isGeneral
   @return Returns:
   @li true, if the SolverEigen is a generalized one
   @li false otherwise
   **

   @fn SolverEigen::getEigenValues
   @param eig A vector
   Allocate and populates eig with the eigenvalues of this SolverEigen
   **

   @fn SolverEigen::setProblem(std::string type)
   @param type A string

   Sets the type of eigen problem
   **

   @fn SolverEigen::setNumberOfEigenValues
   @param nEigenValues A natural number

   Sets the number of eigenvalues computed by
   SolverEigen::solve() to the given number
   **

   @fn SolverEigen::setMaxIteration
   @param maxIt An integer

   Sets the maximum number of iteration of SolverEigen::solve()
   to the given value
   **

   @fn SolverEigen::setTolerance
   @param tol A real value

   Sets the convergence tolerance of SolverEigen::solve() to the given value
   **

   @fn SolverEigen::setTarget(Complex target)
   @param target A complex value

   Sets the eigenvalue target of SolverEigen::solve() to the given value
   **

  @fn SolverEigen::setWhichEigenpairs(std::string type)
  @param type A string

   Sets the type of eigenpairs that SolverEigen::solve() will look for
*/

#endif
