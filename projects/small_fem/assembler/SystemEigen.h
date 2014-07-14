#ifndef _SYSTEMEIGEN_H_
#define _SYSTEMEIGEN_H_

#include <complex>
#include "SystemAbstract.h"
#include "petscmat.h"

#include "SmallFem.h"

/**
   @class SystemEigen
   @brief This class assembles an eigenvalue system

   This class assembles an eigenvalue system.

   The eigenvalue problem can be generalized or not:
   @li An eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$
   @li A Generalized Eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$

   In the case of a generalized problem, a second set of Formulation%s
   can be used to assemble matrix @f$\mathbf{B}@f$

   The Solver used is <a href="http://www.grycap.upv.es/slepc/">SLEPc</a>.
 */

class SystemEigen: public SystemAbstract<Complex>{
 protected:
  std::list<const FormulationBlock<Complex>*> formulationB;
  bool general;

  std::vector<size_t> procSize;
  std::vector<size_t> procMinRange;
  std::vector<size_t> procMaxRange;

  std::vector<size_t> nNZCountB;

  Mat* A;
  Mat* B;

  PetscInt    nEigenValues;
  PetscInt    maxIt;
  PetscReal   tol;
  PetscScalar target;
  std::string whichEigenpair;

  fullVector<Complex>* eigenValue;
  std::vector<fullVector<Complex> >* eigenVector;

 public:
  SystemEigen(void);
  virtual ~SystemEigen(void);

  void addFormulationB(const Formulation<Complex>& formulation);

  virtual size_t getNComputedSolution(void)                             const;
  virtual void   getSolution(fullVector<Complex>& sol, size_t nSol)     const;
  virtual void   getSolution(std::map<Dof, Complex >& sol, size_t nSol) const;
  virtual void   getSolution(FEMSolution<Complex>& feSol,
                             const FunctionSpace& fs,
                             const GroupOfElement& domain)              const;

  bool isGeneral(void) const;
  void getEigenValues(fullVector<Complex>& eig) const;

  void setNumberOfEigenValues(size_t nEigenValues);
  void setMaxIteration(size_t maxIt);
  void setTolerance(double tol);
  void setTarget(Complex target);
  void setWhichEigenpairs(std::string type);

  virtual void assemble(void);
  virtual void solve(void);

  virtual void writeMatrix(std::string fileName,
                           std::string matrixName) const;
 private:
  void countCom(std::list<const FormulationBlock<Complex>*>::iterator it,
                std::list<const FormulationBlock<Complex>*>::iterator end,
                std::vector<size_t>& nNZCount);

  void assembleCom(std::list<const FormulationBlock<Complex>*>::iterator it,
                   std::list<const FormulationBlock<Complex>*>::iterator end,
                   SolverMatrix<Complex>& tmpMat,
                   SolverVector<Complex>& tmpRHS);

  Mat* toPetsc(SolverMatrix<Complex>* tmp, size_t size, int myProc);
};


/**
   @fn SystemEigen::SystemEigen
    Instantiates a new SystemEigen
   ***

   @fn SystemEigen::~SystemEigen
   Deletes this SystemEigen
   **

   @fn SystemEigen::addFormulationB
   @param formulation A Formulation

   Adds the given Formulation to the Formulation%s that will be assembled
   for matrix @f$\mathbf{B}@f$

   If at least of Formulation is added with SystemEigen::addFormulationB(),
   this SystemEigen becomes generalized
   **

   @fn SystemEigen::isGeneral
   @return Returns:
   @li true, if the SystemEigen is a generalized one
   @li false otherwise
   **

   @fn SystemEigen::getEigenValues
   @param eig A vector
   Allocate and populates eig with the eigenvalues of this SystemEigen
   **

   @fn SystemEigen::setNumberOfEigenValues
   @param nEigenValues A natural number

   Sets the number of eigenvalues computed by
   SystemEigen::solve() to the given number
   **

   @fn SystemEigen::setMaxIteration
   @param maxIt An integer

   Sets the maximum number of iteration of SystemEigen::solve()
   to the given value
   **

   @fn SystemEigen::setTolerance
   @param tol A real value

   Sets the convergence tolerance of SystemEigen::solve() to the given value
   **

   @fn SystemEigen::setTarget(Complex target)
   @param target A complex value

   Sets the eigenvalue target of SystemEigen::solve() to the given value
   **

  @fn SystemEigen::setWhichEigenpairs(std::string type)
   @param type A string

   Sets the type of eigenpairs that SystemEigen::solve() will look for
*/

#endif
