#ifndef _SYSTEMEIGEN_H_
#define _SYSTEMEIGEN_H_

#include <complex>
#include "SystemAbstract.h"
#include "petscmat.h"

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

class SystemEigen: public SystemAbstract<std::complex<double> >{
 protected:
  std::list<const FormulationBlock<std::complex<double> >*> formulationB;
  bool general;

  Mat* A;
  Mat* B;

  PetscInt nEigenValues;
  fullVector<std::complex<double> >* eigenValue;
  std::vector<fullVector<std::complex<double> > >* eigenVector;

 public:
  SystemEigen(void);
  virtual ~SystemEigen(void);

  void addFormulationB(const Formulation<std::complex<double> >& formulation);

  virtual size_t getNComputedSolution(void)                             const;
  virtual void   getSolution(fullVector<std::complex<double> >& sol,
                             size_t nSol)                               const;
  virtual void   getSolution(std::map<Dof, std::complex<double> >& sol,
                             size_t nSol)                               const;
  virtual void   getSolution(FEMSolution<std::complex<double> >& feSol,
                             const FunctionSpace& fs,
                             const GroupOfElement& domain)              const;

  bool isGeneral(void) const;
  void getEigenValues(fullVector<std::complex<double> >& eig) const;
  void setNumberOfEigenValues(size_t nEigenValues);

  virtual void assemble(void);
  virtual void solve(void);

  virtual void writeMatrix(std::string fileName,
                           std::string matrixName) const;
 private:
  void assembleCom(SolverMatrix<std::complex<double> >& tmpMat,
                   SolverVector<std::complex<double> >& tmpRHS,
                   const FormulationBlock<std::complex<double> >& formulation,
                   formulationPtr term);
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
*/

#endif
