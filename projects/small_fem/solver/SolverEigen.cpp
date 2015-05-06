#include "Exception.h"
#include "SolverEigen.h"

using namespace std;

SolverEigen::SolverEigen(void){
  // Default //
  general        = false;
  nEigenValues   = 10;
  maxIt          = 100;
  tol            = 1E-6;
  target         = 0;
  whichProblem   = EPS_NHEP;               // Non-hermitian problem (Ax=lx)
  whichEigenpair = EPS_SMALLEST_MAGNITUDE; // Smallest in magnitude
  eigenValue     = NULL;
  eigenVector    = NULL;
}

SolverEigen::~SolverEigen(void){
  if(eigenVector)
    delete eigenVector;

  if(eigenValue)
    delete eigenValue;
}

bool SolverEigen::isGeneral(void) const{
  return general;
}

void SolverEigen::getEigenValues(fullVector<Complex>& val) const{
  val.setAsProxy(*eigenValue, 0, eigenValue->size());
}

void SolverEigen::getEigenVector(fullVector<Complex>& vec, int eig) const{
  if(eig < 0)
    throw Exception("SolverEigen::getEigenVector(): eig must be greater "
                    "or equal to 0");

  if(eig >= getNComputedSolution())
    throw Exception("SolverEigen::getEigenVector(): eig must be smaller "
                    "than SolverEigen::getNComputedSolution()");

  vec.setAsProxy((*eigenVector)[eig], 0, (*eigenVector)[eig].size());
}

int SolverEigen::getNComputedSolution(void) const{
  return nEigenValues;
}

void SolverEigen::setMatrixA(Mat A){
  int nRow;
  int nCol;
  MatGetSize(A, &nRow, &nCol);

  if(nRow != nCol)
    throw Exception("SolverEigen::setMatrixA(): matrix is not square (%d, %d)",
                    nRow, nCol);

  this->A = A;
}

void SolverEigen::setMatrixB(Mat B){
  int nRow;
  int nCol;
  MatGetSize(B, &nRow, &nCol);

  if(nRow != nCol)
    throw Exception("SolverEigen::setMatrixB(): matrix is not square (%d, %d)",
                    nRow, nCol);

  general = true;
  this->B = B;
}

void SolverEigen::setProblem(std::string type){
  if(!type.compare("gen_non_hermitian"))
    this->whichProblem = EPS_GNHEP;

  else if(!type.compare("non_hermitian"))
    this->whichProblem = EPS_NHEP;

  else if(!type.compare("pos_gen_non_hermitian"))
    this->whichProblem = EPS_PGNHEP;

  else
    throw Exception("SolverEigen::setProblem(): legal values are: "
                    "gen_non_hermitian; non_hermitian; pos_gen_non_hermitian "
                    "(current value is %s)",
                    type.c_str());
}

void SolverEigen::setWhichEigenpair(std::string type){
  if(!type.compare("smallest_magnitude"))
    this->whichEigenpair = EPS_SMALLEST_MAGNITUDE;

  else if(!type.compare("target_real"))
    this->whichEigenpair = EPS_TARGET_REAL;

  else if(!type.compare("target_magnitude"))
    this->whichEigenpair = EPS_TARGET_MAGNITUDE;

  else
    throw Exception("SolverEigen::setWhichEigenpairs(): legal values are: "
                    "smallest_magnitude; target_real; target_magnitude "
                    "(current value is %s)",
                    type.c_str());
}

void SolverEigen::setTarget(Complex target){
  this->target = target.real() + PETSC_i * target.imag();
}

void SolverEigen::setNumberOfEigenValues(int nEigenValues){
  if(nEigenValues <= 0)
    throw Exception("SolverEigen::setNumberOfEigenValues(): "
                    "I cannot compute less than zero eigenvalues");

  this->nEigenValues = nEigenValues;
}

void SolverEigen::setMaxIteration(int maxIt){
  if(maxIt <= 0)
    throw Exception("SolverEigen::setMaxIteration(): "
                    "the maximum number of iterations "
                    "must be bigger than zero");

  this->maxIt = (PetscInt)(maxIt);
}

void SolverEigen::setTolerance(double tol){
  this->tol = (PetscReal)(tol);
}

void SolverEigen::solve(void){
  // Checks //
  // nEigenValues
  if(nEigenValues == 0)
    throw
      Exception("SolverEigen: the number of eigenvalues to compute is zero");

  // Matrix A
  if(A == NULL)
    throw Exception("SolverEigen::solver(): no A matrix given");

  // Matrices sizes
  if(general){
    int sizeA;
    int sizeB;

    MatGetSize(A, &sizeA, NULL);
    MatGetSize(B, &sizeB, NULL);

    if(sizeA != sizeB)
      throw Exception("SolverEigen: cannot solve generalized problem "
                      "with matrices with different sizes");
  }

  // Asking for too much eigenvalues ?
  int size;
  MatGetSize(A, &size, NULL);

  if(nEigenValues > size)
    nEigenValues = size;

  // Set up problem //
  EPS solver;
  EPSCreate(MPI_COMM_WORLD, &solver);

  // Operator(s)
  if(general)
    EPSSetOperators(solver, A, B);
  else
    EPSSetOperators(solver, A, NULL);

  // Set problem type
  EPSSetProblemType(solver, whichProblem);

  // Set Options
  EPSSetDimensions(solver, nEigenValues, PETSC_DECIDE, PETSC_DECIDE);
  EPSSetTolerances(solver, tol, maxIt);

  // Which Eigenpair
  EPSSetWhichEigenpairs(solver, whichEigenpair);
  EPSSetTarget(solver, target);

  // Use Krylov Schur Solver //
  EPSSetType(solver, "krylovschur");

  // Spectral transform (if needed) //
  if(target != 0.0){
    KSP ksp; // Krylov subspace solver
    PC  pc;  // Preconditioner
    ST  st;  // Spectral transform

    EPSGetST(solver, &st);
    STSetType(st, "sinvert");
    STGetKSP(st, &ksp);

    STSetMatMode(st, ST_MATMODE_INPLACE);
    PetscOptionsSetValue("-mat_mumps_icntl_28", "2"); // MUMPS Parallel analysis
    PetscOptionsSetValue("-mat_mumps_icntl_7",  "5"); // METIS
    PetscOptionsSetValue("-mat_mumps_icntl_29", "2"); // ParMETIS

    KSPSetType(ksp, "preonly");
    KSPGetPC(ksp, &pc);
    PCSetType(pc, "lu");
    PCFactorSetMatSolverPackage(pc, "mumps");

    STSetFromOptions(st);
  }

  // Override with PETSc Database & Solve //
  EPSSetFromOptions(solver);
  EPSSolve(solver);

  // Wait for everything to be ok //
  MPI_Barrier(MPI_COMM_WORLD);

  // Get Solution //
  VecScatter   scat;
  PetscScalar  lambda;
  PetscScalar* x;
  Vec          xPetscDist;
  Vec          xPetscSeq;

  MatGetVecs(A, PETSC_NULL, &xPetscDist);
  VecCreate(MPI_COMM_SELF, &xPetscSeq);
  VecSetType(xPetscSeq, VECSEQ);
  VecSetSizes(xPetscSeq, size, size);
  VecScatterCreateToAll(xPetscDist, &scat, NULL);

  EPSGetConverged(solver, &nEigenValues);

  eigenValue  = new fullVector<Complex>(nEigenValues);
  eigenVector = new vector<fullVector<Complex> >(nEigenValues);

  for(PetscInt i = 0; i < nEigenValues; i++){
    EPSGetEigenpair(solver, i, &lambda, NULL, xPetscDist, NULL);

    VecScatterBegin(scat, xPetscDist, xPetscSeq, INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(  scat, xPetscDist, xPetscSeq, INSERT_VALUES,SCATTER_FORWARD);

    VecGetArray(xPetscSeq, &x);

    (*eigenVector)[i].resize(size);
    for(int j = 0; j < size; j++)
      (*eigenVector)[i](j) = x[j];

    (*eigenValue)(i) = lambda;
  }

  // Clear //
  VecDestroy(&xPetscDist);
  VecDestroy(&xPetscSeq);
  VecScatterDestroy(&scat);
  EPSDestroy(&solver);

  // Wait for everything to be ok //
  MPI_Barrier(MPI_COMM_WORLD);
}
