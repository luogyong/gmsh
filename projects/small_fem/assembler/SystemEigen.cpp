#include "slepceps.h"
#include "SystemEigen.h"

#include "SmallFem.h"

using namespace std;

SystemEigen::SystemEigen(void){
  // Is the Problem a General EigenValue Problem ? //
  general = false;

  // Dof Manager (distributed) //
  dofM = new DofManager<Complex>(false);

  // Init //
  A           = NULL;
  B           = NULL;
  eigenValue  = NULL;
  eigenVector = NULL;

  // Per-thread non-zero term //
  #pragma omp parallel
  {
    #pragma omp master
    {
      this->nNZCount.resize(omp_get_num_threads());
      this->nNZCountB.resize(omp_get_num_threads());
    }
  }

  for(size_t i = 0; i < this->nNZCount.size(); i++)
    this->nNZCount[i] = 0;

  for(size_t i = 0; i < this->nNZCountB.size(); i++)
    this->nNZCountB[i] = 0;

  // Default //
  nEigenValues   = 10;
  maxIt          = 100;
  tol            = 1E-6;
  target         = 1;
  whichEigenpair = string("smallest_magnitude");

  // The SystemEigen is not assembled and not solved//
  assembled    = false;
  solved       = false;
}

SystemEigen::~SystemEigen(void){
  delete dofM;

  if(eigenVector)
    delete eigenVector;

  if(eigenValue)
    delete eigenValue;

  if(A){
    MatDestroy(A);
    delete A;
  }

  if(B){
    MatDestroy(B);
    delete B;
  }
}

void SystemEigen::
addFormulationB(const Formulation<Complex>& formulation){
  // Add FormulationBlocks to list B (this->formulationB) //

  // Is this a FormulationBlock ?
  if(formulation.isBlock())
    addFormulationBlock
      (static_cast<const FormulationBlock<Complex>&>(formulation),
       formulationB);

  else
    addFormulationCoupled
      (static_cast<const FormulationCoupled<Complex>&>(formulation),
       formulationB);

  // This EigenSystem is general
  general = true;
}

void SystemEigen::
countCom(std::list<const FormulationBlock<Complex>*>::iterator it,
         std::list<const FormulationBlock<Complex>*>::iterator end,
         std::vector<size_t>& nNZCount){

  // Iterate on Formulations //
  for(; it != end; it++){
    // Get All Dofs (Field & Test) per Element
    const vector<vector<Dof> >& dofField =
      (*it)->field().getKeys((*it)->domain());
    const vector<vector<Dof> >& dofTest  =
      (*it)->test().getKeys((*it)->domain());

    // Count
    const size_t E = dofField.size(); // Should be equal to dofTest.size().?.

    #pragma omp parallel for
    for(size_t i = 0; i < E; i++)
      nNZCount[omp_get_thread_num()] =
        SystemAbstract<Complex>::countTerms(nNZCount[omp_get_thread_num()],
                                            i, dofField[i], dofTest[i], **it);
  }
}

void SystemEigen::
assembleCom(std::list<const FormulationBlock<Complex>*>::iterator it,
            std::list<const FormulationBlock<Complex>*>::iterator end,
            SolverMatrix<Complex>& tmpMat){

  // Iterate on Formulations //
  for(; it != end; it++){
    // Get All Dofs (Field & Test) per Element
    const vector<vector<Dof> >& dofField =
      (*it)->field().getKeys((*it)->domain());
    const vector<vector<Dof> >& dofTest  =
      (*it)->test().getKeys((*it)->domain());

    // Assemble Systems
    const size_t E = dofField.size();   // Should be equal to dofTest.size().?.

    #pragma omp parallel for
    for(size_t i = 0; i < E; i++)
      SystemAbstract<Complex>::
        assembleLHSOnly(tmpMat, i, dofField[i], dofTest[i], **it);
  }
}

Mat* SystemEigen::toPetsc(SolverMatrix<Complex>* tmp, size_t size){
  // MPI range //
  int nProc;
  int myProc;
  int myProcSize;

  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);

  myProcSize = procSize[myProc];

  // Get tmp data //
  int*       row;
  int*       col;
  Complex* value;
  size_t tmpSize = tmp->get(&row, &col, &value, true);

  // Get matrix ownership //
  vector<size_t> owner;
  getOwnership(procSize, owner);

  // Full sparsity of PETSc matrix //
  int* nonZeroDiagFull    = new int[size];
  int* nonZeroOffDiagFull = new int[size];

  for(size_t i = 0; i < size; i++)
    nonZeroDiagFull[i] = 0;

  for(size_t i = 0; i < size; i++)
    nonZeroOffDiagFull[i] = 0;

  petscSparsity(nonZeroDiagFull, row, col, tmpSize,
                procMinRange, procMaxRange, owner, true);

  petscSparsity(nonZeroOffDiagFull, row, col, tmpSize,
                procMinRange, procMaxRange, owner, false);

  // Local sparcity of PETSc matrix //
  int* nonZeroDiag    = new int[myProcSize];
  int* nonZeroOffDiag = new int[myProcSize];

  for(int i = 0; i < nProc; i++)
    MPI_Reduce(&nonZeroDiagFull[procMinRange[i]]   , nonZeroDiag   ,procSize[i],
               MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);

  for(int i = 0; i < nProc; i++)
    MPI_Reduce(&nonZeroOffDiagFull[procMinRange[i]], nonZeroOffDiag,procSize[i],
               MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);

  // Clear //
  delete[] nonZeroDiagFull;
  delete[] nonZeroOffDiagFull;

  // Maximum row size (on and off diagonal) //
  int offSize = size - myProcSize;          // Off diag maximum size

  for(int i = 0; i < myProcSize; i++)
    if(nonZeroDiag[i] > myProcSize)         // Diag is limited to proc size
      nonZeroDiag[i] = myProcSize;

  for(int i = 0; i < myProcSize; i++)
    if(nonZeroOffDiag[i] > offSize)         // Off diag is limited to offSize
      nonZeroDiag[i] = offSize;

  // Copy tmp into PETSc matrix //
  Mat* M = new Mat;

  MatCreateAIJ(MPI_COMM_WORLD,
               myProcSize, myProcSize, size, size,
               42, nonZeroDiag, 42, nonZeroOffDiag, M);

  petscSerialize(row, col, value, tmpSize, *M);

  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (*M, MAT_FINAL_ASSEMBLY);
  /*
  PetscViewer viewer;

  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "B2.m", &viewer);
  PetscObjectSetName((PetscObject)(*M), "B2");
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

  MatView(*M, viewer);

  PetscViewerDestroy(&viewer);
  */
  // Free //
  delete[] nonZeroDiag;
  delete[] nonZeroOffDiag;

  // Wait for everything to be ok //
  MPI_Barrier(MPI_COMM_WORLD);

  // Retrun //
  return M;
}

void SystemEigen::assemble(void){
  // Enumerate Dofs in DofManager //
  dofM->generateGlobalIdSpace();

  // Count FE terms //
  countCom(formulation.begin() , formulation.end() , nNZCount);
  countCom(formulationB.begin(), formulationB.end(), nNZCountB);

  // Alloc Temp Sparse Matrices (not with PETSc) //
  const size_t sizeLocal  = dofM->getLocalSize();
  const size_t sizeGlobal = dofM->getGlobalSize();

  SolverMatrix<Complex>* tmpA;
  SolverMatrix<Complex>* tmpB;

  tmpA = new SolverMatrix<Complex>(sizeLocal ,sizeLocal, nNZCount);
  tmpB = new SolverMatrix<Complex>(sizeLocal ,sizeLocal, nNZCountB);

  // MPI size //
  int nProc;
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  getProcSize(sizeGlobal, nProc, procSize);
  getProcMinRange(procSize, procMinRange);
  getProcMaxRange(procSize, procMaxRange);

  // Assemble Formulations A //
  assembleCom(formulation.begin(), formulation.end(), *tmpA);
  A = toPetsc(tmpA, sizeGlobal); // Allocates A
  delete tmpA;

  // Assemble Formulations B //
  if(general){
    assembleCom(formulationB.begin(), formulationB.end(), *tmpB);
    B = toPetsc(tmpB, sizeGlobal); // Allocates B
    delete tmpB;
  }

  else{
    delete tmpB;
  }

  // The SystemEigen is assembled //
  assembled = true;

  // Wait for everything to be ok //
  MPI_Barrier(MPI_COMM_WORLD);
}

void SystemEigen::solve(void){
  // Check nEigenValues
  if(!nEigenValues)
    throw
      Exception("The number of eigenvalues to compute is zero");

  // Is the SystemEigen assembled ? //
  if(!assembled)
    assemble();

  // Build Solver //
  EPS solver;
  EPSCreate(MPI_COMM_WORLD, &solver);

  if(general)
    EPSSetOperators(solver, *A, *B);
  else
    EPSSetOperators(solver, *A, NULL);

  if(general)
    EPSSetProblemType(solver, EPS_GNHEP);
  else
    EPSSetProblemType(solver, EPS_NHEP);

  // Set Options //
  EPSSetDimensions(solver, nEigenValues, PETSC_DECIDE, PETSC_DECIDE);
  EPSSetTolerances(solver, tol, maxIt);

  // Which Eigenpair //
  if(!whichEigenpair.compare("smallest_magnitude")){
    EPSSetWhichEigenpairs(solver, EPS_SMALLEST_MAGNITUDE);
    EPSSetTarget(solver, target);
  }

  else if(!whichEigenpair.compare("target_real")){
    EPSSetWhichEigenpairs(solver, EPS_TARGET_REAL);
    EPSSetTarget(solver, target);
  }

  else if(!whichEigenpair.compare("target_magnitude")){
    EPSSetWhichEigenpairs(solver, EPS_TARGET_MAGNITUDE);
    EPSSetTarget(solver, target);
  }

  else
    throw
      Exception("%s -- %s %s %s", "SystemEigen::solve",
                                  "set SystemEigen::setWhichEigenpair()",
                                  "to a legal value: set to",
                                  whichEigenpair.c_str());

  // Use Krylov Schur Solver and MUMPS //
  KSP ksp; // Krylov subspace solver
  PC  pc;  // Preconditioner
  ST  st;  // Spectral transform

  EPSSetType(solver, "krylovschur");

  EPSGetST(solver, &st);
  STSetType(st, STSINVERT);
  STGetKSP(st, &ksp);

  KSPSetType(ksp, "preonly");
  KSPGetPC(ksp, &pc);
  PCSetType(pc, "lu");
  PCFactorSetMatSolverPackage(pc, "mumps");

  // Override with PETSc Database //
  EPSSetFromOptions(solver);
  STSetFromOptions(st);

  // Solve //
  EPSSolve(solver);

  // Wait for everything to be ok //
  MPI_Barrier(MPI_COMM_WORLD);

  // Get Solution //
  const size_t size = dofM->getGlobalSize();

  VecScatter   scat;
  PetscScalar  lambda;
  PetscScalar* x;
  Vec          xPetscDist;
  Vec          xPetscSeq;

  MatGetVecs(*A, PETSC_NULL, &xPetscDist);
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
    VecScatterEnd(scat, xPetscDist, xPetscSeq, INSERT_VALUES,SCATTER_FORWARD);

    VecGetArray(xPetscSeq, &x);

    (*eigenVector)[i].resize(size);
    for(size_t j = 0; j < size; j++)
      (*eigenVector)[i](j) = x[j];

    (*eigenValue)(i) = lambda;
  }

  VecDestroy(&xPetscDist);
  VecDestroy(&xPetscSeq);
  VecScatterDestroy(&scat);
  EPSDestroy(&solver);

  // System solved ! //
  solved = true;

  // Wait for everything to be ok //
  MPI_Barrier(MPI_COMM_WORLD);
}

bool SystemEigen::isGeneral(void) const{
  return general;
}

void SystemEigen::getEigenValues(fullVector<Complex>& eig) const{
  eig.setAsProxy(*eigenValue, 0, eigenValue->size());
}

void SystemEigen::
setNumberOfEigenValues(size_t nEigenValues){
  const size_t nDof = dofM->getGlobalSize();

  if(nEigenValues > nDof)
    throw
      Exception
      ("I can't compute more Eigenvalues (%d) than the number of unknowns (%d)",
       nEigenValues, nDof);

  else
    this->nEigenValues = nEigenValues;
}

void SystemEigen::setMaxIteration(size_t maxIt){
  this->maxIt = (PetscInt)(maxIt);
}

void SystemEigen::setTolerance(double tol){
  this->tol = (PetscReal)(tol);
}

void SystemEigen::setTarget(Complex target){
  this->target = target.real() + PETSC_i * target.imag();
}

void SystemEigen::setWhichEigenpairs(std::string type){
  this->whichEigenpair = type;
}

size_t SystemEigen::getNComputedSolution(void) const{
  return nEigenValues;
}

void SystemEigen::getSolution(fullVector<Complex>& sol,
                              size_t nSol) const{
  sol.setAsProxy((*eigenVector)[nSol], 0, (*eigenVector)[nSol].size());
}

void SystemEigen::getSolution(std::map<Dof, Complex>& sol,
                              size_t nSol) const{
  // Get All Dofs
  map<Dof, Complex>::iterator it  = sol.begin();
  map<Dof, Complex>::iterator end = sol.end();

  // Loop on Dofs and set Values
  for(; it != end; it++){
    size_t gId = dofM->getGlobalId(it->first);

    if(gId == DofManager<Complex>::isFixedId())
      it->second = dofM->getValue(it->first);

    else
      it->second = (*eigenVector)[nSol](gId);
  }
}

void SystemEigen::getSolution(FEMSolution<Complex>& feSol,
                              const FunctionSpace& fs,
                              const GroupOfElement& domain) const{
  // Solved ?
  if(!solved)
    throw Exception("System: addSolution -- System not solved");

  // Coefficients //
  // Get Dofs
  set<Dof> dof;
  fs.getKeys(domain, dof);

  // Get Coefficient
  const set<Dof>::iterator   end = dof.end();
  set<Dof>::iterator         it  = dof.begin();
  map<Dof, Complex> coef;

  for(; it != end; it++)
    coef.insert(pair<Dof, Complex>(*it, 0));

  // Iterate on Solutions //
  for(int i = 0; i < nEigenValues; i++){
    // Populate Map
    getSolution(coef, i);

    // FEMSolution
    feSol.addCoefficients(i, 0, domain, fs, coef);
  }
}

void SystemEigen::
getSolution(FEMSolution<Complex>& feSol,
            const FunctionSpace& fs,
            const std::vector<const GroupOfElement*>& domain) const{
  // Get size
  const size_t size = domain.size();

  // Get solution for each domain
  for(size_t i = 0; i < size; i++)
    getSolution(feSol, fs, *domain[i]);
}

void SystemEigen::writeMatrix(string fileName,
                              string matrixName) const{
  throw Exception("SystemEigen::writeMatrix -- not implemented");
}
