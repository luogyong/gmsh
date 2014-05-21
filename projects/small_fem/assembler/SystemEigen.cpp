#include "slepceps.h"
#include "SystemEigen.h"

using namespace std;

SystemEigen::SystemEigen(void){
  // Is the Problem a General EigenValue Problem ? //
  general = false;

  // Init //
  A           = NULL;
  B           = NULL;
  eigenValue  = NULL;
  eigenVector = NULL;

  // Default //
  whichEigenpair = string("smallest_magnitude");
  target         = 0;
  nEigenValues   = 10;

  // The SystemEigen is not assembled and not solved//
  assembled    = false;
  solved       = false;
}

SystemEigen::~SystemEigen(void){
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
addFormulationB(const Formulation<complex<double> >& formulation){
  // Add FormulationBlocks to list B (this->formulationB) //

  // Is this a FormulationBlock ?
  if(formulation.isBlock())
    addFormulationBlock
      (static_cast<const FormulationBlock<complex<double> >&>(formulation),
       formulationB);

  else
    addFormulationCoupled
      (static_cast<const FormulationCoupled<complex<double> >&>(formulation),
       formulationB);

  // This EigenSystem is general
  general = true;
}

void SystemEigen::assembleCom(SolverMatrix<complex<double> >& tmpMat,
                              SolverVector<complex<double> >& tmpRHS,
                              const FormulationBlock<complex<double> >&
                                                                    formulation,
                              formulationPtr term){
  // Get All Dofs (Field & Test) per Element //
  vector<vector<Dof> > dofField;
  vector<vector<Dof> > dofTest;
  formulation.field().getKeys(formulation.domain(), dofField);
  formulation.test().getKeys(formulation.domain(), dofTest);

  // Assemble Systems (tmpA and tmpB) //
  const size_t E = dofField.size();   // Should be equal to dofTest.size().?.

  #pragma omp parallel for
  for(size_t i = 0; i < E; i++)
    SystemAbstract::assemble
      (tmpMat, tmpRHS, i, dofField[i], dofTest[i], term, formulation);
}

void SystemEigen::assemble(void){
  // Enumerate Dofs in DofManager //
  dofM.generateGlobalIdSpace();

  // Alloc Temp Sparse Matrices (not with PETSc) //
  const size_t size = dofM.getUnfixedDofNumber();

  SolverVector<complex<double> > tmpRHS(size);
  SolverMatrix<complex<double> > tmpA(size, size);
  SolverMatrix<complex<double> > tmpB(size, size);

  // Get Formulation Terms //
  formulationPtr term = &FormulationBlock<complex<double> >::weak;

  // Iterate on Formulations A //
  list<const FormulationBlock<complex<double> >*>::iterator it;
  list<const FormulationBlock<complex<double> >*>::iterator end;

  it  = formulation.begin();
  end = formulation.end();

  for(; it != end; it++)
    assembleCom(tmpA, tmpRHS, **it, term);

  // Iterate on Formulations B //
  if(general){
    it  = formulationB.begin();
    end = formulationB.end();

    for(; it != end; it++)
      assembleCom(tmpB, tmpRHS, **it, term);
  }

  // Copy tmpA into Assembled PETSc matrix //
  // Data
  vector<int>              row;
  vector<int>              col;
  vector<complex<double> > value;
  int                      nNZ;

  // Serialize (CStyle) tmpA & Copy
  nNZ = tmpA.serializeCStyle(row, col, value);
  A   = new Mat;

  MatCreateSeqAIJFromTriple(MPI_COMM_SELF, size, size,
                            row.data(), col.data(), value.data(),
                            A, nNZ, PETSC_FALSE);

  // Copy tmpB (CStyle) into Assembled PETSc matrix (if needed) //
  if(general){
    nNZ = tmpB.serializeCStyle(row, col, value);
    B   = new Mat;

    MatCreateSeqAIJFromTriple(MPI_COMM_SELF, size, size,
                              row.data(), col.data(), value.data(),
                              B, nNZ, PETSC_FALSE);
  }

  /*
  tmpA.writeToMatlabFile("Asf_mat.m", "Asf");
  tmpB.writeToMatlabFile("Bsf_mat.m", "Bsf");
  */

  // The SystemEigen is assembled //
  assembled = true;
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
  EPSCreate(MPI_COMM_SELF, &solver);

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
  EPSSetTolerances(solver, 1E-6, 100);

  // Which Eigenpair //
  if(!whichEigenpair.compare("smallest_magnitude")){
    EPSSetWhichEigenpairs(solver, EPS_SMALLEST_MAGNITUDE);
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

  // Get Solution //
  const size_t size = dofM.getUnfixedDofNumber();

  PetscScalar  lambda;
  PetscScalar* x;
  Vec          xPetsc;

  MatGetVecs(*A, PETSC_NULL, &xPetsc);

  EPSGetConverged(solver, &nEigenValues);

  eigenValue  = new fullVector<complex<double> >(nEigenValues);
  eigenVector = new vector<fullVector<complex<double> > >(nEigenValues);

  for(PetscInt i = 0; i < nEigenValues; i++){
    EPSGetEigenpair(solver, i, &lambda, NULL, xPetsc, NULL);

    VecGetArray(xPetsc, &x);

    (*eigenVector)[i].resize(size);
    for(size_t j = 0; j < size; j++)
      (*eigenVector)[i](j) = x[j];

    (*eigenValue)(i) = lambda;
  }

  VecDestroy(&xPetsc);
  EPSDestroy(&solver);

  // System solved ! //
  solved = true;
}

bool SystemEigen::isGeneral(void) const{
  return general;
}

void SystemEigen::getEigenValues(fullVector<complex<double> >& eig) const{
  eig.setAsProxy(*eigenValue, 0, eigenValue->size());
}

void SystemEigen::setWhichEigenpairs(std::string type){
  this->whichEigenpair = type;
}

void SystemEigen::setTarget(Complex target){
  this->target = target.real() + PETSC_i * target.imag();
}

void SystemEigen::
setNumberOfEigenValues(size_t nEigenValues){
  const size_t nDof = dofM.getUnfixedDofNumber();

  if(nEigenValues > nDof)
    throw
      Exception
      ("I can't compute more Eigenvalues (%d) than the number of unknowns (%d)",
       nEigenValues, nDof);

  else
    this->nEigenValues = nEigenValues;
}

size_t SystemEigen::getNComputedSolution(void) const{
  return nEigenValues;
}

void SystemEigen::getSolution(fullVector<std::complex<double> >& sol,
                              size_t nSol) const{
  sol.setAsProxy((*eigenVector)[nSol], 0, (*eigenVector)[nSol].size());
}

void SystemEigen::getSolution(std::map<Dof, std::complex<double> >& sol,
                              size_t nSol) const{
  // Get All Dofs
  map<Dof, complex<double> >::iterator it  = sol.begin();
  map<Dof, complex<double> >::iterator end = sol.end();

  // Loop on Dofs and set Values
  for(; it != end; it++)
    it->second = (*eigenVector)[nSol](dofM.getGlobalId(it->first));
}

void SystemEigen::getSolution(FEMSolution<std::complex<double> >& feSol,
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
  map<Dof, complex<double> > coef;

  for(; it != end; it++)
    coef.insert(pair<Dof, complex<double> >(*it, 0));

  // Iterate on Solutions //
  for(int i = 0; i < nEigenValues; i++){
    // Populate Map
    getSolution(coef, i);

    // FEMSolution
    feSol.addCoefficients(i, 0, domain, fs, coef);
  }
}

void SystemEigen::writeMatrix(string fileName,
                              string matrixName) const{
  throw Exception("SystemEigen::writeMatrix -- not implemented");
}
