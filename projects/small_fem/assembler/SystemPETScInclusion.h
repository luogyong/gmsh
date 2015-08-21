///////////////////////////////////////////
// Templates Implementations for SystemPETSc: //
// Inclusion compilation model           //
//                                       //
// Damn you gcc: we want 'export' !      //
///////////////////////////////////////////

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "petscksp.h"

template<typename scalar>
SystemPETSc<scalar>::SystemPETSc(void){
  // Init //
  A      = NULL;
  b      = NULL;
  xPetsc = NULL;
  x      = NULL;

  // Dof Manager //
  this->dofM = new DofManager<scalar>;

  // Per-thread non-zero term //
  #pragma omp parallel
  {
    #pragma omp master
    this->nNZCount.resize(omp_get_num_threads());
  }

  for(size_t i = 0; i < this->nNZCount.size(); i++)
    this->nNZCount[i] = 0;

  // The system is not assembled and not solved //
  this->assembled = false;
  this->solved    = false;
}

template<typename scalar>
SystemPETSc<scalar>::~SystemPETSc(void){
  delete this->dofM;

  MatDestroy(&A);
  VecDestroy(&b);
  VecDestroy(&xPetsc);

  if(x)
    delete x;
}

template<typename scalar>
void SystemPETSc<scalar>::assemble(void){
  // Enumerate Dofs in DofManager //
  this->dofM->generateGlobalIdSpace();

  // Formulations Iterators //
  typename std::list<const FormulationBlock<scalar>*>::iterator it;
  typename std::list<const FormulationBlock<scalar>*>::iterator end;

  // Count FE terms //
  it  = this->formulation.begin();
  end = this->formulation.end();

  for(; it != end; it++){
    // Get All Dofs (Field & Test) per Element
    const std::vector<std::vector<Dof> >& dofField =
      (*it)->field().getKeys((*it)->domain());
    const std::vector<std::vector<Dof> >& dofTest  =
      (*it)->test().getKeys((*it)->domain());

    // Count
    const size_t E = dofField.size(); // Should be equal to dofTest.size().?.

    #pragma omp parallel for
    for(size_t i = 0; i < E; i++)
      this->nNZCount[omp_get_thread_num()] =
        SystemAbstract<scalar>::countTerms(this->nNZCount[omp_get_thread_num()],
                                           i, dofField[i], dofTest[i], **it);
  }

  // Alloc //
  const size_t size = this->dofM->getLocalSize();

  SolverMatrix<scalar>* tmpA;
  SolverVector<scalar>* tmpb;

  tmpA = new SolverMatrix<scalar>(size, size, this->nNZCount);
  tmpb = new SolverVector<scalar>(size);

  // Assemble //
  it  = this->formulation.begin();
  end = this->formulation.end();

  for(; it != end; it++){
    // Get All Dofs (Field & Test) per Element
    const std::vector<std::vector<Dof> >& dofField =
      (*it)->field().getKeys((*it)->domain());
    const std::vector<std::vector<Dof> >& dofTest  =
      (*it)->test().getKeys((*it)->domain());

    // Assemble
    const size_t E = dofField.size(); // Should be equal to dofTest.size().?.

    #pragma omp parallel for
    for(size_t i = 0; i < E; i++)
      SystemAbstract<scalar>::
        assemble(*tmpA, *tmpb, i, dofField[i], dofTest[i], **it);
  }

  // To Petsc //
  // Matrix
  int*      row;
  int*      col;
  scalar* value;
  int nNZ = tmpA->get(&row, &col, &value);

  PetscScalar* petsc = toPetscScalar(value, nNZ);
  MatCreateSeqAIJFromTriple(PETSC_COMM_SELF, size, size, row, col, petsc, &A,
                            nNZ, PETSC_TRUE);

  // RHS
  VecCreate(PETSC_COMM_SELF, &b);
  VecSetSizes(b, size, size);
  VecSetType(b, "seq");

  scalar* rhs = tmpb->getData();
  for(size_t i = 0; i < size; i++)
    VecSetValue(b, i, rhs[i], INSERT_VALUES);

  // X
  VecCreate(PETSC_COMM_SELF, &xPetsc);
  VecSetSizes(xPetsc, size, size);
  VecSetType(xPetsc, "seq");

  // Assemble
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(b);

  // Clear //
  delete   tmpA;
  delete   tmpb;
  delete[] petsc;

  // The system is assembled //
  this->assembled = true;
}

template<typename scalar>
void SystemPETSc<scalar>::solve(void){
  // Is the System assembled ? //
  if(!this->assembled)
    assemble();

  // Max It //
  const int maxIt = this->dofM->getLocalSize();

  // Use Petsc //
  KSP solver;
  PC  precond;

  // Create solver
  KSPCreate(MPI_COMM_WORLD, &solver);
  KSPSetOperators(solver, A, A); // For PETSc 3.4 add DIFFERENT_NONZERO_PATTERN)
  KSPGetPC(solver, &precond);

  // GMRES
  KSPSetType(solver, "gmres");
  KSPSetTolerances(solver, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, maxIt);
  KSPGMRESSetRestart(solver, maxIt); // No restart
  PCSetType(precond, PCNONE); // No preconditioner

  // MUMPS
  //PCSetType(precond, PCLU);
  //PCFactorSetMatSolverPackage(precond, MATSOLVERMUMPS);

  // Override with PETSc Database
  KSPSetFromOptions(solver);
  PCSetFromOptions(precond);

  // Solve and delete
  KSPSolve(solver, b, xPetsc);
  KSPDestroy(&solver);

  // Get Solution //
  getSolution();

  // System solved ! //
  this->solved = true;
}

template<typename scalar>
size_t SystemPETSc<scalar>::getNComputedSolution(void) const{
  return 1;
}

template<typename scalar>
void SystemPETSc<scalar>::
getSolution(fullVector<scalar>& sol, size_t nSol) const{
  sol.setAsProxy(*x, 0, x->size());
}

template<typename scalar>
void SystemPETSc<scalar>::
getSolution(std::map<Dof, scalar>& sol, size_t nSol) const{
  // Get All Dofs
  typename std::map<Dof, scalar>::iterator it  = sol.begin();
  typename std::map<Dof, scalar>::iterator end = sol.end();

  // Loop on Dofs and set Values
  for(; it != end; it++){
    size_t gId = this->dofM->getGlobalId(it->first);

    if(gId == DofManager<scalar>::isFixedId())
      it->second = this->dofM->getValue(it->first);

    else
      it->second = (*x)(gId);
  }
}

template<typename scalar>
void SystemPETSc<scalar>::getSolution(FEMSolution<scalar>& feSol,
                                      const FunctionSpace& fs,
                                      const GroupOfElement& domain) const{
  // Solved ?
  if(!this->solved)
    throw Exception("SystemPETSc: addSolution -- System not solved");

  // Coefficients //
  // Get Dofs
  std::set<Dof> dof;
  fs.getKeys(domain, dof);

  // Get Coefficient
  const std::set<Dof>::iterator  end = dof.end();
  std::set<Dof>::iterator        it  = dof.begin();
  typename std::map<Dof, scalar> coef;

  for(; it != end; it++)
    coef.insert(std::pair<Dof, scalar>(*it, 0));

  // Populate Map
  getSolution(coef, 0);

  // FEMSolution //
  feSol.addCoefficients(0, 0, domain, fs, coef);
}

template<typename scalar>
void SystemPETSc<scalar>::
getSolution(FEMSolution<scalar>& feSol,
            const FunctionSpace& fs,
            const std::vector<const GroupOfElement*>& domain) const{
  // Get size
  const size_t size = domain.size();

  // Get solution for each domain
  for(size_t i = 0; i < size; i++)
    getSolution(feSol, fs, *domain[i]);
}

template<typename scalar>
void SystemPETSc<scalar>::writeMatrix(std::string fileName,
                                      std::string matrixName) const{
  // Names //
  std::string name(matrixName);
  std::string file(fileName);

  file.append(std::string(".m"));
  PetscObjectSetName((PetscObject)(A), name.c_str());

  // Viewers //
  PetscViewer viewer;

  // Binary
  // PetscViewerBinaryOpen(PETSC_COMM_WORLD,file.c_str(),FILE_MODE_WRITE,&viewer);
  // PetscViewerSetFormat(viewer, PETSC_VIEWER_NATIVE);

  // ASCII
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, file.c_str(), &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

  // Do your job //
  MatView(A, viewer);

  // Clean & coherence //
  PetscViewerDestroy(&viewer);
}
