#include <iostream>

#include "petscksp.h"

#include "FormulationHelper.h"
#include "SystemHelper.h"
#include "Exception.h"
#include "DDMSolver.h"

using namespace std;

DDMSolver::DDMSolver(const Formulation<Complex>& wave,
                     const Formulation<Complex>& sommerfeld,
                     const GroupOfElement& dirichlet,
                     DDMContext& context,
                     Formulation<Complex>& ddm,
                     Formulation<Complex>& update,
                     map<Dof, Complex>& rhs){
  // MPI //
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  if(numProcs != 2)
    throw Exception("I just do two MPI Processes");

  // DDM //
  this->wave       = &wave;
  this->sommerfeld = &sommerfeld;
  this->context    = &context;
  this->ddm        = &ddm;
  this->upDdm      = &update;

  // Dirichlet Domain & FunctionSpace //
  this->dirichlet = &dirichlet;
  this->fs        = &context.getFunctionSpace();

  // DDM Dof values //
  this->ddmG = &context.getDDMDofs();
  this->rhs  = &rhs;


  // DDM Unknowns number //
  size_t size = context.getDDMDofs().size();

  // MPI out //
  outEntity.resize(size, 0);
  outType.resize(size, 0);
  outValue.resize(size, 0);

  // MPI in //
  inEntity.resize(size, 0);
  inType.resize(size, 0);
  inValue.resize(size, 0);

  // FullContext //
  this->fullCtx.myId       = this->myId;
  this->fullCtx.DDMctx     = this->context;
  this->fullCtx.dirichlet  = this->dirichlet;
  this->fullCtx.fs         = this->fs;
  this->fullCtx.wave       = this->wave;
  this->fullCtx.sommerfeld = this->sommerfeld;
  this->fullCtx.ddm        = this->ddm;
  this->fullCtx.upDdm      = this->upDdm;
  this->fullCtx.outValue   = &this->outValue;
  this->fullCtx.inValue   = &this->inValue;

  // PETSc //
  // Unknown vector
  VecCreate(MPI_COMM_WORLD, &x);
  VecSetSizes(x, size, PETSC_DECIDE);
  VecSetFromOptions(x);

  // RHS vector
  VecCreate(MPI_COMM_WORLD, &b);
  VecSetSizes(b, size, PETSC_DECIDE);
  VecSetFromOptions(b);

  // Matrix
  MatCreateShell(MPI_COMM_WORLD, size, size, PETSC_DECIDE, PETSC_DECIDE,
                 (void*)(&this->fullCtx), &A);

  MatShellSetOperation(A, MATOP_MULT, (void(*)(void))(matMult));
}

DDMSolver::~DDMSolver(void){
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);
}

void DDMSolver::solve(int nStep){
  // Create Solver
  KSP solver;

  KSPCreate(MPI_COMM_WORLD, &solver);
  KSPSetOperators(solver, A, A, DIFFERENT_NONZERO_PATTERN);
  KSPSetType(solver, "gmres");

  KSPSetTolerances(solver, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, nStep);

  // Override with PETSc Database //
  KSPSetFromOptions(solver);

  // Set RHS (b) //

  // Exchange

  serialize(*rhs, outValue);
  exchange(myId, outValue, inValue);
  unserialize(*rhs, inValue);
  PetscBarrier(PETSC_NULL);

  setVecFromDof(b, *rhs);
  //VecScale(b, -1);

  //VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  // Solve and Delete Solver //
  PetscBarrier(PETSC_NULL);
  KSPSolve(solver, b, x);
  PetscBarrier(PETSC_NULL);
  KSPDestroy(&solver);
}

void DDMSolver::getSolution(map<Dof, Complex>& ddm){
  setDofFromVec(x, ddm);

  //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
}

void DDMSolver::see(const map<Dof, Complex>& dof, int id, size_t step, string name){
  int myId;
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  if(myId == id){
    cout << "After Iteration: " << step + 1 << endl;
    map<Dof, Complex>::const_iterator it  = dof.begin();
    map<Dof, Complex>::const_iterator end = dof.end();

    for(; it != end; it++)
      cout << name << myId << ": " << it->first.toString()
           << ": " << it->second << endl;
    cout << " --- " << endl;
  }
}

PetscErrorCode DDMSolver::matMult(Mat A, Vec x, Vec y){
  //VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  // Get FullContext //
  FullContext* fullCtx;
  MatShellGetContext(A, (void**)(&fullCtx));

  // Formulations //
  const Formulation<Complex>& wave       = *fullCtx->wave;
  const Formulation<Complex>& sommerfeld = *fullCtx->sommerfeld;
  Formulation<Complex>&       ddm        = *fullCtx->ddm;
  Formulation<Complex>&       upDdm      = *fullCtx->upDdm;

  // Vec x is now the DDM Dof //
  DDMContext&     context = *fullCtx->DDMctx;
  map<Dof, Complex>& ddmG = context.getDDMDofs();
  setDofFromVec(x, ddmG);

  // Exchange ddmG //
  //serialize(ddmG, *fullCtx->outValue);
  //exchange(fullCtx->myId, *fullCtx->outValue, *fullCtx->inValue);
  //unserialize(ddmG, *fullCtx->inValue);

  // Update ddm Dofs in DDM Context //
  context.setDDMDofs(ddmG);
  ddm.update();

  // Volume Problem //
  System<Complex> volume;
  volume.addFormulation(wave);
  volume.addFormulation(sommerfeld);
  volume.addFormulation(ddm);

  // Constraint
  int                   myId      =  fullCtx->myId;
  const GroupOfElement& dirichlet = *fullCtx->dirichlet;
  const FunctionSpace&  fs        = *fullCtx->fs;

  if(myId == 0)
    SystemHelper<Complex>::dirichlet(volume, fs, dirichlet, fZero);

  // Assemble & Solve
  volume.assemble();
  volume.solve();

  // Put new System in DDM Context //
  context.setSystem(volume);

  // Update G //
  upDdm.update(); // update volume solution (at DDM border)

  System<Complex> update;
  update.addFormulation(upDdm);

  update.assemble();
  update.solve();
  update.getSolution(ddmG, 0);

  // Exchange ddmG //
  serialize(ddmG, *fullCtx->outValue);
  exchange(fullCtx->myId, *fullCtx->outValue, *fullCtx->inValue);
  unserialize(ddmG, *fullCtx->inValue);

  context.setDDMDofs(ddmG);

  // Exchange Vec x //
  //serialize(x, *fullCtx->outValue);
  //exchange(fullCtx->myId, *fullCtx->outValue, *fullCtx->inValue);
  //unserialize(x, *fullCtx->inValue);

  PetscBarrier(PETSC_NULL);
  // Y = X - ddmG //
  //see(ddmG, 1, 0, "g");
  setVecFromDof(y, ddmG);
  //VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  VecAYPX(y, -1, x);

  // Done //
  PetscBarrier(PETSC_NULL);

  PetscFunctionReturn(0);
}

void DDMSolver::setVecFromDof(Vec& v, std::map<Dof, Complex>& dof){
  // Pointer to PETSc Vec data //
  Complex* ptr;
  VecGetArray(v, &ptr);

  // Serialize Dof values into PETSc Vec //
  map<Dof, Complex>::iterator it  = dof.begin();
  map<Dof, Complex>::iterator end = dof.end();

  for(size_t i = 0; it != end; it++, i++)
    ptr[i] = it->second;

  // PETSc Coherence //
  VecRestoreArray(v, &ptr);
}

void DDMSolver::setDofFromVec(Vec& v, std::map<Dof, Complex>& dof){
  // Pointer to PETSc Vec data //
  Complex* ptr;
  VecGetArray(v, &ptr);

  // Copy PETSc Vec values into Map //
  map<Dof, Complex>::iterator it  = dof.begin();
  map<Dof, Complex>::iterator end = dof.end();

  for(size_t i = 0; it != end; it++, i++)
    it->second = ptr[i];

  // PETSc Coherence //
  VecRestoreArray(v, &ptr);
}

Complex DDMSolver::fZero(fullVector<double>& xyz){
  return Complex(0, 0);
}

void DDMSolver::see(Vec& v){
  int myId;
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  Complex* ptr;
  VecGetArray(v, &ptr);

  int size;
  VecGetLocalSize(v, &size);

  if(myId == 0){
    cout << "---" << myId << "---" << endl;
    for(int i = 0; i < size; i++)
      cout << ptr[i] << endl;
    cout << "---" << myId << "---" << endl;
  }

  VecRestoreArray(v, &ptr);
}

void DDMSolver::serialize(const map<Dof, Complex>& data,
                          vector<Complex>& value){

  map<Dof, Complex>::const_iterator it  = data.begin();
  map<Dof, Complex>::const_iterator end = data.end();

  for(size_t i = 0; it != end; i++, it++)
    value[i]  = it->second;
}

void DDMSolver::serialize(const Vec& x,
                          vector<Complex>& value){
  // Pointer to PETSc Vec data //
  Complex* ptr;
  VecGetArray(x, &ptr);

  // Size //
  int size;
  VecGetLocalSize(x, &size);

  // Copy //
  for(int i = 0; i < size; i++)
    value[i]  = ptr[i];

  // PETSc Coherence //
  VecRestoreArray(x, &ptr);
}

void DDMSolver::unserialize(map<Dof, Complex>& data,
                            const vector<Complex>& value){

  map<Dof, Complex>::iterator it  = data.begin();
  map<Dof, Complex>::iterator end = data.end();

  for(size_t i = 0; it != end; it++, i++)
    it->second = value[i];
}

void DDMSolver::unserialize(Vec& x,
                            const vector<Complex>& value){
  // Pointer to PETSc Vec data //
  Complex* ptr;
  VecGetArray(x, &ptr);

  // Size //
  int size;
  VecGetLocalSize(x, &size);

  // Copy //
  for(int i = 0; i < size; i++)
    ptr[i]  = value[i];

  // PETSc Coherence //
  VecRestoreArray(x, &ptr);
}

void DDMSolver::exchange(int myId,
                         vector<Complex>& outValue,
                         vector<Complex>& inValue){
  MPI_Status s;
  size_t     size = outValue.size();

  if(myId == 0){
    // Send to 1
    MPI_Ssend((void*)(outValue.data()),  size,
              MPI_DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD);

    // Recv from 1
    MPI_Recv((void*)(inValue.data()),  size,
             MPI_DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD, &s);
  }

  if(myId == 1){
    // Recv from 0
    MPI_Recv((void*)(inValue.data()),  size,
             MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, &s);

    // Send to 0
    MPI_Ssend((void*)(outValue.data()),  size,
              MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
  }
}

void DDMSolver::displaySolution(const System<Complex>& system,
                                const FunctionSpace& fs,
                                const GroupOfElement& goe,
                                int id, int step, string name){
  int myId;
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  if(myId == id){
    // Get unknonws
    map<Dof, Complex> solution;
    FormulationHelper::initDofMap(fs, goe, solution);

    // Get Solution
    system.getSolution(solution, 0);

    // Print it
    cout << "After Iteration: " << step + 1 << endl;
    map<Dof, Complex>::iterator it  = solution.begin();
    map<Dof, Complex>::iterator end = solution.end();

    for(; it != end; it++)
      cout << name << myId << ": " << it->first.toString()
           << ": " << it->second << endl;
    cout << " --- " << endl;
  }
}
