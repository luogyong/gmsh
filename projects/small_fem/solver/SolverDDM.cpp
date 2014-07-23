#include <iostream>

#include "petscksp.h"

#include "FormulationHelper.h"
#include "SystemHelper.h"
#include "Exception.h"
#include "SolverDDM.h"

using namespace std;

SolverDDM::SolverDDM(const Formulation<Complex>& wave,
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
  this->fullCtx.once       = false;
  this->fullCtx.myId       = this->myId;
  this->fullCtx.DDMctx     = this->context;
  this->fullCtx.dirichlet  = this->dirichlet;
  this->fullCtx.fs         = this->fs;
  this->fullCtx.wave       = this->wave;
  this->fullCtx.sommerfeld = this->sommerfeld;
  this->fullCtx.ddm        = this->ddm;
  this->fullCtx.upDdm      = this->upDdm;
  this->fullCtx.volume     = new System<Complex>;
  this->fullCtx.update     = NULL; // new System<Complex>;
  this->fullCtx.outValue   = &this->outValue;
  this->fullCtx.inValue    = &this->inValue;

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

SolverDDM::~SolverDDM(void){
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  delete this->fullCtx.volume;
  // delete this->fullCtx.update;
}

void SolverDDM::solve(int nStep){
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

  // Set PETSc Vector and wait for MPI coherence
  setVecFromDof(b, *rhs);
  PetscBarrier(PETSC_NULL);

  // Solve and Delete Solver //
  KSPSolve(solver, b, x);
  PetscBarrier(PETSC_NULL);

  KSPDestroy(&solver);
}

void SolverDDM::getSolution(map<Dof, Complex>& ddm){
  setDofFromVec(x, ddm);
}

PetscErrorCode SolverDDM::matMult(Mat A, Vec x, Vec y){
  // Get FullContext //
  FullContext* fullCtx;
  MatShellGetContext(A, (void**)(&fullCtx));

  // Formulations //
  const Formulation<Complex>& wave       = *fullCtx->wave;
  const Formulation<Complex>& sommerfeld = *fullCtx->sommerfeld;
  Formulation<Complex>&       ddm        = *fullCtx->ddm;
  Formulation<Complex>&       upDdm      = *fullCtx->upDdm;

  System<Complex>&            volume     = *fullCtx->volume;
  // System<Complex>&            update     = *fullCtx->update;

  // Vec x is now the DDM Dof //
  DDMContext&     context = *fullCtx->DDMctx;
  map<Dof, Complex>& ddmG = context.getDDMDofs();
  setDofFromVec(x, ddmG);

  // Update ddm Dofs in DDM Context //
  context.setDDMDofs(ddmG);
  ddm.update();

  // Solve Full Volume Problem & Prepare Update Problem (once) //
  if(!fullCtx->once){
    // Prepare Volume Problem
    volume.addFormulation(wave);
    volume.addFormulation(sommerfeld);
    volume.addFormulation(ddm);

    // Constraint
    const GroupOfElement& dirichlet = *fullCtx->dirichlet;
    const FunctionSpace&  fs        = *fullCtx->fs;

    SystemHelper<Complex>::dirichlet(volume, fs, dirichlet, fZero);

    // Assemble & Solve
    volume.assemble();
    volume.solve();

    // Put new System in DDM Context
    context.setSystem(volume);

    // Prepare DDM Update Formulation
    // update.addFormulation(upDdm);

    // Once
    fullCtx->once = true;
  }

  // Reassemble Volume RHS and Solve //
  else{
    volume.assembleAgainRHS();
    volume.solve();
  }

  // Update G & Solve Update Problem //
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

  // Y = X - ddmG //
  setVecFromDof(y, ddmG);
  VecAYPX(y, -1, x);

  // Wait for MPI coherence and return //
  PetscBarrier(PETSC_NULL);
  PetscFunctionReturn(0);
}

void SolverDDM::setVecFromDof(Vec& v, std::map<Dof, Complex>& dof){
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

void SolverDDM::setDofFromVec(Vec& v, std::map<Dof, Complex>& dof){
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

Complex SolverDDM::fZero(fullVector<double>& xyz){
  return Complex(0, 0);
}

void SolverDDM::serialize(const map<Dof, Complex>& data,
                          vector<Complex>& value){

  map<Dof, Complex>::const_iterator it  = data.begin();
  map<Dof, Complex>::const_iterator end = data.end();

  for(size_t i = 0; it != end; i++, it++)
    value[i]  = it->second;
}

void SolverDDM::unserialize(map<Dof, Complex>& data,
                            const vector<Complex>& value){

  map<Dof, Complex>::iterator it  = data.begin();
  map<Dof, Complex>::iterator end = data.end();

  for(size_t i = 0; it != end; it++, i++)
    it->second = value[i];
}

void SolverDDM::exchange(int myId,
                         vector<Complex>& outValue,
                         vector<Complex>& inValue){
  MPI_Status  status;
  MPI_Request request;
  size_t      size  = outValue.size();
  int         nxtId = (myId + 1) % 2;

  // Asynchornous exchange //
  MPI_Isend((void*)(outValue.data()), size, MPI_DOUBLE_COMPLEX,
            nxtId, 0, MPI_COMM_WORLD, &request);

  MPI_Recv((void*)(inValue.data()), size, MPI_DOUBLE_COMPLEX,
           nxtId, 0, MPI_COMM_WORLD, &status);

  // Buffer Coherence & MPI internals Free //
  MPI_Wait(&request, &status);
}
