#include <iostream>

#include "petscksp.h"

#include "FormulationHelper.h"
#include "SystemHelper.h"
#include "MPIDofMap.h"
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
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myProc);

  // DDM //
  this->wave       = &wave;
  this->sommerfeld = &sommerfeld;
  this->ddm        = &ddm;
  this->upDdm      = &update;
  this->context    = &context;

  // Volume System //
  this->volume     = new System<Complex>;
  this->update     = new System<Complex>;
  this->once       = false;

  // Dirichlet Domain & FunctionSpace //
  this->dirichlet = &dirichlet;
  this->fs        = &context.getFunctionSpace();

  // DDM Dof values //
  this->ddmG = &context.getDDMDofs();
  this->rhs  = &rhs;

  // Build DDM Dof neighbourhood //
  buildNeighbourhood(*this->ddmG);

  // DDM Unknowns number //
  size_t           sizeFull = context.getDDMDofs().size();
  pair<size_t, size_t> size = splitSize();

  // MPI out //
  outEntityOne.resize(size.first, 0);
  outTypeOne.resize(size.first,   0);
  outValueOne.resize(size.first,  0);

  outEntityTwo.resize(size.second, 0);
  outTypeTwo.resize(size.second,   0);
  outValueTwo.resize(size.second,  0);

  // MPI in //
  inEntityOne.resize(size.first, 0);
  inTypeOne.resize(size.first,   0);
  inValueOne.resize(size.first,  0);

  inEntityTwo.resize(size.second, 0);
  inTypeTwo.resize(size.second,   0);
  inValueTwo.resize(size.second,  0);

  // PETSc //
  // Unknown vector
  VecCreate(MPI_COMM_WORLD, &x);
  VecSetSizes(x, sizeFull, PETSC_DECIDE);
  VecSetFromOptions(x);

  // RHS vector
  VecCreate(MPI_COMM_WORLD, &b);
  VecSetSizes(b, sizeFull, PETSC_DECIDE);
  VecSetFromOptions(b);

  // Matrix
  MatCreateShell(MPI_COMM_WORLD, sizeFull, sizeFull, PETSC_DECIDE, PETSC_DECIDE,
                 (void*)(this), &A);

  MatShellSetOperation(A, MATOP_MULT, (void(*)(void))(matMult));
}

SolverDDM::~SolverDDM(void){
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  delete volume;
  delete update;
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
  exchange(*rhs);

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

void SolverDDM::buildNeighbourhood(const map<Dof, Complex>& local){
  // Get owner map //
  multimap<Dof, int> allOwners;
  MPIDofMap<Complex>::getDofOwners(local, allOwners);

  // Init neighbour //
  neighbour.resize(local.size());

  // Look for local Dof neighbour //
  int ownerOne;
  int ownerTwo;
  int myNeighbour;

  map<Dof, Complex>::const_iterator it  = local.begin();
  map<Dof, Complex>::const_iterator end = local.end();

  pair<multimap<Dof, int>::iterator, multimap<Dof, int>::iterator> range;

  for(size_t i = 0; it != end; it++, i++){
    // Look for this Dof owners
    range = allOwners.equal_range(it->first);

    // It shall be owned by only two domains
    if(distance(range.first, range.second) == 2){
      // Get owners of current Dof
      ownerOne = range.first->second;
      range.first++;
      ownerTwo = range.first->second;

      // Get neighbour of current Dof
      if     (ownerOne != myProc && ownerTwo == myProc)
        myNeighbour = ownerOne;

      else if(ownerTwo != myProc && ownerOne == myProc)
        myNeighbour = ownerTwo;

      else
        throw Exception("SolverDDM::buildOwnerMap(): "
                        "no neighbourg found for Dof %s",
                        it->first.toString().c_str());

      // Add it to neighbour
      neighbour[i] =  myNeighbour;
    }

    else
      throw Exception("SolverDDM::buildOwnerMap(): "
                      "a Dof can be shared by only two domains");
  }
}

pair<size_t, size_t> SolverDDM::splitSize(void){
  size_t sizeOne  =  1;
  size_t sizeTwo  =  0;

  neighbourOne = neighbour[0];
  neighbourTwo = -1;

  size_t size = neighbour.size();
  for(size_t i = 1; i < size; i++){
    if(neighbour[i] == neighbourOne){
      sizeOne++;
    }

    else if(neighbour[i] != neighbourOne && neighbourTwo == -1){
      neighbourTwo = neighbour[i];
      sizeTwo++;
    }

    else if(neighbour[i] == neighbourTwo && sizeTwo != 0){
      sizeTwo++;
    }

    else
      throw Exception("SolverDDM::splitDof: unknown neighbour %d",
                      neighbour[i]);
  }

  return pair<size_t, size_t>(sizeOne, sizeTwo);
}

void SolverDDM::serialize(map<Dof, Complex>& data, int neighbour,
                          vector<int>& entity, vector<int>& type,
                          vector<Complex>& value){

  map<Dof, Complex>::const_iterator it  = data.begin();
  map<Dof, Complex>::const_iterator end = data.end();

  for(size_t i = 0, j = 0; it != end; it++, j++){
    if(this->neighbour[j] == neighbour){
      entity[i] = it->first.getEntity();
      type[i]   = it->first.getType();
      value[i]  = it->second;

      i++;
    }
  }
}

void SolverDDM::unserialize(map<Dof, Complex>& data,
                            vector<int>& entity, vector<int>& type,
                            vector<Complex>& value){

  map<Dof, Complex>::iterator end = data.end();
  map<Dof, Complex>::iterator it;


  size_t size = entity.size(); // equal to type.size() and value.size();
  for(size_t i = 0; i < size; i++){
    it = data.find(Dof(entity[i], type[i]));
    if(it != end)
      it->second = value[i];

    else
      throw Exception("SolverDDM::unserialize() unknown Dof %s",
                      Dof(entity[i], type[i]).toString().c_str());
  }
}

template<>
void SolverDDM::exchange(int target, vector<Complex>& out, vector<Complex>& in){
  MPI_Status  status;
  MPI_Request request;
  size_t      size = out.size();

  // Asynchornous exchange //
  MPI_Isend((void*)(out.data()), size, MPI_DOUBLE_COMPLEX,
            target, 0, MPI_COMM_WORLD, &request);

  MPI_Recv((void*)(in.data()), size, MPI_DOUBLE_COMPLEX,
           target, 0, MPI_COMM_WORLD, &status);

  // Buffer Coherence & MPI internals Free //
  MPI_Wait(&request, &status);
}

template<>
void SolverDDM::exchange(int target, vector<int>& out, vector<int>& in){
  MPI_Status  status;
  MPI_Request request;
  size_t      size = out.size();

  // Asynchornous exchange //
  MPI_Isend((void*)(out.data()), size, MPI_INT,
            target, 0, MPI_COMM_WORLD, &request);

  MPI_Recv((void*)(in.data()), size, MPI_INT,
           target, 0, MPI_COMM_WORLD, &status);

  // Buffer Coherence & MPI internals Free //
  MPI_Wait(&request, &status);
}

void SolverDDM::exchange(map<Dof, Complex>& data){
  // Neighbour one //
  if(neighbourOne != -1){
    serialize(data, neighbourOne, outEntityOne, outTypeOne, outValueOne);

    exchange<int>    (neighbourOne, outEntityOne, inEntityOne);
    exchange<int>    (neighbourOne, outTypeOne  , inTypeOne);
    exchange<Complex>(neighbourOne, outValueOne , inValueOne);

    unserialize(data, inEntityOne, inTypeOne, inValueOne);
  }

  // Neighbour two //
  if(neighbourTwo != -1){
    serialize(data, neighbourTwo, outEntityTwo, outTypeTwo, outValueTwo);

    exchange<int>    (neighbourTwo, outEntityTwo, inEntityTwo);
    exchange<int>    (neighbourTwo, outTypeTwo  , inTypeTwo);
    exchange<Complex>(neighbourTwo, outValueTwo , inValueTwo);

    unserialize(data, inEntityTwo, inTypeTwo, inValueTwo);
  }

  PetscBarrier(PETSC_NULL);
}

void SolverDDM::setVecFromDof(Vec& v, map<Dof, Complex>& dof){
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

void SolverDDM::setDofFromVec(Vec& v, map<Dof, Complex>& dof){
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

Complex SolverDDM::fZeroScal(fullVector<double>& xyz){
  return Complex(0, 0);
}

fullVector<Complex> SolverDDM::fZeroVect(fullVector<double>& xyz){
  fullVector<Complex> tmp(3);

  tmp(0) = Complex(0, 0);
  tmp(1) = Complex(0, 0);
  tmp(2) = Complex(0, 0);

  return tmp;
}

PetscErrorCode SolverDDM::matMult(Mat A, Vec x, Vec y){
  // Get SolverDDM Object //
  SolverDDM* solver;
  MatShellGetContext(A, (void**)(&solver));

  // Formulations //
  const Formulation<Complex>& wave       = *solver->wave;
  const Formulation<Complex>& sommerfeld = *solver->sommerfeld;
  Formulation<Complex>&       ddm        = *solver->ddm;
  Formulation<Complex>&       upDdm      = *solver->upDdm;

  // Systems //
  System<Complex>&            volume     = *solver->volume;
  System<Complex>&            update     = *solver->update;

  // Vec x is now the DDM Dof //
  DDMContext&     context = *solver->context;
  map<Dof, Complex>& ddmG = context.getDDMDofs();
  setDofFromVec(x, ddmG);

  // Update ddm Dofs in DDM Context //
  context.setDDMDofs(ddmG);
  ddm.update();

  // Solve Full Volume & Update Problems Once //
  if(!solver->once){
    // Prepare Volume Problem
    volume.addFormulation(wave);
    volume.addFormulation(sommerfeld);
    volume.addFormulation(ddm);

    // Constraint
    const GroupOfElement& dirichlet = *solver->dirichlet;
    const FunctionSpace&  fs        = *solver->fs;

    if(fs.isScalar())
      SystemHelper<Complex>::dirichlet(volume, fs, dirichlet, fZeroScal);
    else
      SystemHelper<Complex>::dirichlet(volume, fs, dirichlet, fZeroVect);

    // Assemble & Solve
    volume.assemble();
    volume.solve();

    // Put new System in DDM Context
    context.setSystem(volume);

    // Update G in Update formulation
    upDdm.update();

    // Update system: add formulation, assemble & solve
    update.addFormulation(upDdm);
    update.assemble();
    update.solve();

    // Once
    solver->once = true;
  }

  // Reassemble Volume & Update RHS and Solve //
  else{
    volume.assembleAgainRHS();  // Volume problem
    volume.solveAgain();        // --------------

    upDdm.update();             // Recompute upDdm terms

    update.assembleAgainRHS();  // Update problem
    update.solveAgain();        // --------------
  }

  // Exchange ddmG //
  update.getSolution(ddmG, 0);
  solver->exchange(ddmG);
  context.setDDMDofs(ddmG);

  // Y = X - ddmG //
  setVecFromDof(y, ddmG);
  VecAYPX(y, -1, x);

  // Wait for MPI coherence and return //
  PetscBarrier(PETSC_NULL);
  PetscFunctionReturn(0);
}
