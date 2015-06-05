#include "SmallFem.h"
#include <iostream>

#include "petscksp.h"

#include "FormulationHelper.h"
#include "SystemHelper.h"
#include "MPIDofMap.h"
#include "Exception.h"
#include "MPIOStream.h"
#include "Timer.h"

#include "SolverDDM.h"

using namespace std;

SolverDDM::SolverDDM(const Formulation<Complex>& wave,
                     const Formulation<Complex>& sommerfeld,
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

  // FunctionSpace //
  this->fs = &context.getFunctionSpace();

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
  VecSetType(x, "mpi");
  VecSetSizes(x, sizeFull, PETSC_DECIDE);
  VecSetFromOptions(x);

  // RHS vector
  VecCreate(MPI_COMM_WORLD, &b);
  VecSetType(b, "mpi");
  VecSetSizes(b, sizeFull, PETSC_DECIDE);
  VecSetFromOptions(b);

  // Matrix
  MatCreateShell(MPI_COMM_WORLD, sizeFull, sizeFull, PETSC_DECIDE, PETSC_DECIDE,
                 (void*)(this), &A);

  // Set opertation //
  MatShellSetOperation(A, MATOP_MULT, (void(*)(void))(SolverDDM::matMult));

  // Iterative solver data //
  useRandGuess = false;
  maxIt        = 100;
  restart      = 100;

  // History //
  hist  = NULL;
  nHist = 0;
}

SolverDDM::~SolverDDM(void){
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  delete volume;
  delete update;

  if(nHist != 0)
    PetscFree(hist);
}

void SolverDDM::setMaximumIteration(int maxIt){
  if(maxIt >= 0)
    this->maxIt = maxIt;
  else
    throw Exception("SolverDDM::setMaximumIteration: "
                    "a negative value was given!");
}

void SolverDDM::setRestart(int restart){
  if(restart >= 0)
    this->restart = restart;
  else
    throw Exception("SolverDDM::setRestart: "
                    "a negative value was given!");
}

void SolverDDM::setRandGuess(void){
  useRandGuess = true;
}

void SolverDDM::unsetRandGuess(void){
  useRandGuess = false;
}

void SolverDDM::solve(void){
  // Create Solver
  KSP solver;

  KSPCreate(MPI_COMM_WORLD, &solver);
  KSPSetOperators(solver, A, A);// For PETSc 3.4 add DIFFERENT_NONZERO_PATTERN);
  KSPSetType(solver, "gmres");

  KSPSetTolerances(solver, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, maxIt);
  KSPGMRESSetRestart(solver, restart);

  // Override with PETSc Database //
  KSPSetFromOptions(solver);

  // Set RHS (b) //
  // Exchange
  exchange(*rhs);

  // Set PETSc Vector and wait for MPI coherence
  setVecFromDof(b, *rhs);
  PetscBarrier(PETSC_NULL);

  // Set History
  PetscMalloc(maxIt * sizeof(PetscReal), (void**)(&hist));
  KSPSetResidualHistory(solver, hist, maxIt, PETSC_TRUE);

  // Randomize initial guess if needed
  if(useRandGuess){
    PetscRandom rand;
    PetscRandomCreate(MPI_COMM_WORLD, &rand);
    PetscRandomSetType(rand, PETSCRAND48);
    PetscRandomSetFromOptions(rand);

    KSPSetInitialGuessNonzero(solver, PETSC_TRUE);
    VecSetRandom(x, rand);
    PetscRandomDestroy(&rand);
  }
  else{
    KSPSetInitialGuessNonzero(solver, PETSC_FALSE);
  }

  // Solve, Get History and Delete Solver //
  KSPSolve(solver, b, x);
  KSPGetResidualHistory(solver, &hist, &nHist);

  PetscBarrier(PETSC_NULL);
  KSPDestroy(&solver);
}

void SolverDDM::constructIterationMatrix(Mat* I){
  // Tell //
  MPIOStream cout(0, std::cout);
  cout << "  " << "SolverDDM::constructIterationMatrix" << endl
       << "  " << "-----------------------------------" << endl;

  // Get Size //
  PetscInt size;
  PetscInt check;
  MatGetSize(A, &size, &check);

  if(size != check)
    throw Exception("SolverDDM::constructIterationMatrix(): "
                    "Shell Matrix is not square");

  // Get Local Size //
  size_t sizeLoc = context->getDDMDofs().size();

  // Allocate Inspection Vector //
  Vec v;
  VecCreate(PETSC_COMM_WORLD, &v);
  VecSetType(v, "mpi");
  VecSetSizes(v, sizeLoc, size);

  // Allocate Resulting Vector //
  Vec w;
  VecDuplicate(v, &w);

  // Get ranges //
  PetscInt start;
  PetscInt stop;
  VecGetOwnershipRange(v, &start, &stop);

  // Sanity //
  PetscInt sanityVecStart;
  PetscInt sanityVecStop;
  VecGetOwnershipRange(w, &sanityVecStart, &sanityVecStop);

  if(start != sanityVecStart || stop != sanityVecStop)
    throw Exception("SolverDDM::constructIterationMatrix(): "
                    "parallel vectors don't have the same parallel pattern");

  // Allocate Iteration Matrix I //
  ///////////////////////////////////////////////////////////////////////////
  //  w0 w1 ... wn                                                         //
  // +--+--+---+--+                                                        //
  // |     p0     |   Each (group of) row is owned by a processor.         //
  // +------------+   It has access to every column of its (group of) row. //
  // |     p1     |                                                        //
  // +------------+   Thus we need PETSc mpiaij format.                    //
  // |     ...    |                                                        //
  // +------------+   Each resulting vector wi is the ith column of I.     //
  // |     pn     |                                                        //
  // +--+--+---+--+                                                        //
  ///////////////////////////////////////////////////////////////////////////
  /*
  MatCreateDense(PETSC_COMM_WORLD, sizeLoc, sizeLoc, size, size, NULL, I);
  */
  MatCreate(PETSC_COMM_WORLD, I);
  MatSetType(*I, "mpiaij");
  MatSetSizes(*I, sizeLoc, sizeLoc, size, size);
  MatSetUp(*I);
  MatMPIAIJSetPreallocation(*I, sizeLoc, PETSC_NULL, size - sizeLoc,PETSC_NULL);

  // Sanity //
  PetscInt sanityMatRowStart;
  PetscInt sanityMatRowStop;
  MatGetOwnershipRange(*I, &sanityMatRowStart, &sanityMatRowStop);

  if(start != sanityMatRowStart || stop != sanityMatRowStop)
    throw Exception("SolverDDM::constructIterationMatrix(): "
                    "parallel vector and matrix don't have "
                    "the same parallel pattern");

  // Local vector value and index //
  PetscScalar*     tmp;
  vector<PetscInt> idx(sizeLoc);

  for(size_t i = 0; i < sizeLoc; i++)
    idx[i] = start + i;

  // Set v = [0..1..0]; apply A; store w in I //
  for(PetscInt i = 0; i < size; i++){
    // Tell
    cout << "  " << "Constructing column: " << i << "/" << size << endl;

    // Set v
    VecScale(v, 0);
    if(start <= i && i < stop)
      VecSetValue(v, i, 1, INSERT_VALUES);

    VecAssemblyBegin(v);
    VecAssemblyEnd(v);

    // Apply A
    matMult(A, v, w);

    // Store in I
    VecGetArray(w, &tmp);
    MatSetValues(*I, sizeLoc, idx.data(), 1, &i, tmp, INSERT_VALUES);
    VecRestoreArray(w, &tmp);
  }

  // Assemble I //
  MatAssemblyBegin(*I, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*I, MAT_FINAL_ASSEMBLY);

  // Clear //
  VecDestroy(&v);
  VecDestroy(&w);
}

void SolverDDM::constructIterationMatrix(string name, string filename){
  // Create Iteration Matrix //
  Mat I;
  constructIterationMatrix(&I);

  // Dump I //
  //////////////////////////////////////////////////////////////
  // Binary format is used                                    //
  //                                                          //
  // Add the petscfolder/bin/matlab path to your matlab paths //
  // and type the following command in matlab.                //
  //                                                          //
  // For real arithmetic:                                     //
  // A = PetscBinaryRead(filename);                           //
  //                                                          //
  // For complex arithmetic:                                  //
  // A = PetscBinaryRead(filename , 'complex', true);         //
  //////////////////////////////////////////////////////////////
  PetscViewer viewer;
  PetscObjectSetName((PetscObject)(I), name.c_str());

  //PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &viewer);
  //PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(),
                        FILE_MODE_WRITE, &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_NATIVE);
  MatView(I, viewer);

  // Clear //
  MatDestroy(&I);
  PetscViewerDestroy(&viewer);
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
                      "Dof %s is shared by more than two domains",
                      it->first.toString().c_str());
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
void SolverDDM::exchange(int target, vector<Complex>& out, vector<Complex>& in,
                         MPI_Request* send, MPI_Request* recv){
  size_t size = out.size();

  // Asynchornous exchange //
  MPI_Isend((void*)(out.data()), size, MPI_DOUBLE_COMPLEX,
            target, 0, MPI_COMM_WORLD, send);

  MPI_Irecv((void*)(in.data()), size, MPI_DOUBLE_COMPLEX,
            target, 0, MPI_COMM_WORLD, recv);
}

template<>
void SolverDDM::exchange(int target, vector<int>& out, vector<int>& in,
                         MPI_Request* send, MPI_Request* recv){
  //MPI_Status status;
  size_t size = out.size();

  // Asynchornous exchange //
  MPI_Isend((void*)(out.data()), size, MPI_INT,
            target, 0, MPI_COMM_WORLD, send);

  MPI_Irecv((void*)(in.data()), size, MPI_INT,
            target, 0, MPI_COMM_WORLD, recv);
}

void SolverDDM::exchange(map<Dof, Complex>& data){
  Timer exch;
  Timer ser;
  Timer user;

  MPI_Request  req[12];
  MPI_Status  stat[12];

  // Exchange with neighbour one //
  if(neighbourOne != -1){
    ser.start();
    serialize(data, neighbourOne, outEntityOne, outTypeOne, outValueOne);
    ser.stop();

    exch.start();
    exchange<int>    (neighbourOne, outEntityOne, inEntityOne, &req[0],&req[1]);
    exchange<int>    (neighbourOne, outTypeOne  , inTypeOne  , &req[2],&req[3]);
    exchange<Complex>(neighbourOne, outValueOne , inValueOne , &req[4],&req[5]);
    exch.stop();
  }

  // Exchange with neighbour two //
  if(neighbourTwo != -1){
    ser.start();
    serialize(data, neighbourTwo, outEntityTwo, outTypeTwo, outValueTwo);
    ser.stop();

    exch.start();
    exchange<int>    (neighbourTwo, outEntityTwo, inEntityTwo, &req[6],&req[7]);
    exchange<int>    (neighbourTwo, outTypeTwo  , inTypeTwo  , &req[8],&req[9]);
    exchange<Complex>(neighbourTwo, outValueTwo , inValueTwo ,
                      &req[10],&req[11]);
    exch.stop();
  }

  // Reconstruct data from neighbour one //
  if(neighbourOne != -1){
    MPI_Waitall(6, &req[0], &stat[0]);
    user.start();
    unserialize(data, inEntityOne, inTypeOne, inValueOne);
    user.stop();
  }

  // Reconstruct data from neighbour two //
  if(neighbourTwo != -1){
    MPI_Waitall(6, &req[6], &stat[6]);
    user.start();
    unserialize(data, inEntityTwo, inTypeTwo, inValueTwo);
    user.stop();
  }
  /*
  cout << "Exchange    (" << myProc << "): " << exch.time() << endl << flush;
  cout << "Serialize   (" << myProc << "): " <<  ser.time() << endl << flush;
  cout << "Unserialize (" << myProc << "): " << user.time() << endl << flush;
  */
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
  // Timer //
  Timer all;
  Timer sys;

  all.start();
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
  sys.start();
  if(!solver->once){
    // Prepare Volume Problem
    volume.addFormulation(wave);
    volume.addFormulation(sommerfeld);
    volume.addFormulation(ddm);

    // Function Space //
    const FunctionSpace& fs = *solver->fs;

    // Constraint
    const vector<const GroupOfElement*>&
      dirichlet = context.getDirichletDomain();

    size_t sizeDirichlet = dirichlet.size();
    if(fs.isScalar())
      for(size_t i = 0; i < sizeDirichlet; i++)
        SystemHelper<Complex>::dirichlet(volume, fs, *dirichlet[i], fZeroScal);
    else
      for(size_t i = 0; i < sizeDirichlet; i++)
        SystemHelper<Complex>::dirichlet(volume, fs, *dirichlet[i], fZeroVect);

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
  sys.stop();

  update.getSolution(ddmG, 0);
  solver->exchange(ddmG);
  context.setDDMDofs(ddmG);

  // Y = X - ddmG //
  setVecFromDof(y, ddmG);
  VecAYPX(y, -1, x);

  // Wait for MPI coherence and return //
  PetscBarrier(PETSC_NULL);
  all.stop();

  int myProc;
  MPI_Comm_rank(MPI_COMM_WORLD,&myProc);
  cout << "Total  (" << myProc << "): " << all.time() << endl << flush;
  cout << "System (" << myProc << "): " << sys.time() << endl << flush;

  PetscFunctionReturn(0);
}

void SolverDDM::getHistory(vector<double>& history) const{
  history.resize(nHist);

  for(PetscInt i = 0; i < nHist; i++)
    history[i] = hist[i];
}
