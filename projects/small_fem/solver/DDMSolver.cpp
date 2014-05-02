#include <iostream>

#include "FormulationHelper.h"
#include "SystemHelper.h"
#include "Exception.h"
#include "DDMSolver.h"

using namespace std;

DDMSolver::DDMSolver(const Formulation<Complex>& wave,
                     const Formulation<Complex>& sommerfeld,
                     DDMContext& context,
                     Formulation<Complex>& ddm,
                     Formulation<Complex>& update,
                     map<Dof, Complex>& ddmG,
                     const GroupOfElement& source,
                     Complex (*fSource)(fullVector<double>& xyz)){
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

  // FunctionSpace & Domain //
  this->fs        = &context.getFunctionSpace();
  this->ddmBorder = &context.getDomain();
  this->source    = &source;

  // fSource //
  this->fSource = fSource;

  // DDM Dof values //
  this->ddmG = &ddmG;

  // Systems //
  this->volume = NULL;
  this->update = NULL;

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
}

DDMSolver::~DDMSolver(void){
}

void DDMSolver::solve(int nStep){
  for(int step = 0; step < nStep; step++){
    // Volume //
    // Terms
    volume = new System<Complex>;
    volume->addFormulation(*wave);
    volume->addFormulation(*sommerfeld);
    volume->addFormulation(*ddm);

    // Constraint
    if(myId == 0)
      SystemHelper<Complex>::dirichlet(*volume, *fs, *source, fSource);

    // Assemble
    volume->assemble();

    // Solve
    volume->solve();

    // Put new System in DDM Context //
    context->setSystem(*volume);

    // Display //
    displaySolution(*volume, *fs, *ddmBorder, 0, step, "u");

    // Update G //
    upDdm->update(); // update volume solution (at DDM border)

    update = new System<Complex>;
    update->addFormulation(*upDdm);

    update->assemble();
    update->solve();
    update->getSolution(*ddmG, 0);

    // Serialize ddmG
    serialize(*ddmG, outEntity, outType, outValue);

    // Exchange
    exchange(myId, outEntity, outType, outValue, inEntity, inType, inValue);

    // ddmG is new ddmG : unserialize
    unserialize(*ddmG, inEntity, inType, inValue);

    // Get new DDM Dofs //
    context->setDDMDofs(*ddmG);

    // Update DDM Formulations //
    ddm->update();

    // Write Solution //
    /*
    if(step == maxIt - 1){
      stringstream feSolName;
      FEMSolution<Complex> feSol;
      system->getSolution(feSol, fs, volume);

      feSolName << "ddm" << myId;
      feSol.write(feSolName.str());
    }
    */
    // Clean //
    delete volume;
    delete update;
  }
}

void DDMSolver::serialize(const map<Dof, Complex>& data,
                          vector<int>& entity,
                          vector<int>& type,
                          vector<Complex>& value){

  map<Dof, Complex>::const_iterator it  = data.begin();
  map<Dof, Complex>::const_iterator end = data.end();

  for(size_t i = 0; it != end; i++, it++){
    entity[i] = (int)(it->first.getEntity());
    type[i]   = (int)(it->first.getType());
    value[i]  = it->second;
  }
}

void DDMSolver::unserialize(map<Dof, Complex>& data,
                            const vector<int>& entity,
                            const vector<int>& type,
                            const vector<Complex>& value){

  map<Dof, Complex>::iterator it  = data.begin();
  map<Dof, Complex>::iterator end = data.end();

  for(size_t i = 0; it != end; it++, i++){
    if((int)(it->first.getType()) != type[i])
      throw Exception("Snif");

    if((int)(it->first.getEntity()) != entity[i])
      throw Exception("Snif");

    it->second = value[i];
  }
}

void DDMSolver::exchange(int myId,
                         vector<int>& outEntity,
                         vector<int>& outType,
                         vector<Complex>& outValue,
                         vector<int>& inEntity,
                         vector<int>& inType,
                         vector<Complex>& inValue){
  MPI_Status s;
  size_t     size = outEntity.size();

  if(myId == 0){
    // Send to 1
    MPI_Ssend((void*)(outEntity.data()), size, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Ssend((void*)(outType.data()),   size, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Ssend((void*)(outValue.data()),  size,
              MPI::DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD);

    // Recv from 1
    MPI_Recv((void*)(inEntity.data()), size, MPI_INT, 1, 0, MPI_COMM_WORLD, &s);
    MPI_Recv((void*)(inType.data()),   size, MPI_INT, 1, 0, MPI_COMM_WORLD, &s);
    MPI_Recv((void*)(inValue.data()),  size,
             MPI::DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD, &s);
  }

  if(myId == 1){
    // Recv from 0
    MPI_Recv((void*)(inEntity.data()), size, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
    MPI_Recv((void*)(inType.data()),   size, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
    MPI_Recv((void*)(inValue.data()),  size,
             MPI::DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, &s);

    // Send to 0
    MPI_Ssend((void*)(outEntity.data()), size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Ssend((void*)(outType.data()),   size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Ssend((void*)(outValue.data()),  size,
              MPI::DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
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
