#include <cmath>
#include <iostream>

#include "SmallFem.h"

#include "System.h"
#include "SystemHelper.h"

#include "FormulationOO2.h"
#include "FormulationEMDA.h"
#include "FormulationSommerfeld.h"
#include "FormulationSteadyWave.h"
#include "FormulationUpdateEMDA.h"
#include "FormulationUpdateOO2.h"

using namespace std;

Complex fSource(fullVector<double>& xyz){
  return Complex(1, 0);
}

void initMap(FunctionSpace& fs, GroupOfElement& goe, map<Dof, Complex>& data){

  set<Dof> dSet;
  fs.getKeys(goe, dSet);

  set<Dof>::iterator it  = dSet.begin();
  set<Dof>::iterator end = dSet.end();

  for(; it != end; it++)
    data.insert(pair<Dof, Complex>(*it, 0));
}

void displaySolution(map<Dof, Complex>& solution,
                     int id, size_t step, string name){
  int myId;
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  if(myId == id){
    cout << "After Iteration: " << step + 1 << endl;
    map<Dof, Complex>::iterator it  = solution.begin();
    map<Dof, Complex>::iterator end = solution.end();

    for(; it != end; it++)
      cout << name << myId << ": " << it->first.toString()
           << ": " << it->second << endl;
    cout << " --- " << endl;
  }
}

void serialize(const map<Dof, Complex>& data,
               vector<int>& entity, vector<int>& type, vector<Complex>& value){

  map<Dof, Complex>::const_iterator it  = data.begin();
  map<Dof, Complex>::const_iterator end = data.end();

  for(size_t i = 0; it != end; i++, it++){
    entity[i] = (int)(it->first.getEntity());
    type[i]   = (int)(it->first.getType());
    value[i]  = it->second;
  }
}

void unserialize(map<Dof, Complex>& data,
                 const vector<int>& entity, const vector<int>& type,
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

void exchange(int myId,
              vector<int>& outEntity, vector<int>& outType,
              vector<Complex>& outValue,
              vector<int>& inEntity, vector<int>& inType,
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

void compute(const Options& option){
  // MPI //
  int numProcs;
  int myId;
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  if(numProcs != 2)
    throw Exception("I just do two MPI Processes");

  // Get Parameters //
  const string ddmType = option.getValue("-ddm")[1];
  const double k       = atof(option.getValue("-k")[1].c_str());
  const size_t order   = atoi(option.getValue("-o")[1].c_str());
  const size_t maxIt   = atoi(option.getValue("-max")[1].c_str());

  // DDM Formulations //
  const string emdaType("emda");
  const string oo2Type("oo2");

  // Variables
  const double Pi = atan(1.0) * 4;
  double lc       = 0;
  double chi      = 0;
  Complex ooA     = 0;
  Complex ooB     = 0;

  // EMDA Stuff
  if(ddmType == emdaType)
    chi = atof(option.getValue("-chi")[1].c_str());

  // OO2 Stuff
  if(ddmType == oo2Type){
    lc = atof(option.getValue("-lc")[1].c_str());

    double ooXsiMin = 0;
    double ooXsiMax = Pi / lc;
    double ooDeltaK = Pi / .06;

    double tmp0 =
      (k * k - ooXsiMin * ooXsiMin) * (k * k - (k - ooDeltaK) * (k - ooDeltaK));

    double tmp1 =
      (ooXsiMax * ooXsiMax - k * k) * ((k + ooDeltaK) * (k + ooDeltaK) - k * k);

    Complex ooAlpha = pow(Complex(tmp0, 0), 0.25) * Complex(0, 1);
    Complex ooBeta  = pow(Complex(tmp1, 0), 0.25);

    ooA = -(ooAlpha * ooBeta - k * k) / (ooAlpha + ooBeta);
    ooB = Complex(-1, 0) / (ooAlpha + ooBeta);
  }

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume(msh);
  GroupOfElement source(msh);
  GroupOfElement infinity(msh);
  GroupOfElement ddmBorder(msh);

  if(myId == 0){
    volume.add(msh.getFromPhysical(7));
    source.add(msh.getFromPhysical(5));
    infinity.add(msh.getFromPhysical(61));
    ddmBorder.add(msh.getFromPhysical(4));
  }

  else{
    volume.add(msh.getFromPhysical(8));
    //No source//
    infinity.add(msh.getFromPhysical(62));
    ddmBorder.add(msh.getFromPhysical(4));
  }

  // Full Domain //
  GroupOfElement domain(msh);
  domain.add(volume);
  domain.add(source);
  domain.add(infinity);
  domain.add(ddmBorder);

  // Function Space //
  FunctionSpaceScalar fs(domain, order);

  // Steady Wave Formulation //
  FormulationSteadyWave<Complex> wave(volume, fs, k);
  FormulationSommerfeld          sommerfeld(infinity, fs, k);

  // Ddm Formulation Pointers //
  Formulation<Complex>* ddm;
  Formulation<Complex>* upDdm;

  // System Pointers //
  System<Complex>* system;
  System<Complex>* update;

  // Solution Maps //
  map<Dof, Complex> ddmG;
  map<Dof, Complex> solution;

  initMap(fs, ddmBorder, solution);
  initMap(fs, ddmBorder, ddmG);

  // MPI Buffers //
  vector<int>     outEntity(ddmG.size(), 0);
  vector<int>     outType(ddmG.size(), 0);
  vector<Complex> outValue(ddmG.size(), 0);

  vector<int>     inEntity(ddmG.size(), 0);
  vector<int>     inType(ddmG.size(), 0);
  vector<Complex> inValue(ddmG.size(), 0);

  for(size_t step = 0; step < maxIt; step++){
    // Formulations for DDM //
    if(ddmType == emdaType)
      ddm = new FormulationEMDA(ddmBorder, fs, k, chi, ddmG);
    else if(ddmType == oo2Type)
      ddm = new FormulationOO2(ddmBorder, fs, ooA, ooB, ddmG);
    else
      throw Exception("Unknown %s DDM border term", ddmType.c_str());

    // System //
    // Terms
    system = new System<Complex>;
    system->addFormulation(wave);
    system->addFormulation(sommerfeld);
    system->addFormulation(*ddm);

    // Constraint
    if(myId == 0)
      SystemHelper<Complex>::dirichlet(*system, fs, source, fSource);

    // Assemble
    system->assemble();

    // Solve
    system->solve();

    // Get DDM Border Solution //
    system->getSolution(solution, 0);

    try{
      if(option.getValue("-disp").size() > 1)
        displaySolution(solution, atoi(option.getValue("-disp")[1].c_str()),
                        step, "u");
    }
    catch(...){
    }

    // Update G //
    if(ddmType == emdaType)
      upDdm  =
        new FormulationUpdateEMDA(ddmBorder, fs, k, chi, solution, ddmG);

   else if(ddmType == oo2Type)
      upDdm  =
        new FormulationUpdateOO2(ddmBorder, fs, ooA, ooB, solution, ddmG);

    else
      throw Exception("Unknown %s DDM border term", ddmType.c_str());

    update = new System<Complex>;
    update->addFormulation(*upDdm);

    update->assemble();
    update->solve();
    update->getSolution(ddmG, 0);

    // Serialize ddmG
    serialize(ddmG, outEntity, outType, outValue);

    // Exchange
    exchange(myId, outEntity, outType, outValue, inEntity, inType, inValue);

    // ddmG is new ddmG : unserialize
    unserialize(ddmG, inEntity, inType, inValue);

    // Write Solution //
    if(step == maxIt - 1){
      stringstream feSolName;
      FEMSolution<Complex> feSol;
      system->getSolution(feSol, fs, volume);

      feSolName << "ddm" << myId;
      feSol.write(feSolName.str());
    }

    // Clean //
    delete upDdm;
    delete ddm;
    delete system;
    delete update;
  }
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-max,-ddm,-chi,-lc,-disp");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
