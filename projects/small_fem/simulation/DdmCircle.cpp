#include <cmath>
#include <iostream>

#include "SmallFem.h"

#include "SolverDDM.h"

#include "DDMContextEMDA.h"
#include "DDMContextOO2.h"
#include "DDMContextOSRC.h"

#include "System.h"
#include "SystemHelper.h"

#include "FormulationHelper.h"

#include "FormulationOO2.h"
#include "FormulationEMDA.h"
#include "FormulationOSRC.h"
#include "FormulationDummy.h"
#include "FormulationSommerfeld.h"
#include "FormulationSteadyWave.h"
#include "FormulationUpdateEMDA.h"
#include "FormulationUpdateOO2.h"
#include "FormulationUpdateOSRC.h"

using namespace std;

static const int scal = 0;
static const int vect = 1;

static double k; // Need to be more sexy !
static double theta = 0;

Complex fSourceScal(fullVector<double>& xyz){
  double p = xyz(0) * cos(theta) + xyz(1) * sin(theta);

  return Complex(cos(k * p), sin(k * p));
  //return Complex(1, 0);
}

fullVector<Complex> fSourceVect(fullVector<double>& xyz){
  fullVector<Complex> tmp(3);

  tmp(0) = Complex(0, 0);
  tmp(1) = Complex(1, 0);
  tmp(2) = Complex(0, 0);

  return tmp;
}

void compute(const Options& option){
  // MPI //
  int nProcs;
  int myProc;
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myProc);

 // Get Type //
  int type;
  if(option.getValue("-type")[1].compare("scalar") == 0)
    type = scal;

  else if(option.getValue("-type")[1].compare("vector") == 0)
    type = vect;

  else
    throw Exception("Bad -type: %s", option.getValue("-type")[1].c_str());

  // Get Parameters //
  const string ddmType = option.getValue("-ddm")[1];
               k       = atof(option.getValue("-k")[1].c_str());
  const size_t order   = atoi(option.getValue("-o")[1].c_str());
  const size_t maxIt   = atoi(option.getValue("-max")[1].c_str());

  // DDM Formulations //
  const string emdaType("emda");
  const string oo2Type("oo2");
  const string osrcType("osrc");

  // Variables
  const double Pi = atan(1.0) * 4;
  double lc       = 0;
  double chi      = 0;
  Complex ooA     = 0;
  Complex ooB     = 0;
  int NPade       = 0;
  Complex keps    = 0;

  // EMDA Stuff
  if(ddmType == emdaType)
    chi = atof(option.getValue("-chi")[1].c_str()) * k;

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

  // OSRC Stuff
  if(ddmType == osrcType){
    double ck = atof(option.getValue("-ck")[1].c_str());
    NPade     = atoi(option.getValue("-pade")[1].c_str());
    keps      = k + Complex(0, k * ck);
  }

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume(msh);
  GroupOfElement source(msh);
  GroupOfElement infinity(msh);
  GroupOfElement ddmBorder(msh);

  // Source
  if(myProc == 0)
    source.add(msh.getFromPhysical(1000));

  // Infinity
  if(myProc == nProcs - 1)
    infinity.add(msh.getFromPhysical(2000 + nProcs - 1));

  // Volume
  volume.add(msh.getFromPhysical(100 + myProc));

  // DDM border
  if(myProc > 0)
    ddmBorder.add(msh.getFromPhysical(4000 + myProc - 1));

  if(myProc < nProcs - 1)
    ddmBorder.add(msh.getFromPhysical(4000 + myProc));

  // Full Domain //
  vector<const GroupOfElement*> domain(4);
  domain[0] = &volume;
  domain[1] = &source;
  domain[2] = &infinity;
  domain[3] = &ddmBorder;

  // Function Space //
  FunctionSpace* fs = NULL;

  if(type == scal)
    fs = new FunctionSpaceScalar(domain, order);
  else
    fs = new FunctionSpaceVector(domain, order);

  vector<const FunctionSpaceScalar*> phi(NPade); // OSRC

  for(int j = 0; j < NPade; j++)
    phi[j] = new FunctionSpaceScalar(ddmBorder, order);

  // Steady Wave Formulation //
  FormulationSteadyWave<Complex> wave(volume, *fs, k);

  // Sommerfeld
  Formulation<Complex>* sommerfeld;

  if(myProc == nProcs - 1)
    sommerfeld = new FormulationSommerfeld(infinity, *fs, k);
  else
    sommerfeld = new FormulationDummy<Complex>;

  // DDM Solution Map //
  map<Dof, Complex> ddmG;
  map<Dof, Complex> rhsG;
  FormulationHelper::initDofMap(*fs, ddmBorder, ddmG);
  FormulationHelper::initDofMap(*fs, ddmBorder, rhsG);

  // Ddm Formulation //
  DDMContext*         context = NULL;
  Formulation<Complex>*   ddm = NULL;
  Formulation<Complex>* upDdm = NULL;

  if(ddmType == emdaType){
    context = new DDMContextEMDA(ddmBorder, *fs, k, chi);
    context->setDDMDofs(ddmG);

    ddm     = new FormulationEMDA(static_cast<DDMContextEMDA&>(*context));
    upDdm   = new FormulationUpdateEMDA(static_cast<DDMContextEMDA&>(*context));
  }

  else if(ddmType == oo2Type){
    context = new DDMContextOO2(ddmBorder, *fs, ooA, ooB);
    context->setDDMDofs(ddmG);

    ddm     = new FormulationOO2(static_cast<DDMContextOO2&>(*context));
    upDdm   = new FormulationUpdateOO2(static_cast<DDMContextOO2&>(*context));
  }

  else if(ddmType == osrcType){
    context = new DDMContextOSRC(ddmBorder, *fs, phi, k, keps, NPade);
    context->setDDMDofs(ddmG);

    ddm     = new FormulationOSRC(static_cast<DDMContextOSRC&>(*context));
    upDdm   = new FormulationUpdateOSRC(static_cast<DDMContextOSRC&>(*context));
  }

  else
    throw Exception("Unknown %s DDM border term", ddmType.c_str());

  // Solve Non homogenous problem //
  System<Complex> nonHomogenous;
  nonHomogenous.addFormulation(wave);
  nonHomogenous.addFormulation(*sommerfeld);
  nonHomogenous.addFormulation(*ddm);

  // Constraint
  if(fs->isScalar())
    SystemHelper<Complex>::dirichlet(nonHomogenous, *fs, source, fSourceScal);
  else
    SystemHelper<Complex>::dirichlet(nonHomogenous, *fs, source, fSourceVect);

  // Assemble & Solve
  nonHomogenous.assemble();
  nonHomogenous.solve();

  // Solve Non homogenous DDM problem //
  context->setSystem(nonHomogenous);
  upDdm->update(); // update volume solution (at DDM border)

  System<Complex> nonHomogenousDDM;
  nonHomogenousDDM.addFormulation(*upDdm);

  nonHomogenousDDM.assemble();
  nonHomogenousDDM.solve();
  nonHomogenousDDM.getSolution(rhsG, 0);

  // DDM Solver //
  SolverDDM solver(wave, *sommerfeld, source, *context, *ddm, *upDdm, rhsG);

  solver.solve(maxIt);

  // Full Problem //
  solver.getSolution(ddmG);

  context->setDDMDofs(ddmG);
  ddm->update();

  System<Complex> full;
  full.addFormulation(wave);
  full.addFormulation(*sommerfeld);
  full.addFormulation(*ddm);

  // Constraint
  if(fs->isScalar())
    SystemHelper<Complex>::dirichlet(full, *fs, source, fSourceScal);
  else
    SystemHelper<Complex>::dirichlet(full, *fs, source, fSourceVect);

  full.assemble();
  full.solve();

  // Draw Solution //
  stringstream stream;
  stream << "circle" << myProc;

  FEMSolution<Complex> feSol;
  full.getSolution(feSol, *fs, volume);
  feSol.write(stream.str());

  // Clean //
  delete ddm;
  delete upDdm;
  delete context;
  delete sommerfeld;
  delete fs;

  for(int j = 0; j < NPade; j++)
    delete phi[j];
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-type,-max,-ddm,-chi,-lc,-ck,-pade");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
