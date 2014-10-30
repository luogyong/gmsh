#include <cmath>
#include <iostream>

#include "SmallFem.h"

#include "SolverDDM.h"

#include "DDMContextEMDA.h"
#include "DDMContextOO2.h"
#include "DDMContextJFLee.h"
#include "DDMContextOSRCScalar.h"
#include "DDMContextOSRCVector.h"

#include "System.h"
#include "SystemHelper.h"

#include "FormulationHelper.h"

#include "FormulationOO2.h"
#include "FormulationEMDA.h"
#include "FormulationJFLee.h"
#include "FormulationOSRCScalar.h"
#include "FormulationOSRCVector.h"

#include "FormulationSommerfeld.h"
#include "FormulationSteadyWave.h"

#include "FormulationUpdateEMDA.h"
#include "FormulationUpdateOO2.h"
#include "FormulationUpdateJFLee.h"
#include "FormulationUpdateOSRCScalar.h"
#include "FormulationUpdateOSRCVector.h"

using namespace std;

static const int scal = 0;
static const int vect = 1;

Complex fSourceScal(fullVector<double>& xyz){
  return Complex(1, 0);
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
  const double k       = atof(option.getValue("-k")[1].c_str());
  const size_t order   = atoi(option.getValue("-o")[1].c_str());
  const size_t maxIt   = atoi(option.getValue("-max")[1].c_str());

  // DDM Formulations //
  const string emdaType("emda");
  const string oo2Type("oo2");
  const string osrcType("osrc");
  const string jflType("jfl");

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
  else if(ddmType == oo2Type){
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
  else if(ddmType == osrcType){
    double ck = atof(option.getValue("-ck")[1].c_str());
    NPade     = atoi(option.getValue("-pade")[1].c_str());
    keps      = k + Complex(0, k * ck);
  }

  // Jin Fa Lee Stuff
  else if(ddmType == jflType){
    lc = atof(option.getValue("-lc")[1].c_str());
  }

  // Unknown Stuff
  else
    throw Exception("DDM: Formulation %s is not known", ddmType.c_str());

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume(msh);
  GroupOfElement source(msh);
  GroupOfElement infinity(msh);
  GroupOfElement ddmBorder(msh);

  if(myProc == 0){
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
  vector<const GroupOfElement*> domain(4);
  domain[0] = &volume;
  domain[1] = &source;
  domain[2] = &infinity;
  domain[3] = &ddmBorder;

  // Function Space //
  FunctionSpace* fs =  NULL;

  if(type == scal)
    fs = new FunctionSpaceScalar(domain, order);
  else
    fs = new FunctionSpaceVector(domain, order);

  // OSRC
  vector<const FunctionSpaceScalar*> OSRCScalPhi;
  vector<const FunctionSpaceVector*> OSRCVectPhi;
  vector<const FunctionSpaceScalar*> OSRCVectRho;
  FunctionSpaceVector*               OSRCVectR = NULL;

  if(ddmType == osrcType && type == scal){
    OSRCScalPhi.resize(NPade);

    for(int j = 0; j < NPade; j++)
      OSRCScalPhi[j] = new FunctionSpaceScalar(ddmBorder, order);
  }

  if(ddmType == osrcType && type == vect){
    OSRCVectPhi.resize(NPade);
    OSRCVectRho.resize(NPade);

    for(int j = 0; j < NPade; j++)
      OSRCVectPhi[j] = new FunctionSpaceVector(ddmBorder, order);

    if(order == 0)
      for(int j = 0; j < NPade; j++)
        OSRCVectRho[j] = new FunctionSpaceScalar(ddmBorder, 1);
    else
      for(int j = 0; j < NPade; j++)
        OSRCVectRho[j] = new FunctionSpaceScalar(ddmBorder, order);

    OSRCVectR = new FunctionSpaceVector(ddmBorder, order);
  }

  // Jin Fa Lee
  FunctionSpaceVector* JFPhi = NULL;
  FunctionSpaceScalar* JFRho = NULL;

  if(ddmType == jflType){
    JFPhi = new FunctionSpaceVector(ddmBorder, order);

    if(order == 0)
      JFRho = new FunctionSpaceScalar(ddmBorder, 1);
    else
      JFRho = new FunctionSpaceScalar(ddmBorder, order);
  }

  // Steady Wave Formulation //
  FormulationSteadyWave<Complex> wave(volume, *fs, k);
  FormulationSommerfeld          sommerfeld(infinity, *fs, k);

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

    ddm   = new FormulationEMDA(static_cast<DDMContextEMDA&>(*context));
    upDdm = new FormulationUpdateEMDA(static_cast<DDMContextEMDA&>(*context));
  }

  else if(ddmType == oo2Type){
    context = new DDMContextOO2(ddmBorder, *fs, ooA, ooB);
    context->setDDMDofs(ddmG);

    ddm   = new FormulationOO2(static_cast<DDMContextOO2&>(*context));
    upDdm = new FormulationUpdateOO2(static_cast<DDMContextOO2&>(*context));
  }

  else if(ddmType == osrcType && type == scal){
    context = new DDMContextOSRCScalar
                                  (ddmBorder, *fs, OSRCScalPhi, k, keps,
                                   NPade, M_PI / 4.);
    context->setDDMDofs(ddmG);

    ddm   = new FormulationOSRCScalar
                                 (static_cast<DDMContextOSRCScalar&>(*context));
    upDdm = new FormulationUpdateOSRCScalar
                                 (static_cast<DDMContextOSRCScalar&>(*context));
  }

  else if(ddmType == osrcType && type == vect){
    context = new DDMContextOSRCVector
                                  (ddmBorder,
                                   *fs, OSRCVectPhi, OSRCVectRho, *OSRCVectR,
                                   k, keps, NPade, M_PI / 2.);
    context->setDDMDofs(ddmG);

    ddm   = new FormulationOSRCVector
                                 (static_cast<DDMContextOSRCVector&>(*context));
    upDdm = new FormulationUpdateOSRCVector
                                 (static_cast<DDMContextOSRCVector&>(*context));
  }

  else if(ddmType == jflType){
    context = new DDMContextJFLee(ddmBorder, *fs, *JFPhi, *JFRho, k, lc);
    context->setDDMDofs(ddmG);

    ddm   = new FormulationJFLee(static_cast<DDMContextJFLee&>(*context));
    upDdm = new FormulationUpdateJFLee(static_cast<DDMContextJFLee&>(*context));
  }

  else
    throw Exception("Unknown %s DDM border term", ddmType.c_str());

  // Solve Non homogenous problem //
  System<Complex> nonHomogenous;
  nonHomogenous.addFormulation(wave);
  nonHomogenous.addFormulation(sommerfeld);
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
  SolverDDM solver(wave, sommerfeld, source, *context, *ddm, *upDdm, rhsG);
  solver.solve(maxIt);

  // Full Problem //
  solver.getSolution(ddmG);

  context->setDDMDofs(ddmG);
  ddm->update();

  System<Complex> full;
  full.addFormulation(wave);
  full.addFormulation(sommerfeld);
  full.addFormulation(*ddm);

  // Constraint
  if(fs->isScalar())
    SystemHelper<Complex>::dirichlet(full, *fs, source, fSourceScal);
  else
    SystemHelper<Complex>::dirichlet(full, *fs, source, fSourceVect);

  full.assemble();
  full.solve();

  full.getSolution(ddmG, 0);

  // Draw Solution //
  try{
    option.getValue("-nopos");
  }
  catch(...){
    stringstream stream;
    stream << "ddm" << myProc;

    FEMSolution<Complex> feSol;
    full.getSolution(feSol, *fs, volume);
    feSol.write(stream.str());
  }

  // Clean //
  delete ddm;
  delete upDdm;
  delete context;
  delete fs;

  if(JFPhi)
    delete JFPhi;

  if(JFRho)
    delete JFRho;

  if((int)(OSRCScalPhi.size()) == NPade)
    for(int j = 0; j < NPade; j++)
      delete OSRCScalPhi[j];

  if((int)(OSRCVectPhi.size()) == NPade)
    for(int j = 0; j < NPade; j++)
      delete OSRCVectPhi[j];

  if((int)(OSRCVectRho.size()) == NPade)
    for(int j = 0; j < NPade; j++)
      delete OSRCVectRho[j];

  if(OSRCVectR)
    delete OSRCVectR;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-type,-max,-ddm,-chi,-lc,-ck,-pade,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
