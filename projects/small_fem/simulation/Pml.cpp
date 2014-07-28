#include <iostream>
#include <complex>
#include <cmath>

#include "Mesh.h"
#include "System.h"
#include "SystemHelper.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "FormulationSteadyWave.h"
#include "FormulationPML.h"

#include "Interpolator.h"
#include "NodeSolution.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

class PMLData{
private:
  static double k;
  static double SigmaMax;
  static double SigmaXmax;
  static double SigmaYmax;
  static double SizePML;
  static double SizePMLX;
  static double SizePMLY;
  static double Xmax;
  static double Ymax;

public:
   PMLData(void);
  ~PMLData(void);

  static void        setK(double k);
  static Complex    fMass(fullVector<double>& xyz);
  static void  fStiffness(fullVector<double>& xyz, fullMatrix<Complex>& tensor);
  static void  fOverStiffness(fullVector<double>& xyz,
                              fullMatrix<Complex>& tensor);

private:
  static double  DampingProfileX(fullVector<double>& xyz);
  static double  DampingProfileY(fullVector<double>& xyz);
  static double  SigmaX(fullVector<double>& xyz);
  static double  SigmaY(fullVector<double>& xyz);
  static Complex Kx(fullVector<double>& xyz);
  static Complex Ky(fullVector<double>& xyz);
};

double PMLData::k         = 0;
double PMLData::SigmaMax  = 1;
double PMLData::SigmaXmax = SigmaMax;
double PMLData::SigmaYmax = SigmaMax;
double PMLData::SizePML   = 8;
double PMLData::SizePMLX  = SizePML;
double PMLData::SizePMLY  = SizePML;
double PMLData::Xmax      = 10;
double PMLData::Ymax      = 10;

fullVector<Complex> fVolVec(fullVector<double>& xyz){
  fullVector<Complex> ret(3);

  ret(0) = Complex(0, 0);
  ret(1) = Complex(0, 0);
  ret(2) = Complex(0, 0);

  return ret;
}

fullVector<Complex> fSourceVec(fullVector<double>& xyz){
  fullVector<Complex> ret(3);

  ret(0) = Complex(0, 0);
  ret(1) = Complex(1, 0);
  ret(2) = Complex(0, 0);

  return ret;
}

fullVector<Complex> fInfinityVec(fullVector<double>& xyz){
  fullVector<Complex> ret(3);

  ret(0) = Complex(0, 0);
  ret(1) = Complex(0, 0);
  ret(2) = Complex(0, 0);

  return ret;
}

Complex fSourceScal(fullVector<double>& xyz){
  return Complex(1, 0);
}

Complex fInfinityScal(fullVector<double>& xyz){
  return Complex(0, 0);
}

void compute(const Options& option){
  // Start Timer //
  Timer timer, assemble, solve;
  timer.start();

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume     = msh.getFromPhysical(7);
  GroupOfElement outerSpace = msh.getFromPhysical(8);
  GroupOfElement source     = msh.getFromPhysical(5);
  GroupOfElement infinity   = msh.getFromPhysical(4);

  // Full Domain //
  vector<const GroupOfElement*> domain(3);
  domain[0] = &volume;
  domain[1] = &outerSpace;
  domain[2] = &source;

  // Get Parameters //
  const double k     = atoi(option.getValue("-k")[1].c_str());
  const size_t order = atoi(option.getValue("-o")[1].c_str());

  // PML Data //
  PMLData::setK(k);
  Complex (*fMass)(fullVector<double>& xyz)     = PMLData::fMass;
  void    (*fStif)(fullVector<double>& xyz,
                   fullMatrix<Complex>& tensor) = PMLData::fStiffness;
  void    (*fOverStif)(fullVector<double>& xyz,
                       fullMatrix<Complex>& tensor) = PMLData::fOverStiffness;

  // Formulation //
  assemble.start();
  FunctionSpaceVector fs(domain, order);
  /*
  FormulationSteadyWave<Complex> wave(volume, fs, k);
  FormulationPML pml(outerSpace, fs, k, fStif, fMass);
  */
  FormulationSteadyWave<Complex> wave(volume,    fs, k);
  FormulationSteadyWave<Complex> pml(outerSpace, fs, k, fOverStif, fStif, fVolVec);

  // System //
  System<Complex> sys;
  sys.addFormulation(wave);
  sys.addFormulation(pml);

  SystemHelper<Complex>::dirichlet(sys, fs, source,   fSourceVec);
  //SystemHelper<Complex>::dirichlet(sys, fs, infinity, fInfinityVec);

  cout << "Free Space (Order: "  << order
       << " --- Wavenumber: "    << k
       << "): " << sys.getSize() << endl;

  // Assemble //
  sys.assemble();
  assemble.stop();

  cout << "Assembled: " << assemble.time() << assemble.unit()
       << endl << flush;

  // Solve //
  solve.start();
  sys.solve();
  solve.stop();

  cout << "Solved: " << solve.time() << solve.unit()
       << endl << flush;

  // Write Sol //
  try{
    option.getValue("-nopos");
  }
  catch(...){

    // Get Dofs In //
    set<Dof> dof;
    fs.getKeys(volume, dof);

    set<Dof>::iterator    end = dof.end();
    set<Dof>::iterator     it = dof.begin();
    map<Dof, Complex>   solIn;

    for(; it != end; it++)
      solIn.insert(pair<Dof, Complex>(*it, 0));

    sys.getSolution(solIn, 0);

    // Get Dofs Out //
    dof.clear();
    fs.getKeys(outerSpace, dof);
    map<Dof, Complex> solOut;

    it  = dof.begin();
    end = dof.end();

    for(; it != end; it++)
      solOut.insert(pair<Dof, Complex>(*it, 0));

    sys.getSolution(solOut, 0);

    // Get Visu Mesh //
    Mesh           visuMsh(option.getValue("-ref")[1]);
    GroupOfElement visuGoeIn  = visuMsh.getFromPhysical(7);
    GroupOfElement visuGoeOut = visuMsh.getFromPhysical(8);

    // Vertex, Value Map //
    std::map<const MVertex*, vector<Complex> > mapOut;
    Interpolator<Complex>::interpolate(outerSpace, visuGoeOut, fs, solOut, mapOut);

    std::map<const MVertex*, vector<Complex> > mapIn;
    Interpolator<Complex>::interpolate(volume, visuGoeIn, fs, solIn, mapIn);

    // Print //
    NodeSolution<Complex> nSolOut;
    nSolOut.addNodeValue(0, 0, visuMsh, mapOut);
    nSolOut.write("pmlRefOut");

    NodeSolution<Complex> nSolIn;
    nSolIn.addNodeValue(0, 0, visuMsh, mapIn);
    nSolIn.write("pmlRefIn");

    FEMSolution<Complex> feSolIn;
    FEMSolution<Complex> feSolOut;
    feSolIn.setModulusPhase();
    feSolOut.setModulusPhase();

    sys.getSolution(feSolIn,  fs, volume);
    sys.getSolution(feSolOut, fs, outerSpace);

    feSolIn.write("pmlIn");
    feSolOut.write("pmlOut");
  }

  // Timer -- Finalize -- Return //
  timer.stop();

  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-ref,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}

PMLData::PMLData(void){
}

PMLData::~PMLData(void){
}

void PMLData::setK(double k){
  PMLData::k = k;
}

double PMLData::DampingProfileX(fullVector<double>& xyz){
  return SigmaXmax / SizePMLX * (fabs(xyz(0)) - Xmax);
}

double PMLData::DampingProfileY(fullVector<double>& xyz){
  return SigmaYmax / SizePMLY * (fabs(xyz(1)) - Ymax);
}

double PMLData::SigmaX(fullVector<double>& xyz){
  return 0.5 * (DampingProfileX(xyz) + fabs(DampingProfileX(xyz)));
}

double PMLData::SigmaY(fullVector<double>& xyz){
  return 0.5 * (DampingProfileY(xyz) + fabs(DampingProfileY(xyz)));
}

Complex PMLData::Kx(fullVector<double>& xyz){
  return Complex(1, SigmaX(xyz) / k);
}

Complex PMLData::Ky(fullVector<double>& xyz){
  return Complex(1, SigmaY(xyz) / k);
}

Complex PMLData::fMass(fullVector<double>& xyz){
  return Kx(xyz) * Ky(xyz);
}

void PMLData::fStiffness(fullVector<double>& xyz, fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = Ky(xyz) / Kx(xyz);
  tensor(1, 1) = Kx(xyz) / Ky(xyz);
  tensor(2, 2) = Complex(0, 0);
}

void PMLData::fOverStiffness(fullVector<double>& xyz,
                             fullMatrix<Complex>& tensor){
  tensor.scale(0);

  tensor(0, 0) = Kx(xyz) / Ky(xyz);
  tensor(1, 1) = Ky(xyz) / Kx(xyz);
  tensor(2, 2) = Complex(0, 0);
}
