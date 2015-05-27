#include "SmallFem.h"
#include "Interpolator.h"
#include <iostream>

using namespace std;

static const Complex I  = Complex(0, 1);
static const double  Pi = M_PI;

static const Complex E0 = Complex(1, 0);
static const double  a  = 1;
static const double  b  = 1;
static const int     m  = 20;
static const int     n  = 20;

static const double  ky = m * Pi / a;
static const double  kz = n * Pi / b;
static       double  k;
static       Complex kx;

void getKx(void){
  kx = sqrt(Complex(k * k, 0) - (ky * ky) - (kz * kz));
}

vector<Complex> fVect(double x, double y, double z){
  vector<Complex> tmp(3);
  /*
  // TMm 2D
  tmp[0] = -E0 * I * ky / k * sin(ky * y) * exp(I * kx * x);
  tmp[1] = +E0 *     kx / k * cos(ky * y) * exp(I * kx * x);
  tmp[2] = Complex(0, 0);
  */

  // TEmn 3D
  tmp[0] = Complex(0, 0);
  tmp[1] = -E0 * cos(ky * y) * sin(kz * z) * exp(I * kx * x);
  tmp[2] = +E0 * sin(ky * y) * cos(kz * z) * exp(I * kx * x);

  /*
  // TMmn 3D
  tmp(0) = E0                                  * sin(ky * y) * sin(kz * z);
  tmp(1) = E0 * (-I * kx * ky) / (k*k - kx*kx) * cos(ky * y) * sin(kz * z);
  tmp(2) = E0 * (-I * kx * kz) / (k*k - kx*kx) * sin(ky * y) * cos(kz * z);
  */
  return tmp;
}

void compute(const Options& option){
  // Get Parameters //
  const size_t nDom  = atoi(option.getValue("-n")[1].c_str());
  k                  = atof(option.getValue("-k")[1].c_str());

  // Compute kx //
  getKx();

  cout << "Wavenumber: " << k     << endl
       << "# Domain:   " << nDom  << endl << flush;

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement volume(msh);

  vector<GroupOfElement*> perVolume(nDom);
  for(size_t i = 0; i < nDom; i++){
    perVolume[i] = new GroupOfElement(msh.getFromPhysical(i + 1));
    volume.add(*perVolume[i]);
  }

  // Grab solution //
  vector<map<const MVertex*, vector<Complex> > > map(nDom);
  for(size_t i = 0; i < nDom; i++){
    MapVertex vertex;
    perVolume[i]->getAllVertex(vertex);

    MapVertex::iterator it  = vertex.begin();
    MapVertex::iterator end = vertex.end();

    for(; it != end; it++){
      double x, y, z;
      x = it->first->x();
      y = it->first->y();
      z = it->first->z();

      pair<const MVertex*, vector<Complex> > tuple(it->first, fVect(x, y, z));
      map[i].insert(tuple);
    }
  }

  // Draw //
  stringstream stream;
  try{
      vector<string> name = option.getValue("-name");
      stream << name[1];
  }
  catch(...){
    stream << "waveguideAna";
  }

  for(size_t i = 0; i < nDom; i++){
    // Interpolation
    stringstream name;
    name << stream.str() << i << ".dat";

    Interpolator<Complex>::write(name.str(), map[i]);
  }

  // Clean //
  for(size_t i = 0; i < nDom; i++)
    delete perVolume[i];
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-k,-n,-name");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
