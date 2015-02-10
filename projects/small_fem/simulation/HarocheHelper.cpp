#include "SmallFem.h"
#include "Exception.h"
#include "HarocheHelper.h"
#include <fstream>

using namespace std;

// PML //
const Complex PML::a   = Complex(1, -1);
const Complex PML::one = Complex(1,  0);

double PML::Xmax;
double PML::Ymax;
double PML::Zmax;
double PML::SizeX;
double PML::SizeY;
double PML::SizeZ;
double PML::kHaroche;

void PML::read(string filename){
  ifstream stream(filename.c_str(), ifstream::in);
  if(!stream.is_open())
    throw Exception("HarocheHelper: cannot open PML file (%s)",
                    filename.c_str());
  stream >> SizeX
         >> SizeY
         >> SizeZ
         >> Xmax
         >> Ymax
         >> Zmax
         >> kHaroche;
  stream.close();
}

Complex PML::dampingX(fullVector<double>& xyz){
  double     f = Xmax + SizeX - fabs(xyz(0));
  double overF = 1 / (f) - 1 / (SizeX);

  return Complex(1, -overF / kHaroche);
  //return a;
}

Complex PML::dampingY(fullVector<double>& xyz){
  double     f = Ymax + SizeY - fabs(xyz(1));
  double overF = 1 / (f) - 1 / (SizeY);

  return Complex(1, -overF / kHaroche);
  //return a;
}

Complex PML::dampingZ(fullVector<double>& xyz){
  double     f = Zmax + SizeZ - fabs(xyz(2));
  double overF = 1 / (f) - 1 / (SizeZ);

  return Complex(1, -overF / kHaroche);
  //return a;
}

// Air
Complex PML::Air::sx(fullVector<double>& xyz){
  return one;
}

Complex PML::Air::sy(fullVector<double>& xyz){
  return one;
}

Complex PML::Air::sz(fullVector<double>& xyz){
  return one;
}

Complex PML::Air::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::Air::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::Air::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// XYZ
Complex PML::XYZ::sx(fullVector<double>& xyz){
  return dampingX(xyz);
}

Complex PML::XYZ::sy(fullVector<double>& xyz){
  return dampingY(xyz);
}

Complex PML::XYZ::sz(fullVector<double>& xyz){
  return dampingZ(xyz);
}

Complex PML::XYZ::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::XYZ::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::XYZ::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// XZ
Complex PML::XZ::sx(fullVector<double>& xyz){
  return dampingX(xyz);
}

Complex PML::XZ::sy(fullVector<double>& xyz){
  return one;
}

Complex PML::XZ::sz(fullVector<double>& xyz){
  return dampingZ(xyz);
}

Complex PML::XZ::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::XZ::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::XZ::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// YZ
Complex PML::YZ::sx(fullVector<double>& xyz){
  return one;
}

Complex PML::YZ::sy(fullVector<double>& xyz){
  return dampingY(xyz);
}

Complex PML::YZ::sz(fullVector<double>& xyz){
  return dampingZ(xyz);
}

Complex PML::YZ::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::YZ::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::YZ::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// XY
Complex PML::XY::sx(fullVector<double>& xyz){
  return dampingX(xyz);
}

Complex PML::XY::sy(fullVector<double>& xyz){
  return dampingY(xyz);
}

Complex PML::XY::sz(fullVector<double>& xyz){
  return one;
}

Complex PML::XY::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::XY::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::XY::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// X
Complex PML::X::sx(fullVector<double>& xyz){
  return dampingX(xyz);
}

Complex PML::X::sy(fullVector<double>& xyz){
  return one;
}

Complex PML::X::sz(fullVector<double>& xyz){
  return one;
}

Complex PML::X::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::X::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::X::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// Y
Complex PML::Y::sx(fullVector<double>& xyz){
  return one;
}

Complex PML::Y::sy(fullVector<double>& xyz){
  return dampingY(xyz);
}

Complex PML::Y::sz(fullVector<double>& xyz){
  return one;
}

Complex PML::Y::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::Y::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::Y::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// Z
Complex PML::Z::sx(fullVector<double>& xyz){
  return one;
}

Complex PML::Z::sy(fullVector<double>& xyz){
  return one;
}

Complex PML::Z::sz(fullVector<double>& xyz){
  return dampingZ(xyz);
}

Complex PML::Z::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::Z::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::Z::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz(xyz);
}

// Material //

// Speed of light
const double Material::cSquare = 2.997924580105029e+08 * 2.997924580105029e+08;

// Epsilon / (c^2) and Nu
void Material::Air::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / cSquare;
  T(1, 1) = Complex(1, 0) / cSquare;
  T(2, 2) = Complex(1, 0) / cSquare;
}

void Material::Air::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0);
  T(1, 1) = Complex(1, 0);
  T(2, 2) = Complex(1, 0);
}

void Material::Air::MuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / cSquare;
  T(1, 1) = Complex(1, 0) / cSquare;
  T(2, 2) = Complex(1, 0) / cSquare;
}

void Material::Air::OverMuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = cSquare;
  T(1, 1) = cSquare;
  T(2, 2) = cSquare;
}

// XYZ
void Material::XYZ::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XYZ::Lxx(xyz) / cSquare;
  T(1, 1) = PML::XYZ::Lyy(xyz) / cSquare;
  T(2, 2) = PML::XYZ::Lzz(xyz) / cSquare;
}

void Material::XYZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XYZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XYZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XYZ::Lzz(xyz);
}

void Material::XYZ::MuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XYZ::Lxx(xyz) * PML::XYZ::Lxx(xyz) / cSquare;
  T(1, 1) = PML::XYZ::Lyy(xyz) * PML::XYZ::Lyy(xyz) / cSquare;
  T(2, 2) = PML::XYZ::Lzz(xyz) * PML::XYZ::Lzz(xyz) / cSquare;
}

void Material::XYZ::OverMuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = cSquare / (PML::XYZ::Lxx(xyz) * PML::XYZ::Lxx(xyz));
  T(1, 1) = cSquare / (PML::XYZ::Lyy(xyz) * PML::XYZ::Lyy(xyz));
  T(2, 2) = cSquare / (PML::XYZ::Lzz(xyz) * PML::XYZ::Lzz(xyz));
}

// XZ
void Material::XZ::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XZ::Lxx(xyz) / cSquare;
  T(1, 1) = PML::XZ::Lyy(xyz) / cSquare;
  T(2, 2) = PML::XZ::Lzz(xyz) / cSquare;
}

void Material::XZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XZ::Lzz(xyz);
}

void Material::XZ::MuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XZ::Lxx(xyz) * PML::XZ::Lxx(xyz) / cSquare;
  T(1, 1) = PML::XZ::Lyy(xyz) * PML::XZ::Lyy(xyz) / cSquare;
  T(2, 2) = PML::XZ::Lzz(xyz) * PML::XZ::Lzz(xyz) / cSquare;
}

void Material::XZ::OverMuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = cSquare / (PML::XZ::Lxx(xyz) * PML::XZ::Lxx(xyz));
  T(1, 1) = cSquare / (PML::XZ::Lyy(xyz) * PML::XZ::Lyy(xyz));
  T(2, 2) = cSquare / (PML::XZ::Lzz(xyz) * PML::XZ::Lzz(xyz));
}

// YZ
void Material::YZ::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::YZ::Lxx(xyz) / cSquare;
  T(1, 1) = PML::YZ::Lyy(xyz) / cSquare;
  T(2, 2) = PML::YZ::Lzz(xyz) / cSquare;
}

void Material::YZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::YZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::YZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::YZ::Lzz(xyz);
}

void Material::YZ::MuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::YZ::Lxx(xyz) * PML::YZ::Lxx(xyz) / cSquare;
  T(1, 1) = PML::YZ::Lyy(xyz) * PML::YZ::Lyy(xyz) / cSquare;
  T(2, 2) = PML::YZ::Lzz(xyz) * PML::YZ::Lzz(xyz) / cSquare;
}

void Material::YZ::OverMuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = cSquare / (PML::YZ::Lxx(xyz) * PML::YZ::Lxx(xyz));
  T(1, 1) = cSquare / (PML::YZ::Lyy(xyz) * PML::YZ::Lyy(xyz));
  T(2, 2) = cSquare / (PML::YZ::Lzz(xyz) * PML::YZ::Lzz(xyz));
}

// XY
void Material::XY::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XY::Lxx(xyz) / cSquare;
  T(1, 1) = PML::XY::Lyy(xyz) / cSquare;
  T(2, 2) = PML::XY::Lzz(xyz) / cSquare;
}

void Material::XY::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XY::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XY::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XY::Lzz(xyz);
}

void Material::XY::MuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XY::Lxx(xyz) * PML::XY::Lxx(xyz) / cSquare;
  T(1, 1) = PML::XY::Lyy(xyz) * PML::XY::Lyy(xyz) / cSquare;
  T(2, 2) = PML::XY::Lzz(xyz) * PML::XY::Lzz(xyz) / cSquare;
}

void Material::XY::OverMuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = cSquare / (PML::XY::Lxx(xyz) * PML::XY::Lxx(xyz));
  T(1, 1) = cSquare / (PML::XY::Lyy(xyz) * PML::XY::Lyy(xyz));
  T(2, 2) = cSquare / (PML::XY::Lzz(xyz) * PML::XY::Lzz(xyz));
}

// X
void Material::X::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::X::Lxx(xyz) / cSquare;
  T(1, 1) = PML::X::Lyy(xyz) / cSquare;
  T(2, 2) = PML::X::Lzz(xyz) / cSquare;
}

void Material::X::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::X::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::X::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::X::Lzz(xyz);
}

void Material::X::MuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::X::Lxx(xyz) * PML::X::Lxx(xyz) / cSquare;
  T(1, 1) = PML::X::Lyy(xyz) * PML::X::Lyy(xyz) / cSquare;
  T(2, 2) = PML::X::Lzz(xyz) * PML::X::Lzz(xyz) / cSquare;
}

void Material::X::OverMuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = cSquare / (PML::X::Lxx(xyz) * PML::X::Lxx(xyz));
  T(1, 1) = cSquare / (PML::X::Lyy(xyz) * PML::X::Lyy(xyz));
  T(2, 2) = cSquare / (PML::X::Lzz(xyz) * PML::X::Lzz(xyz));
}

// Y
void Material::Y::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::Y::Lxx(xyz) / cSquare;
  T(1, 1) = PML::Y::Lyy(xyz) / cSquare;
  T(2, 2) = PML::Y::Lzz(xyz) / cSquare;
}

void Material::Y::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::Y::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::Y::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::Y::Lzz(xyz);
}

void Material::Y::MuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::Y::Lxx(xyz) * PML::Y::Lxx(xyz) / cSquare;
  T(1, 1) = PML::Y::Lyy(xyz) * PML::Y::Lyy(xyz) / cSquare;
  T(2, 2) = PML::Y::Lzz(xyz) * PML::Y::Lzz(xyz) / cSquare;
}

void Material::Y::OverMuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = cSquare / (PML::Y::Lxx(xyz) * PML::Y::Lxx(xyz));
  T(1, 1) = cSquare / (PML::Y::Lyy(xyz) * PML::Y::Lyy(xyz));
  T(2, 2) = cSquare / (PML::Y::Lzz(xyz) * PML::Y::Lzz(xyz));
}

// Z
void Material::Z::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::Z::Lxx(xyz) / cSquare;
  T(1, 1) = PML::Z::Lyy(xyz) / cSquare;
  T(2, 2) = PML::Z::Lzz(xyz) / cSquare;
}

void Material::Z::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::Z::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::Z::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::Z::Lzz(xyz);
}

void Material::Z::MuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::Z::Lxx(xyz) * PML::Z::Lxx(xyz) / cSquare;
  T(1, 1) = PML::Z::Lyy(xyz) * PML::Z::Lyy(xyz) / cSquare;
  T(2, 2) = PML::Z::Lzz(xyz) * PML::Z::Lzz(xyz) / cSquare;
}

void Material::Z::OverMuEps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = cSquare / (PML::Z::Lxx(xyz) * PML::Z::Lxx(xyz));
  T(1, 1) = cSquare / (PML::Z::Lyy(xyz) * PML::Z::Lyy(xyz));
  T(2, 2) = cSquare / (PML::Z::Lzz(xyz) * PML::Z::Lzz(xyz));
}
