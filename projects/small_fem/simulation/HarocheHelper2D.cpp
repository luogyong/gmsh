#include "HarocheHelper2D.h"

using namespace std;

// PML //
const Complex PML::a   = Complex(1, -1);
const Complex PML::one = Complex(1,  0);
const Complex PML::sz  = Complex(1,  0);

double PML::Xmax  = 0.0308668948122371;
double PML::Ymax  = 0.0210668948122371;
double PML::SizeX = 0.0058668948122371;
double PML::SizeY = 0.0058668948122371;

//double PML::SizeX = 0.0586689481223709;
//double PML::SizeY = 0.0586689481223709;

double PML::kHaroche = 1070.95584773;

Complex PML::dampingX(fullVector<double>& xyz){
  double        f = Xmax + SizeX - fabs(xyz(0));
  double    overF = 1 / (f) - 1 / (SizeX);

  return Complex(1, -overF / kHaroche);
  //return a;
}

Complex PML::dampingY(fullVector<double>& xyz){
  double        f = Ymax + SizeY - fabs(xyz(1));
  double    overF = 1 / (f) - 1 / (SizeY);

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

Complex PML::Air::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz / sx(xyz);
}

Complex PML::Air::Lyy(fullVector<double>& xyz){
  return sz * sx(xyz) / sy(xyz);
}

Complex PML::Air::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz;
}

// XY
Complex PML::XY::sx(fullVector<double>& xyz){
  return dampingX(xyz);
}

Complex PML::XY::sy(fullVector<double>& xyz){
  return dampingY(xyz);
}

Complex PML::XY::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz / sx(xyz);
}

Complex PML::XY::Lyy(fullVector<double>& xyz){
  return sz * sx(xyz) / sy(xyz);
}

Complex PML::XY::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz;
}

// X
Complex PML::X::sx(fullVector<double>& xyz){
  return dampingX(xyz);
}

Complex PML::X::sy(fullVector<double>& xyz){
  return one;
}

Complex PML::X::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz / sx(xyz);
}

Complex PML::X::Lyy(fullVector<double>& xyz){
  return sz * sx(xyz) / sy(xyz);
}

Complex PML::X::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz;
}

// Y
Complex PML::Y::sx(fullVector<double>& xyz){
  return one;
}

Complex PML::Y::sy(fullVector<double>& xyz){
  return dampingY(xyz);
}

Complex PML::Y::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz / sx(xyz);
}

Complex PML::Y::Lyy(fullVector<double>& xyz){
  return sz * sx(xyz) / sy(xyz);
}

Complex PML::Y::Lzz(fullVector<double>& xyz){
  return sx(xyz) * sy(xyz) / sz;
}


// Material //

// Speed of light
const double Material::cSquare = 2.997924580105029e+08 * 2.997924580105029e+08;

// Epsilon and Nu
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
