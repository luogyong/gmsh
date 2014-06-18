#include "HarocheHelper.h"

using namespace std;

// PML //
const Complex PML::a   = Complex(1, -1);
const Complex PML::one = Complex(1,  0);

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
  return a;
}

Complex PML::XYZ::sy(fullVector<double>& xyz){
  return a;
}

Complex PML::XYZ::sz(fullVector<double>& xyz){
  return a;
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
  return a;
}

Complex PML::XZ::sy(fullVector<double>& xyz){
  return one;
}

Complex PML::XZ::sz(fullVector<double>& xyz){
  return a;
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
  return a;
}

Complex PML::YZ::sz(fullVector<double>& xyz){
  return a;
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
  return a;
}

Complex PML::XY::sy(fullVector<double>& xyz){
  return a;
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
  return a;
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
  return a;
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
  return a;
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
const double Material::oneOverCSquare = 1 / (2.997924580105029e+08 *
                                             2.997924580105029e+08);

// Epsilon / (c^2) and Nu
void Material::Air::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) * oneOverCSquare;
  T(1, 1) = Complex(1, 0) * oneOverCSquare;
  T(2, 2) = Complex(1, 0) * oneOverCSquare;
}

void Material::Air::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0);
  T(1, 1) = Complex(1, 0);
  T(2, 2) = Complex(1, 0);
}

// XYZ
void Material::XYZ::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XYZ::Lxx(xyz) * oneOverCSquare;
  T(1, 1) = PML::XYZ::Lyy(xyz) * oneOverCSquare;
  T(2, 2) = PML::XYZ::Lzz(xyz) * oneOverCSquare;
}

void Material::XYZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XYZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XYZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XYZ::Lzz(xyz);
}

// XZ
void Material::XZ::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XZ::Lxx(xyz) * oneOverCSquare;
  T(1, 1) = PML::XZ::Lyy(xyz) * oneOverCSquare;
  T(2, 2) = PML::XZ::Lzz(xyz) * oneOverCSquare;
}

void Material::XZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XZ::Lzz(xyz);
}

// YZ
void Material::YZ::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::YZ::Lxx(xyz) * oneOverCSquare;
  T(1, 1) = PML::YZ::Lyy(xyz) * oneOverCSquare;
  T(2, 2) = PML::YZ::Lzz(xyz) * oneOverCSquare;
}

void Material::YZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::YZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::YZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::YZ::Lzz(xyz);
}

// XY
void Material::XY::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XY::Lxx(xyz) * oneOverCSquare;
  T(1, 1) = PML::XY::Lyy(xyz) * oneOverCSquare;
  T(2, 2) = PML::XY::Lzz(xyz) * oneOverCSquare;
}

void Material::XY::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XY::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XY::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XY::Lzz(xyz);
}

// X
void Material::X::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::X::Lxx(xyz) * oneOverCSquare;
  T(1, 1) = PML::X::Lyy(xyz) * oneOverCSquare;
  T(2, 2) = PML::X::Lzz(xyz) * oneOverCSquare;
}

void Material::X::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::X::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::X::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::X::Lzz(xyz);
}

// Y
void Material::Y::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::Y::Lxx(xyz) * oneOverCSquare;
  T(1, 1) = PML::Y::Lyy(xyz) * oneOverCSquare;
  T(2, 2) = PML::Y::Lzz(xyz) * oneOverCSquare;
}

void Material::Y::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::Y::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::Y::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::Y::Lzz(xyz);
}

// Z
void Material::Z::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::Z::Lxx(xyz) * oneOverCSquare;
  T(1, 1) = PML::Z::Lyy(xyz) * oneOverCSquare;
  T(2, 2) = PML::Z::Lzz(xyz) * oneOverCSquare;
}

void Material::Z::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::Z::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::Z::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::Z::Lzz(xyz);
}
