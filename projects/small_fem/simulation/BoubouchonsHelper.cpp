#include "SmallFem.h"
#include "BoubouchonsHelper.h"

using namespace std;

// PML //
const Complex PML::a   = Complex(1, -1);
const Complex PML::one = Complex(1,  0);

Complex PML::dampingX(fullVector<double>& xyz){
  return a;
}

Complex PML::dampingY(fullVector<double>& xyz){
  return a;
}

Complex PML::dampingZ(fullVector<double>& xyz){
  return a;
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

// Material (Epsilon and Nu) //

// XYZ
void Material::XYZ::Eps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XYZ::Lxx(xyz);
  T(1, 1) = PML::XYZ::Lyy(xyz);
  T(2, 2) = PML::XYZ::Lzz(xyz);
}

void Material::XYZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XYZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XYZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XYZ::Lzz(xyz);
}

// XZ
void Material::XZ::Eps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XZ::Lxx(xyz);
  T(1, 1) = PML::XZ::Lyy(xyz);
  T(2, 2) = PML::XZ::Lzz(xyz);
}

void Material::XZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XZ::Lzz(xyz);
}

// YZ
void Material::YZ::Eps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::YZ::Lxx(xyz);
  T(1, 1) = PML::YZ::Lyy(xyz);
  T(2, 2) = PML::YZ::Lzz(xyz);
}

void Material::YZ::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::YZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::YZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::YZ::Lzz(xyz);
}

// XY
void Material::XY::Eps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::XY::Lxx(xyz);
  T(1, 1) = PML::XY::Lyy(xyz);
  T(2, 2) = PML::XY::Lzz(xyz);
}

void Material::XY::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::XY::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::XY::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::XY::Lzz(xyz);
}

// X
void Material::X::Eps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::X::Lxx(xyz);
  T(1, 1) = PML::X::Lyy(xyz);
  T(2, 2) = PML::X::Lzz(xyz);
}

void Material::X::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::X::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::X::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::X::Lzz(xyz);
}

// Y
void Material::Y::Eps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::Y::Lxx(xyz);
  T(1, 1) = PML::Y::Lyy(xyz);
  T(2, 2) = PML::Y::Lzz(xyz);
}

void Material::Y::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::Y::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::Y::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::Y::Lzz(xyz);
}

// Z
void Material::Z::Eps(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = PML::Z::Lxx(xyz);
  T(1, 1) = PML::Z::Lyy(xyz);
  T(2, 2) = PML::Z::Lzz(xyz);
}

void Material::Z::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::Z::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::Z::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::Z::Lzz(xyz);
}
