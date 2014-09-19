#include "HarocheHelper2D.h"

using namespace std;

// PML //
const Complex PML::a   = Complex(1, -1);
const Complex PML::one = Complex(1,  0);
const Complex PML::sz  = Complex(1,  0);

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
  return a;
}

Complex PML::XY::sy(fullVector<double>& xyz){
  return a;
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
  return a;
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
  return a;
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
const double Material::oneOverCSquare = 1 / (2.997924580105029e+08 *
                                             2.997924580105029e+08);

// Epsilon and Nu
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
