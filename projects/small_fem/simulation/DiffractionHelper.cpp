#include <cmath>

#include "DiffractionHelper.h"

using namespace std;

// Math constant //
const double Math::Pi = atan(1.0) * 4;

// Wave //
const double Wave::lambda0 = 1000000;
const double Wave::theta0  = 0.0 * Math::Pi / 180;
const double Wave::phi0    = 0.0 * Math::Pi / 180;
const double Wave::psi0    = 0.0 * Math::Pi / 180;

// Geometry //
const double Geometry::paramaille = 10;
const double Geometry::nm         = 1000;

const double Geometry::a_lat      = 300 * nm;
const double Geometry::period_x   = a_lat;
const double Geometry::period_y   = a_lat;
const double Geometry::period_z   = a_lat;

const double Geometry::PML_top    = Wave::lambda0;
const double Geometry::PML_bot    = Wave::lambda0;
const double Geometry::PML_lat    = Wave::lambda0;

const double Geometry::ro         = Wave::lambda0 / 10;

// Constant //
const double Constant::cel    = 1.0 / sqrt(Material::epsilon0 * Material::mu0);
const double Constant::Freq   = cel / Wave::lambda0;
const double Constant::omega0 = 2.0 * Math::Pi * Freq;
const double Constant::k0     = 2.0 * Math::Pi / Wave::lambda0;

const double Constant::Ae     = 1.0;
const double Constant::Ah     = Ae * sqrt(Material::epsilon0 / Material::mu0);
const double Constant::alpha0 = k0 * sin(Wave::theta0) * cos(Wave::phi0);
const double Constant::beta0  = k0 * sin(Wave::theta0) * sin(Wave::phi0);
const double Constant::gamma0 = k0 * cos(Wave::theta0);

const double Constant::Ex0    =
  Ae * cos(Wave::psi0) * cos(Wave::theta0) * cos(Wave::phi0) -
  Ae * sin(Wave::psi0) * sin(Wave::phi0);
const double Constant::Ey0    =
  Ae * cos(Wave::psi0) * cos(Wave::theta0) * sin(Wave::phi0) +
  Ae * sin(Wave::psi0) * cos(Wave::phi0);
const double Constant::Ez0    =
  -Ae * cos(Wave::psi0) * sin(Wave::theta0);
const double Constant::Hx0    =
  -1 / (omega0 * Material::mu0) * (beta0  * Ez0 - gamma0 * Ey0);
const double Constant::Hy0    =
  -1 / (omega0 * Material::mu0) * (gamma0 * Ex0 - alpha0 * Ez0);
const double Constant::Hz0    =
  -1 / (omega0 * Material::mu0) * (alpha0 * Ey0 - beta0  * Ex0);
const double Constant::Pinc   =
  0.5 * Ae * Ae * sqrt(Material::epsilon0 / Material::mu0) * cos(Wave::theta0);

// Signal //
Complex Signal::Prop(fullVector<double>& xyz){
  return
    Constant::Ae * Complex(cos(Constant::alpha0 * xyz(0) +
                               Constant::beta0  * xyz(1) +
                               Constant::gamma0 * xyz(2)) ,
                           sin(Constant::alpha0 * xyz(0) +
                               Constant::beta0  * xyz(1) +
                               Constant::gamma0 * xyz(2)));
}

fullVector<Complex> Signal::Einc(fullVector<double>& xyz){
  fullVector<Complex> ret(3);

  ret(0) = Constant::Ex0 * Prop(xyz);
  ret(1) = Constant::Ey0 * Prop(xyz);
  ret(2) = Constant::Ez0 * Prop(xyz);

  return ret;
}

fullVector<Complex> Signal::In::source(fullVector<double>& xyz){
  double kSquare = (Constant::omega0 / Constant::cel) *
                   (Constant::omega0 / Constant::cel);

  fullVector<Complex> e = Einc(xyz);
  fullMatrix<Complex> epsilon(3, 3);
  fullMatrix<Complex> epsilon1(3, 3);

  Material::In::Epsilon(xyz, epsilon);
  Material::In::Epsilon1(xyz, epsilon1);

  fullVector<Complex> ret(3);

  ret(0) = kSquare * (epsilon(0, 0) - epsilon1(0, 0)) * e(0);
  ret(1) = kSquare * (epsilon(1, 1) - epsilon1(1, 1)) * e(1);
  ret(2) = kSquare * (epsilon(2, 2) - epsilon1(2, 2)) * e(2);

  return ret;
}

fullVector<Complex> Signal::Out::source(fullVector<double>& xyz){
  double kSquare = (Constant::omega0 / Constant::cel) *
                   (Constant::omega0 / Constant::cel);

  fullVector<Complex> e = Einc(xyz);
  fullMatrix<Complex> epsilon(3, 3);
  fullMatrix<Complex> epsilon1(3, 3);

  Material::Out::Epsilon(xyz, epsilon);
  Material::Out::Epsilon1(xyz, epsilon1);

  fullVector<Complex> ret(3);

  ret(0) = kSquare * (epsilon(0, 0) - epsilon1(0, 0)) * e(0);
  ret(1) = kSquare * (epsilon(1, 1) - epsilon1(1, 1)) * e(1);
  ret(2) = kSquare * (epsilon(2, 2) - epsilon1(2, 2)) * e(2);

  return ret;
}

fullVector<Complex> Signal::PML::source(fullVector<double>& xyz){
  fullVector<Complex> ret(3);

  ret(0) = Complex(0, 0);
  ret(1) = Complex(0, 0);
  ret(2) = Complex(0, 0);

  return ret;
}

// PML //
const Complex PML::a   = Complex(1, 1);
const Complex PML::one = Complex(1, 0);

// Volume
Complex PML::Volume::sx(fullVector<double>& xyz){
  return one;
}

Complex PML::Volume::sy(fullVector<double>& xyz){
  return one;
}

Complex PML::Volume::sz(fullVector<double>& xyz){
  return one;
}

Complex PML::Volume::Lxx(fullVector<double>& xyz){
  return sy(xyz) * sz(xyz) / sx(xyz);
}

Complex PML::Volume::Lyy(fullVector<double>& xyz){
  return sz(xyz) * sx(xyz) / sy(xyz);
}

Complex PML::Volume::Lzz(fullVector<double>& xyz){
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
// Constant
const double  Material::mu0      =  4 * Math::Pi * 100.0 * Geometry::nm;
const double  Material::epsilon0 =  8.854187817E-3 * Geometry::nm;
const Complex Material::epsIn    =  Complex(9, -1);
const Complex Material::epsOut   =  Complex(1,  0);

// In
void Material::In::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsIn;
  T(1, 1) = epsIn;
  T(2, 2) = epsIn;
}

void Material::In::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut;
  T(1, 1) = epsOut;
  T(2, 2) = epsOut;
}

void Material::In::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0);
  T(1, 1) = Complex(1, 0);
  T(2, 2) = Complex(1, 0);
}

void Material::In::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0);
  T(1, 1) = Complex(1, 0);
  T(2, 2) = Complex(1, 0);
}

// Out
void Material::Out::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut;
  T(1, 1) = epsOut;
  T(2, 2) = epsOut;
}

void Material::Out::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut;
  T(1, 1) = epsOut;
  T(2, 2) = epsOut;
}

void Material::Out::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0);
  T(1, 1) = Complex(1, 0);
  T(2, 2) = Complex(1, 0);
}

void Material::Out::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0);
  T(1, 1) = Complex(1, 0);
  T(2, 2) = Complex(1, 0);
}

// XYZ
void Material::XYZ::Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut * PML::XYZ::Lxx(xyz);
  T(1, 1) = epsOut * PML::XYZ::Lyy(xyz);
  T(2, 2) = epsOut * PML::XYZ::Lzz(xyz);
}

void Material::XYZ::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut * PML::XYZ::Lxx(xyz);
  T(1, 1) = epsOut * PML::XYZ::Lyy(xyz);
  T(2, 2) = epsOut * PML::XYZ::Lzz(xyz);
}

void Material::XYZ::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) * PML::XYZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) * PML::XYZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) * PML::XYZ::Lzz(xyz);
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

  T(0, 0) = epsOut * PML::XZ::Lxx(xyz);
  T(1, 1) = epsOut * PML::XZ::Lyy(xyz);
  T(2, 2) = epsOut * PML::XZ::Lzz(xyz);
}

void Material::XZ::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut * PML::XZ::Lxx(xyz);
  T(1, 1) = epsOut * PML::XZ::Lyy(xyz);
  T(2, 2) = epsOut * PML::XZ::Lzz(xyz);
}

void Material::XZ::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) * PML::XZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) * PML::XZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) * PML::XZ::Lzz(xyz);
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

  T(0, 0) = epsOut * PML::YZ::Lxx(xyz);
  T(1, 1) = epsOut * PML::YZ::Lyy(xyz);
  T(2, 2) = epsOut * PML::YZ::Lzz(xyz);
}

void Material::YZ::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut * PML::YZ::Lxx(xyz);
  T(1, 1) = epsOut * PML::YZ::Lyy(xyz);
  T(2, 2) = epsOut * PML::YZ::Lzz(xyz);
}

void Material::YZ::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) * PML::YZ::Lxx(xyz);
  T(1, 1) = Complex(1, 0) * PML::YZ::Lyy(xyz);
  T(2, 2) = Complex(1, 0) * PML::YZ::Lzz(xyz);
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

  T(0, 0) = epsOut * PML::XY::Lxx(xyz);
  T(1, 1) = epsOut * PML::XY::Lyy(xyz);
  T(2, 2) = epsOut * PML::XY::Lzz(xyz);
}

void Material::XY::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut * PML::XY::Lxx(xyz);
  T(1, 1) = epsOut * PML::XY::Lyy(xyz);
  T(2, 2) = epsOut * PML::XY::Lzz(xyz);
}

void Material::XY::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) * PML::XY::Lxx(xyz);
  T(1, 1) = Complex(1, 0) * PML::XY::Lyy(xyz);
  T(2, 2) = Complex(1, 0) * PML::XY::Lzz(xyz);
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

  T(0, 0) = epsOut * PML::X::Lxx(xyz);
  T(1, 1) = epsOut * PML::X::Lyy(xyz);
  T(2, 2) = epsOut * PML::X::Lzz(xyz);
}

void Material::X::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut * PML::X::Lxx(xyz);
  T(1, 1) = epsOut * PML::X::Lyy(xyz);
  T(2, 2) = epsOut * PML::X::Lzz(xyz);
}

void Material::X::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) * PML::X::Lxx(xyz);
  T(1, 1) = Complex(1, 0) * PML::X::Lyy(xyz);
  T(2, 2) = Complex(1, 0) * PML::X::Lzz(xyz);
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

  T(0, 0) = epsOut * PML::Y::Lxx(xyz);
  T(1, 1) = epsOut * PML::Y::Lyy(xyz);
  T(2, 2) = epsOut * PML::Y::Lzz(xyz);
}

void Material::Y::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut * PML::Y::Lxx(xyz);
  T(1, 1) = epsOut * PML::Y::Lyy(xyz);
  T(2, 2) = epsOut * PML::Y::Lzz(xyz);
}

void Material::Y::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) * PML::Y::Lxx(xyz);
  T(1, 1) = Complex(1, 0) * PML::Y::Lyy(xyz);
  T(2, 2) = Complex(1, 0) * PML::Y::Lzz(xyz);
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

  T(0, 0) = epsOut * PML::Z::Lxx(xyz);
  T(1, 1) = epsOut * PML::Z::Lyy(xyz);
  T(2, 2) = epsOut * PML::Z::Lzz(xyz);
}

void Material::Z::Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = epsOut * PML::Z::Lxx(xyz);
  T(1, 1) = epsOut * PML::Z::Lyy(xyz);
  T(2, 2) = epsOut * PML::Z::Lzz(xyz);
}

void Material::Z::Mu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) * PML::Z::Lxx(xyz);
  T(1, 1) = Complex(1, 0) * PML::Z::Lyy(xyz);
  T(2, 2) = Complex(1, 0) * PML::Z::Lzz(xyz);
}

void Material::Z::Nu(fullVector<double>& xyz, fullMatrix<Complex>& T){
  T.scale(0);

  T(0, 0) = Complex(1, 0) / PML::Z::Lxx(xyz);
  T(1, 1) = Complex(1, 0) / PML::Z::Lyy(xyz);
  T(2, 2) = Complex(1, 0) / PML::Z::Lzz(xyz);
}
