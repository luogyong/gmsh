#include <complex>
#include "Interpolator.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

using namespace std;

// Real Implementation //
// ------------------- //

template<>
void Interpolator<double>::
interpolate(const MElement& element,
            const FunctionSpace& fs,
            const std::vector<double>& coef,
            const fullVector<double>& xyz,
            fullVector<double>& value){

  // Scalar or Vector ?
  const bool isScalar = fs.isScalar();

  if(isScalar){
    // Function Space Scalar
    const FunctionSpaceScalar* fsScalar =
      static_cast<const FunctionSpaceScalar*>(&fs);

    // Alloc value
    value.resize(1);

    // Interpolate
    value(0) = fsScalar->interpolate(element, coef, xyz);
  }

  else{
    // Function Space Vector
    const FunctionSpaceVector* fsVector =
      static_cast<const FunctionSpaceVector*>(&fs);

    // Alloc & Interpolate
    value = fsVector->interpolate(element, coef, xyz);
  }
}

// Complex Implementation //
// ---------------------- //

template<>
void Interpolator<complex<double> >::
interpolate(const MElement& element,
            const FunctionSpace& fs,
            const std::vector<complex<double> >& coef,
            const fullVector<double>& xyz,
            fullVector<complex<double> >& value){

  // Real & Imaginary //
  const size_t size = coef.size();

  vector<double> coefReal(size);
  vector<double> coefImag(size);

  for(size_t i = 0; i < size; i++)
    coefReal[i] = coef[i].real();

  for(size_t i = 0; i < size; i++)
    coefImag[i] = coef[i].imag();

  // Scalar or Vector ?
  const bool isScalar = fs.isScalar();

  if(isScalar){
    // Function Space Scalar
    const FunctionSpaceScalar* fsScalar =
      static_cast<const FunctionSpaceScalar*>(&fs);

    // Alloc value
    value.resize(1);

    // Interpolate
    double real = fsScalar->interpolate(element, coefReal, xyz);
    double imag = fsScalar->interpolate(element, coefImag, xyz);

    // Complex
    value(0) = complex<double>(real, imag);
  }

  else{
    // Function Space Vector
    const FunctionSpaceVector* fsVector =
      static_cast<const FunctionSpaceVector*>(&fs);

    // Alloc value
    value.resize(3);

    // Interpolate
    fullVector<double> real = fsVector->interpolate(element, coefReal, xyz);
    fullVector<double> imag = fsVector->interpolate(element, coefImag, xyz);

    // Complex
    value(0) = complex<double>(real(0), imag(0));
    value(1) = complex<double>(real(1), imag(1));
    value(2) = complex<double>(real(2), imag(2));
  }
}
