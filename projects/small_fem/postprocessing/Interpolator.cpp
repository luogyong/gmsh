#include "SmallFem.h"
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
void Interpolator<Complex>::
interpolate(const MElement& element,
            const FunctionSpace& fs,
            const std::vector<Complex>& coef,
            const fullVector<double>& xyz,
            fullVector<Complex>& value){

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
    value(0) = Complex(real, imag);
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
    value(0) = Complex(real(0), imag(0));
    value(1) = Complex(real(1), imag(1));
    value(2) = Complex(real(2), imag(2));
  }
}

template<>
void Interpolator<double>::write(string filename,
                                 const map<const MVertex*,
                                           vector<double> >& value){
  // Open //
  ofstream stream;
  stream.open(filename.c_str());

  // Header //
  stream << value.size()                 << endl
         << value.begin()->second.size() << endl
         << "real"                       << endl;

  // Iterators //
  map<const MVertex*, vector<double> >::const_iterator  it = value.begin();
  map<const MVertex*, vector<double> >::const_iterator end = value.end();

  // Write //
  map<int, vector<double> > data;

  // Populate
  for(; it != end; it++)
    data.insert(pair<int, vector<double> >(it->first->getNum(), it->second));

  // Dump
  dump(stream, data);

  // Close //
  stream.close();
}

template<>
void Interpolator<Complex>::write(string filename,
                                  const map<const MVertex*,
                                            vector<Complex> >& value){
  // Open //
  ofstream stream;
  stream.open(filename.c_str());

  // Header //
  stream << "complex"                    << endl
         << value.size()                 << endl
         << value.begin()->second.size() << endl;

  // Iterators //
  map<const MVertex*, vector<Complex> >::const_iterator  it = value.begin();
  map<const MVertex*, vector<Complex> >::const_iterator end = value.end();

  // Write //
  map<int, vector<double> > data;

  // Populate
  const int      dim = it->second.size();
  vector<double> tmp(dim * 2);

  for(; it != end; it++){
    for(int i = 0; i < dim; i++)
      tmp[i] = it->second[i].real();

    for(int i = 0, j = dim; i < dim; i++, j++)
      tmp[j] = it->second[i].imag();

    data.insert(pair<int, vector<double> >(it->first->getNum(), tmp));
  }

  // Dump
  dump(stream, data);

  // Close //
  stream.close();
}
