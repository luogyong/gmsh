#include <complex>
#include "FEMSolution.h"

using namespace std;

// Real Implementation //
// ------------------- //

template<>
void FEMSolution<double>::
toLagrange(const MElement& element,
           const vector<BasisLagrange*>& lagrange,
           const vector<double>& fsCoef,
           const FunctionSpace& fs,
           vector<double>& lCoef){

  // Element Type //
  const int eType = element.getType();

  // Projection //
  if(fs.isScalar())
    lCoef =
      lagrange[eType]->project(element, fsCoef,
                               static_cast<const FunctionSpaceScalar&>(fs));
  else
    lCoef =
      lagrange[eType]->project(element, fsCoef,
                               static_cast<const FunctionSpaceVector&>(fs));
}

template<>
void FEMSolution<double>::
toPView(GModel& model, map<int, vector<double> >& data,
        size_t step, double time, int partition, int nComp){

  pView->addData(&model, data, step, time, partition, nComp);
}

// Complex Implementation //
// ---------------------- //

template<>
void FEMSolution<complex<double> >::
toLagrange(const MElement& element,
           const vector<BasisLagrange*>& lagrange,
           const vector<complex<double> >& fsCoef,
           const FunctionSpace& fs,
           vector<complex<double> >& lCoef){

  // Element Type //
  const int eType = element.getType();

  // Real & Imaginary //
  const size_t fsSize = fsCoef.size();

  vector<double> fsCoefReal(fsSize);
  vector<double> fsCoefImag(fsSize);

  vector<double> lCoefReal;
  vector<double> lCoefImag;

  for(size_t i = 0; i < fsSize; i++)
    fsCoefReal[i] = fsCoef[i].real();

  for(size_t i = 0; i < fsSize; i++)
    fsCoefImag[i] = fsCoef[i].imag();

  // Projection //
  if(fs.isScalar()){
    lCoefReal =
      lagrange[eType]->project(element, fsCoefReal,
                               static_cast<const FunctionSpaceScalar&>(fs));
    lCoefImag =
      lagrange[eType]->project(element, fsCoefImag,
                               static_cast<const FunctionSpaceScalar&>(fs));
  }

  else{
    lCoefReal =
      lagrange[eType]->project(element, fsCoefReal,
                               static_cast<const FunctionSpaceVector&>(fs));
    lCoefImag =
      lagrange[eType]->project(element, fsCoefImag,
                               static_cast<const FunctionSpaceVector&>(fs));
  }

  // Complex Number //
  const size_t lSize = lCoefReal.size();
  lCoef.resize(lSize);

  for(size_t i = 0; i < lSize; i++)
    lCoef[i] = complex<double>(lCoefReal[i], lCoefImag[i]);
}

template<>
void FEMSolution<complex<double> >::
toPView(GModel& model, map<int, vector<complex<double> > >& data,
        size_t step, double time, int partition, int nComp){

  // Split Data //
  // New Real / Imag
  map<int, vector<double> > real;
  map<int, vector<double> > imag;

  // Iterate on Data
  const map<int, vector<complex<double> > >::iterator end = data.end();
  map<int, vector<complex<double> > >::iterator        it = data.begin();

  for(; it != end; it++){
    // Element
    int eNum = it->first;

    // Number of coefficients
    size_t nCoef = it->second.size();

    // New vectors
    vector<double> vReal(nCoef);
    vector<double> vImag(nCoef);

    for(size_t i = 0; i < nCoef; i++)
      vReal[i] = it->second[i].real();

    for(size_t i = 0; i < nCoef; i++)
      vImag[i] = it->second[i].imag();

    real.insert(pair<int, vector<double> >(eNum, vReal));
    imag.insert(pair<int, vector<double> >(eNum, vImag));
  }

  // Add to PView //
  pView->addData(&model, real, 2 * step + 0, time, partition, nComp);
  pView->addData(&model, imag, 2 * step + 1, time, partition, nComp);
}
