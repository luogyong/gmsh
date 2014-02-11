#include "TermProjectionGrad.h"
#include "SmallFem.h"

// Real Version //
// ------------ //
template<>
void TermProjectionGrad<double>::computeC(const Basis& basis,
                                          const bFunction& getFunction,
                                          const fullVector<double>& gW,
                                          fullMatrix<double>**& cM){
  const size_t nG = gW.size();
  size_t k;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(3 * nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      (basis.*getFunction)(s);

    // Loop on Gauss Points
    k = 0;

    for(size_t g = 0; g < nG; g++){
      for(size_t a = 0; a < 3; a++){
        for(size_t i = 0; i < nFunction; i++)
          (*cM[s])(k, i) = gW(g) * phi(i, k);

        k++;
      }
    }
  }
}

// Complex Version //
// --------------- //
template<>
void TermProjectionGrad<Complex>::computeC(const Basis& basis,
                                           const bFunction& getFunction,
                                           const fullVector<double>& gW,
                                           fullMatrix<Complex>**& cM){
  const size_t nG = gW.size();
  size_t k;

  // Alloc //
  cM = new fullMatrix<Complex>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<Complex>(3 * nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      (basis.*getFunction)(s);

    // Loop on Gauss Points
    k = 0;

    for(size_t g = 0; g < nG; g++){
      for(size_t a = 0; a < 3; a++){
        for(size_t i = 0; i < nFunction; i++)
          (*cM[s])(k, i) = Complex(gW(g) * phi(i, k), 0);

        k++;
      }
    }
  }
}
