#include "TermProjectionField.h"
#include "SmallFem.h"

// Real Version //
// ------------ //
template<>
void TermProjectionField<double>::computeC(const Basis& basis,
                                           const fullVector<double>& gW,
                                           fullMatrix<double>**& cM){
  const size_t nG = gW.size();

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(size_t g = 0; g < nG; g++)
      for(size_t i = 0; i < nFunction; i++)
        (*cM[s])(g, i) = gW(g) * phi(i, g);
  }
}

// Complex Version //
// --------------- //
template<>
void TermProjectionField<Complex>::computeC(const Basis& basis,
                                            const fullVector<double>& gW,
                                            fullMatrix<Complex>**& cM){
  const size_t nG = gW.size();

  // Alloc //
  cM = new fullMatrix<Complex>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<Complex>(nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(size_t g = 0; g < nG; g++)
      for(size_t i = 0; i < nFunction; i++)
        (*cM[s])(g, i) = Complex(gW(g) * phi(i, g), 0);
  }
}
