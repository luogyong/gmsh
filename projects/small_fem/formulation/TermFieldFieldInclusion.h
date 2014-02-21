///////////////////////////////////////////////////
// Templates Implementations for TermFieldField: //
// Inclusion compilation model                   //
//                                               //
// Damn you gcc: we want 'export' !              //
///////////////////////////////////////////////////

#include "Exception.h"
#include "ReferenceSpaceManager.h"

template<typename scalar>
TermFieldField<scalar>::TermFieldField(const GroupOfJacobian& goj,
                                       const Basis& basis,
                                       const Quadrature& quadrature){
  // Basis Check //
  if(basis.getForm() != 0)
    throw Exception("A Field Field Term must use a 0form basis");

  // Type //
  int eType = basis.getType();

  // Orientations & Functions //
  this->orientationStat = &goj.getAllElements().getOrientationStats(eType);
  this->nOrientation    = ReferenceSpaceManager::getNOrientation(eType);
  this->nFunction       = basis.getNFunction();

  // Get Integration Data
  //const fullMatrix<double>& gC = quadrature.getPoints();
  const fullVector<double>& gW = quadrature.getWeights();

  // Compute //
  fullMatrix<scalar>** cM;
  fullMatrix<scalar>** bM;

  computeC(basis, gW, cM);
  computeB(goj, gW.size(), bM);

  allocA(this->nFunction * this->nFunction);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

template<typename scalar>
TermFieldField<scalar>::~TermFieldField(void){
}

template<typename scalar>
void TermFieldField<scalar>::computeC(const Basis& basis,
                                      const fullVector<double>& gW,
                                      fullMatrix<scalar>**& cM){
  const size_t nG = gW.size();
  size_t l;

  // Alloc //
  cM = new fullMatrix<scalar>*[this->nOrientation];

  for(size_t s = 0; s < this->nOrientation; s++)
    cM[s] = new fullMatrix<scalar>(nG, this->nFunction * this->nFunction);

  // Fill //
  for(size_t s = 0; s < this->nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi = basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(size_t g = 0; g < nG; g++){
      // Loop on Functions
      l = 0;

      for(size_t i = 0; i < this->nFunction; i++){
        for(size_t j = 0; j < this->nFunction; j++){
          (*cM[s])(g, l) = gW(g) * phi(i, g) * phi(j, g);
          l++;
        }
      }
    }
  }
}

template<typename scalar>
void TermFieldField<scalar>::computeB(const GroupOfJacobian& goj,
                                      size_t nG,
                                      fullMatrix<scalar>**& bM){
  size_t offset = 0;
  size_t j;

  // Alloc //
  bM = new fullMatrix<scalar>*[this->nOrientation];

  for(size_t s = 0; s < this->nOrientation; s++)
    bM[s] = new fullMatrix<scalar>((*this->orientationStat)[s], nG);

  // Fill //
  for(size_t s = 0; s < this->nOrientation; s++){
    // Loop On Element
    j = 0;

    for(size_t e = offset; e < offset + (*this->orientationStat)[s]; e++){
      // Get Jacobians
      const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
        jacM = goj.getJacobian(e).getJacobianMatrix();

      // Loop on Gauss Points
      for(size_t g = 0; g < nG; g++)
        (*bM[s])(j, g) = fabs(jacM[g]->second);;

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*this->orientationStat)[s];
  }
}
