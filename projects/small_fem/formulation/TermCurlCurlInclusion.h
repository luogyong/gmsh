/////////////////////////////////////////////////
// Templates Implementations for TermCurlCurl: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include <cmath>
#include "Exception.h"
#include "ReferenceSpaceManager.h"

template<typename scalar>
TermCurlCurl<scalar>::
TermCurlCurl(const GroupOfJacobian& goj,
             const Basis& basis,
             const fullVector<double>& integrationWeights){

  // Basis Check //
  bFunction getFunction;

  switch(basis.getForm()){
  case 1:
    getFunction = &Basis::getPreEvaluatedDerivatives;
    break;

  case 2:
    getFunction = &Basis::getPreEvaluatedFunctions;
    break;

  default:
    throw Exception
      ("A Curl Curl Term takes a 2form basis or the curl of a 1form basis");
  }

  // Type //
  int eType = basis.getType();

  // Orientations & Functions //
  this->orientationStat = &goj.getAllElements().getOrientationStats(eType);
  this->nOrientation    = ReferenceSpaceManager::getNOrientation(eType);
  this->nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<scalar>** cM;
  fullMatrix<scalar>** bM;

  computeC(basis, getFunction, integrationWeights, cM);
  computeB(goj, integrationWeights.size(), bM);

  allocA(this->nFunction * this->nFunction);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

template<typename scalar>
TermCurlCurl<scalar>::~TermCurlCurl(void){
}

template<typename scalar>
void TermCurlCurl<scalar>::computeC(const Basis& basis,
                                    const bFunction& getFunction,
                                    const fullVector<double>& gW,
                                    fullMatrix<scalar>**& cM){
  const size_t nG = gW.size();

  // Alloc //
  cM = new fullMatrix<scalar>*[this->nOrientation];

  for(size_t s = 0; s < this->nOrientation; s++)
    cM[s] = new fullMatrix<scalar>(9 * nG, this->nFunction * this->nFunction);

  // Fill //
  //#pragma omp parallel
  for(size_t s = 0; s < this->nOrientation; s++){

    // Get functions for this Orientation
    const fullMatrix<double>& phi = (basis.*getFunction)(s);

    // fullMatrix is in *Column-major* //
    //  --> iterate on column first    //
    //       --> iterate on functions  //

    // Loop on Functions
    //#pragma omp for
    for(size_t i = 0; i < this->nFunction; i++){
      for(size_t j = 0; j < this->nFunction; j++){

        // Loop on Gauss Points
        for(size_t g = 0; g < nG; g++){
          for(size_t a = 0; a < 3; a++){
            for(size_t b = 0; b < 3; b++){
              (*cM[s])(g * 9 + a * 3 + b, i * this->nFunction + j) =
                gW(g) * phi(i, g * 3 + a) * phi(j, g * 3 + b);
            }
          }
        }
      }
    }
  }
}

template<typename scalar>
void TermCurlCurl<scalar>::computeB(const GroupOfJacobian& goj,
                                    size_t nG,
                                    fullMatrix<scalar>**& bM){
  size_t offset = 0;
  size_t j;
  size_t k;

  // Alloc //
  bM = new fullMatrix<scalar>*[this->nOrientation];

  for(size_t s = 0; s < this->nOrientation; s++)
    bM[s] = new fullMatrix<scalar>((*this->orientationStat)[s], 9 * nG);

  // Fill //

  // Despite that fullMatrix is Column-major  //
  // Row-major fill seems faster for matrix B //

  for(size_t s = 0; s < this->nOrientation; s++){
    // Loop On Element
    j = 0;

    for(size_t e = offset; e < offset + (*this->orientationStat)[s]; e++){
      // Get Jacobians
      const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
        MJac = goj.getJacobian(e).getJacobianMatrix();

      // Loop on Gauss Points
      k = 0;

      for(size_t g = 0; g < nG; g++){
        for(size_t a = 0; a < 3; a++){
          for(size_t b = 0; b < 3; b++){
            (*bM[s])(j, k) = 0;

            for(size_t i = 0; i < 3; i++)
              (*bM[s])(j, k) +=
                (*MJac[g]->first)(a, i) * (*MJac[g]->first)(b, i);

            (*bM[s])(j, k) /= fabs(MJac[g]->second);

            k++;
          }
        }
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*this->orientationStat)[s];
  }
}
