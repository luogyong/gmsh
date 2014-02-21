/////////////////////////////////////////////////
// Templates Implementations for TermGradGrad: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include "Exception.h"
#include "ReferenceSpaceManager.h"

template<typename scalar>
TermGradGrad<scalar>::TermGradGrad(const GroupOfJacobian& goj,
                                   const Basis& basis,
                                   const Quadrature& quadrature){
  // Pre-eval //
  preEvalF(goj, quadrature, 1);

  // Init //
  init(goj, basis, quadrature);
}

template<typename scalar>
TermGradGrad<scalar>::TermGradGrad(const GroupOfJacobian& goj,
                                   const Basis& basis,
                                   const Quadrature& quadrature,
                                   scalar (*f)(const fullVector<double>& xyz)){
  // Pre-eval //
  preEvalF(goj, quadrature, f);

  // Init //
  init(goj, basis, quadrature);
}

template<typename scalar>
void TermGradGrad<scalar>::init(const GroupOfJacobian& goj,
                                const Basis& basis,
                                const Quadrature& quadrature){
  // Basis Check //
  bFunction getFunction;

  switch(basis.getForm()){
  case 0:
    getFunction = &Basis::getPreEvaluatedDerivatives;
    break;

  case 1:
    getFunction = &Basis::getPreEvaluatedFunctions;
    break;

  default:
    throw Exception
      ("A Grad Grad Term takes a 1form basis or the gradient of a 0form basis");
  }

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

  computeC(basis, getFunction, gW, cM);
  computeB(goj, gW.size(), bM);

  allocA(this->nFunction * this->nFunction);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

template<typename scalar>
void TermGradGrad<scalar>::preEvalF(const GroupOfJacobian& goj,
                                    const Quadrature& quadrature,
                                    scalar f){
  // Data //
  const size_t nPoint   = quadrature.getPoints().size1();
  const size_t nElement = goj.getAllElements().getAll().size();
  const size_t size     = nPoint * nElement;

  // Alloc //
  alpha.resize(size);

  // Populate //
  for(size_t i = 0; i < size; i++)
    alpha[i] = f;
}

template<typename scalar>
void TermGradGrad<scalar>::preEvalF(const GroupOfJacobian& goj,
                                    const Quadrature& quadrature,
                                    scalar (*f)(const fullVector<double>& xyz)){
  // Data //
  const fullMatrix<double>&                gC = quadrature.getPoints();
  const std::vector<const MElement*>& element = goj.getAllElements().getAll();
  const size_t nPoint                         = gC.size1();
  const size_t nElement                       = element.size();

  // Tmp //
  scalar             fxyz;
  double             pxyz[3];
  fullVector<double> xyz(3);

  // Alloc //
  alpha.resize(nPoint * nElement);

  // Populate //
  for(size_t e = 0; e < nElement; e++){
    for(size_t g = 0; g < nPoint; g++){
      // Compute f in the *physical* coordinate
      ReferenceSpaceManager::
        mapFromABCtoXYZ(*element[e], gC(g, 0), gC(g, 1), gC(g, 2), pxyz);

        xyz(0) = pxyz[0];
        xyz(1) = pxyz[1];
        xyz(2) = pxyz[2];

        alpha[e * nPoint + g] = f(xyz);
    }
  }
}

template<typename scalar>
TermGradGrad<scalar>::~TermGradGrad(void){
}

template<typename scalar>
void TermGradGrad<scalar>::computeC(const Basis& basis,
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
void TermGradGrad<scalar>::computeB(const GroupOfJacobian& goj,
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
        invJac = goj.getJacobian(e).getInvertJacobianMatrix();

      // Loop on Gauss Points
      k = 0;

      for(size_t g = 0; g < nG; g++){
        for(size_t a = 0; a < 3; a++){
          for(size_t b = 0; b < 3; b++){
            // Jacobians
            (*bM[s])(j, k) = 0;

            for(size_t i = 0; i < 3; i++)
              (*bM[s])(j, k) +=
                (*invJac[g]->first)(i, a) * (*invJac[g]->first)(i, b);

            (*bM[s])(j, k) *= alpha[e * nG + g] * fabs(invJac[g]->second);

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
