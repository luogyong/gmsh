/////////////////////////////////////////////////
// Templates Implementations for TermCurlCurl: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

/////////////////////////////////////////////////
// WARNING: Jacobian matrices are transposed ! //
/////////////////////////////////////////////////

#include "Exception.h"
#include "ReferenceSpaceManager.h"

template<typename scalar>
TermCurlCurl<scalar>::TermCurlCurl(const GroupOfJacobian& goj,
                                   const Basis& basis,
                                   const Quadrature& quadrature){
  // Dummy Pre evaluation //
  preEvalDummy(goj, quadrature);

  // Init //
  init(goj, basis, quadrature);
}

template<typename scalar>
TermCurlCurl<scalar>::TermCurlCurl(const GroupOfJacobian& goj,
                                   const Basis& basis,
                                   const Quadrature& quadrature,
                                   void  (*f)(fullVector<double>& xyz,
                                              fullMatrix<scalar>&   T)){
  // Tensor Pre evaluation //
  preEvalT(goj, quadrature, f);

  // Init //
  init(goj, basis, quadrature);
}

template<typename scalar>
void TermCurlCurl<scalar>::init(const GroupOfJacobian& goj,
                                const Basis& basis,
                                const Quadrature& quadrature){
  // Basis Check //
  BFunction getFunction;

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
  TJac.clear();
}

template<typename scalar>
TermCurlCurl<scalar>::~TermCurlCurl(void){
}

template<typename scalar>
void TermCurlCurl<scalar>::computeC(const Basis& basis,
                                    const BFunction& getFunction,
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
        // Get TJac
        const fullMatrix<scalar>& myTJac = this->TJac[e * nG + g];

        for(size_t a = 0; a < 3; a++){
          for(size_t b = 0; b < 3; b++){
            (*bM[s])(j, k) = 0;

            for(size_t i = 0; i < 3; i++)
              // Waring Jacobian matrices are transposed in Gmsh
              (*bM[s])(j, k) += myTJac(a, i) * (*MJac[g]->first)(b, i);
                //(*MJac[g]->first)(a, i) * (*MJac[g]->first)(b, i);

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

template<typename scalar>
void TermCurlCurl<scalar>::preEvalDummy(const GroupOfJacobian& goj,
                                        const Quadrature& quadrature){
  // Data //
  const size_t nPoint   = quadrature.getPoints().size1();
  const size_t nElement = goj.getAllElements().getAll().size();
  const size_t     size = nElement * nPoint;

  // Alloc //
  TJac.resize(size);

  #pragma omp parallel for // First touch
  for(size_t i = 0; i < size; i++)
    TJac[i].resize(3, 3);

  // Populate //
  #pragma omp parallel for
  for(size_t e = 0; e < nElement; e++){
    for(size_t g = 0; g < nPoint; g++){
      // Get Jacobians for element 'e' at point 'g'
      const fullMatrix<double>& J =
        *goj.getJacobian(e).getJacobianMatrix()[g]->first;

      // Get Reference to TJac[e][g]
      fullMatrix<scalar>& myTJac = TJac[e * nPoint + g];

      // Copy jacobian matrix for element 'e' at point 'g' in TJac
      // Hand done because of template programming (column first)
      myTJac(0, 0) = J(0, 0); // Column 0
      myTJac(1, 0) = J(1, 0); //  |
      myTJac(2, 0) = J(2, 0); //  _

      myTJac(0, 1) = J(0, 1); // Column 1
      myTJac(1, 1) = J(1, 1); //  |
      myTJac(2, 1) = J(2, 1); //  _

      myTJac(0, 2) = J(0, 2); // Column 2
      myTJac(1, 2) = J(1, 2); //  |
      myTJac(2, 2) = J(2, 2); //  _
    }
  }
}

template<typename scalar>
void TermCurlCurl<scalar>::preEvalT(const GroupOfJacobian& goj,
                                    const Quadrature& quadrature,
                                    void  (*f)(fullVector<double>& xyz,
                                               fullMatrix<scalar>& T)){
  // Data //
  const fullMatrix<double>&                gC = quadrature.getPoints();
  const std::vector<const MElement*>  element = goj.getAllElements().getAll();
  const size_t                       nPoint   = gC.size1();
  const size_t                       nElement = element.size();
  const size_t                           size = nElement * nPoint;

  // Alloc //
  TJac.resize(size);

  #pragma omp parallel for // First touch
  for(size_t i = 0; i < size; i++)
    TJac[i].resize(3, 3, true); // Alloc 3x3 null matrices

  // Populate //
  #pragma omp parallel
  {
    // Temps (one for each thread) //
    fullVector<double>  xyz(3);
    double             pxyz[3];

    // Tensor (one for each thread) //
    fullMatrix<scalar> T(3, 3);

    #pragma omp for
    for(size_t e = 0; e < nElement; e++){
      for(size_t g = 0; g < nPoint; g++){
        // Compute Tensor in the *physical* coordinate
        ReferenceSpaceManager::
          mapFromABCtoXYZ(*element[e], gC(g, 0), gC(g, 1), gC(g, 2), pxyz);

        xyz(0) = pxyz[0];
        xyz(1) = pxyz[1];
        xyz(2) = pxyz[2];

        f(xyz, T);

        // Get Jacobians for element 'e' at point 'g'
        const fullMatrix<double>& J =
          *goj.getJacobian(e).getJacobianMatrix()[g]->first;

        // Get Reference to TJac[e][g]
        fullMatrix<scalar>& myTJac = TJac[e * nPoint + g];

        // A = T * J -- (Taking that Jacobian matrix is transposed in Gmsh)
        // Hand done since template and 3x3 matrices (BLAS not needed)
        for(int i = 0; i < 3; i++)
          for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
              myTJac(i, j) += J(i, k) * T(j, k);
      }
    }
  }
}
