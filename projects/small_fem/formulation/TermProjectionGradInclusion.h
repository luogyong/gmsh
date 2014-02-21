///////////////////////////////////////////////////////
// Templates Implementations for TermProjectionGrad: //
// Inclusion compilation model                       //
//                                                   //
// Damn you gcc: we want 'export' !                  //
///////////////////////////////////////////////////////

#include "Exception.h"
#include "ReferenceSpaceManager.h"

template<typename scalar>
TermProjectionGrad<scalar>::
TermProjectionGrad(const GroupOfJacobian& goj,
                   const Basis& basis,
                   const Quadrature& quadrature,
                   fullVector<scalar> (*f)(fullVector<double>& xyz)){

  // Save F and Evaluator //
  Eval evaluator = &TermProjectionGrad<scalar>::fContainer;
  this->f        = f;

  // Init //
  init(goj, basis, quadrature, evaluator);
}

template<typename scalar>
TermProjectionGrad<scalar>::
TermProjectionGrad(const GroupOfJacobian& goj,
                   const Basis& basis,
                   const Quadrature& quadrature,
                   const FunctionSpaceScalar& fs,
                   const std::map<Dof, scalar>& dof){

  // Save FunctionSpace, Dof values and Evaluator //
  Eval evaluator = &TermProjectionGrad<scalar>::interpolateGrad;
  this->fsScalar = &fs;
  this->dofValue = &dof;

  // Init //
  init(goj, basis, quadrature, evaluator);
}

template<typename scalar>
void TermProjectionGrad<scalar>::init(const GroupOfJacobian& goj,
                                      const Basis& basis,
                                      const Quadrature& quadrature,
                                      const Eval& evaluator){
  // Basis Check //
  BFunction getFunction;

  switch(basis.getForm()){
  case 0:
    getFunction = &Basis::getPreEvaluatedDerivatives;
    break;

  case 1:
    getFunction = &Basis::getPreEvaluatedFunctions;
    break;

  default:
    throw
      Exception
      ("A Grad Term must use a 1form basis or the gradient of a 0form basis");
  }

  // Type //
  int eType = basis.getType();

  // Orientations & Functions //
  this->orientationStat = &goj.getAllElements().getOrientationStats(eType);
  this->nOrientation    = ReferenceSpaceManager::getNOrientation(eType);
  this->nFunction       = basis.getNFunction();

  // Get Integration Data
  const fullMatrix<double>& gC = quadrature.getPoints();
  const fullVector<double>& gW = quadrature.getWeights();

  // Compute //
  fullMatrix<scalar>** cM;
  fullMatrix<scalar>** bM;

  computeC(basis, getFunction, gW, cM);
  computeB(goj, gC, evaluator, bM);

  this->allocA(this->nFunction);
  this->computeA(bM, cM);

  // Clean up //
  this->clean(bM, cM);
}

template<typename scalar>
TermProjectionGrad<scalar>::~TermProjectionGrad(void){
}

template<typename scalar>
void TermProjectionGrad<scalar>::computeC(const Basis& basis,
                                          const BFunction& getFunction,
                                          const fullVector<double>& gW,
                                          fullMatrix<scalar>**& cM){
  const size_t nG = gW.size();
  size_t k;

  // Alloc //
  cM = new fullMatrix<scalar>*[this->nOrientation];

  for(size_t s = 0; s < this->nOrientation; s++)
    cM[s] = new fullMatrix<scalar>(3 * nG, this->nFunction);

  // Fill //
  for(size_t s = 0; s < this->nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi = (basis.*getFunction)(s);

    // Loop on Gauss Points
    k = 0;

    for(size_t g = 0; g < nG; g++){
      for(size_t a = 0; a < 3; a++){
        for(size_t i = 0; i < this->nFunction; i++)
          (*cM[s])(k, i) = gW(g) * phi(i, k);

        k++;
      }
    }
  }
}

template<typename scalar>
void TermProjectionGrad<scalar>::computeB(const GroupOfJacobian& goj,
                                          const fullMatrix<double>& gC,
                                          const Eval& evaluator,
                                          fullMatrix<scalar>**& bM){
  const size_t nG = gC.size1();
  size_t offset = 0;
  size_t j;

  fullVector<double> xyz(3);
  double             pxyz[3];
  fullVector<scalar> fxyz;

  // Alloc //
  bM = new fullMatrix<scalar>*[this->nOrientation];

  for(size_t s = 0; s < this->nOrientation; s++)
    bM[s] = new fullMatrix<scalar>((*this->orientationStat)[s], 3 * nG);

  // Fill //
  for(size_t s = 0; s < this->nOrientation; s++){
    // Loop On Element
    j = 0;

    for(size_t e = offset; e < offset + (*this->orientationStat)[s]; e++){
      // Get Jacobians
      const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
        invJac = goj.getJacobian(e).getInvertJacobianMatrix();

      // Get Element
      const MElement& element = goj.getAllElements().get(e);

      // Loop on Gauss Points
      for(size_t g = 0; g < nG; g++){
        // Compute f in the *physical* coordinate
        ReferenceSpaceManager::
          mapFromABCtoXYZ(element, gC(g, 0), gC(g, 1), gC(g, 2), pxyz);

        xyz(0) = pxyz[0];
        xyz(1) = pxyz[1];
        xyz(2) = pxyz[2];

        fxyz = (this->*evaluator)(element, xyz);

        // Compute B
        (*bM[s])(j, g * 3)     = 0;
        (*bM[s])(j, g * 3 + 1) = 0;
        (*bM[s])(j, g * 3 + 2) = 0;

        for(size_t i = 0; i < 3; i++){
          (*bM[s])(j, g * 3)     += fxyz(i) * (*invJac[g]->first)(i, 0);
          (*bM[s])(j, g * 3 + 1) += fxyz(i) * (*invJac[g]->first)(i, 1);
          (*bM[s])(j, g * 3 + 2) += fxyz(i) * (*invJac[g]->first)(i, 2);
        }

        (*bM[s])(j, g * 3)     = (*bM[s])(j, g * 3)    *fabs(invJac[g]->second);
        (*bM[s])(j, g * 3 + 1) = (*bM[s])(j, g * 3 + 1)*fabs(invJac[g]->second);
        (*bM[s])(j, g * 3 + 2) = (*bM[s])(j, g * 3 + 2)*fabs(invJac[g]->second);
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*this->orientationStat)[s];
  }
}
