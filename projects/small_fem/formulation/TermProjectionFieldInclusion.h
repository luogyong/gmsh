////////////////////////////////////////////////////////
// Templates Implementations for TermProjectionField: //
// Inclusion compilation model                        //
//                                                    //
// Damn you gcc: we want 'export' !                   //
////////////////////////////////////////////////////////

#include "Exception.h"
#include "ReferenceSpaceManager.h"

template<typename scalar>
TermProjectionField<scalar>::
TermProjectionField(const GroupOfJacobian& goj,
                    const Basis& basis,
                    const fullVector<double>& integrationWeights,
                    const fullMatrix<double>& integrationPoints,
                    scalar (*f)(fullVector<double>& xyz)){

  // Save F and Evaluator //
  Eval evaluator = &TermProjectionField<scalar>::fContainer;
  this->f        = f;

  // Init //
  init(goj, basis, integrationWeights, integrationPoints, evaluator);
}

template<typename scalar>
TermProjectionField<scalar>::
TermProjectionField(const GroupOfJacobian& goj,
                    const Basis& basis,
                    const fullVector<double>& integrationWeights,
                    const fullMatrix<double>& integrationPoints,
                    const FunctionSpaceScalar& fs,
                    const std::map<Dof, scalar>& dof){

  // Save FunctionSpace, Dof values and Evaluator //
  Eval evaluator = &TermProjectionField<scalar>::interpolate;
  this->fsScalar = &fs;
  this->dofValue = &dof;

  // Init //
  init(goj, basis, integrationWeights, integrationPoints, evaluator);
}

template<typename scalar>
void TermProjectionField<scalar>::
init(const GroupOfJacobian& goj,
     const Basis& basis,
     const fullVector<double>& integrationWeights,
     const fullMatrix<double>& integrationPoints,
     const Eval& evaluator){

  // Basis Check //
  if(basis.getForm() != 0)
    throw Exception("A Field Term must use a 0form basis");

  // Type //
  int eType = basis.getType();

  // Orientations & Function //
  this->orientationStat = &goj.getAllElements().getOrientationStats(eType);
  this->nOrientation    = ReferenceSpaceManager::getNOrientation(eType);
  this->nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<scalar>** cM;
  fullMatrix<scalar>** bM;

  computeC(basis, integrationWeights, cM);
  computeB(goj, basis, integrationPoints, evaluator, bM);

  this->allocA(this->nFunction);
  this->computeA(bM, cM);

  // Clean up //
  this->clean(bM, cM);
}

template<typename scalar>
TermProjectionField<scalar>::~TermProjectionField(void){
}

template<typename scalar>
void TermProjectionField<scalar>::computeC(const Basis& basis,
                                           const fullVector<double>& gW,
                                           fullMatrix<scalar>**& cM){
  const size_t nG = gW.size();

  // Alloc //
  cM = new fullMatrix<scalar>*[this->nOrientation];

  for(size_t s = 0; s < this->nOrientation; s++)
    cM[s] = new fullMatrix<scalar>(nG, this->nFunction);

  // Fill //
  for(size_t s = 0; s < this->nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi = basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(size_t g = 0; g < nG; g++)
      for(size_t i = 0; i < this->nFunction; i++)
        (*cM[s])(g, i) = gW(g) * phi(i, g);
  }
}

template<typename scalar>
void TermProjectionField<scalar>::computeB(const GroupOfJacobian& goj,
                                           const Basis& basis,
                                           const fullMatrix<double>& gC,
                                           const Eval& evaluator,
                                           fullMatrix<scalar>**& bM){
  const size_t nG = gC.size1();
  size_t offset = 0;
  size_t j;

  fullVector<double> xyz(3);
  double             pxyz[3];
  scalar             fxyz;

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

      // Get Element
      const MElement& element = goj.getAllElements().get(e);

      for(size_t g = 0; g < nG; g++){
        // Compute f in the *physical* coordinate
        ReferenceSpaceManager::
          mapFromABCtoXYZ(element, gC(g, 0), gC(g, 1), gC(g, 2), pxyz);

        xyz(0) = pxyz[0];
        xyz(1) = pxyz[1];
        xyz(2) = pxyz[2];

        fxyz = (this->*evaluator)(element, xyz);

        // Compute B
        (*bM[s])(j, g) = fxyz * fabs(jacM[g]->second);
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*this->orientationStat)[s];
  }
}
