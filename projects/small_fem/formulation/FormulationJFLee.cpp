#include <cmath>

#include "FormulationJFLee.h"
#include "FormulationJFLeeOne.h"
#include "FormulationJFLeeTwo.h"
#include "FormulationJFLeeThree.h"
#include "FormulationJFLeeFour.h"
#include "FormulationJFLeeFive.h"
#include "FormulationJFLeeSix.h"
#include "FormulationJFLeeSeven.h"
#include "FormulationJFLeeEight.h"

#include "Exception.h"

using namespace std;

FormulationJFLee::FormulationJFLee(DDMContextJFLee& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and auxiliary FunctionSpaces from DDMContext //
  const GroupOfElement&    domain = context.getDomain();
  const FunctionSpaceVector& fPhi = context.getPhi();
  const FunctionSpaceScalar& fRho = context.getRho();

  // Save field FunctionSpace: saved for update() //
  if(context.getFunctionSpace().isScalar())
    throw Exception
      ("FormulationJFLee needs a vectorial FunctionSpace for primary unknown");

  field = static_cast<const FunctionSpaceVector*>(&context.getFunctionSpace());

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationJFLee needs a uniform mesh");

  // Get Basis from primary space //
  basis = &field->getBasis(eType); // Saved for update()
  const size_t order = basis->getOrder();

  // Get auxiliary bases //
  const Basis& basisPhi = fPhi.getBasis(eType);
  const Basis& basisRho = fRho.getBasis(eType);

  // Wavenumber and JF Lee Coef //
  double   k = context.getK();
  Complex C1 = context.getC1();
  Complex C2 = context.getC2();

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Gaussian Quadrature //
  gauss = new Quadrature(eType, order, 2); // Saved for update()
  const fullMatrix<double>& gC = gauss->getPoints();

  // Pre evaluation //
  basis->preEvaluateFunctions(gC);
  basis->preEvaluateDerivatives(gC);

  basisPhi.preEvaluateFunctions(gC);
  basisPhi.preEvaluateDerivatives(gC);

  basisRho.preEvaluateFunctions(gC);
  basisRho.preEvaluateDerivatives(gC);

  // Jacobians //
  jac = new GroupOfJacobian(domain, gC, "both"); // Saved for update()

  // Local Terms //
  proj   = new TermProjectionGrad<Complex>(*jac, *basis, *gauss, *field, ddm);
  term11 = new TermGradGrad<double>       (*jac, *basis,              *gauss);
  term01 = new TermGradGrad<double>       (*jac,  basisRho, basisPhi, *gauss);
  term22 = new TermCurlCurl<double>       (*jac, *basis,              *gauss);
  term00 = new TermFieldField<double>     (*jac,  basisRho          , *gauss);
  term10 = new TermGradGrad<double>       (*jac,  basisPhi, basisRho, *gauss);

  // Formulations //
  // NB: FormulationJFLee is a friend of FormulationJFLee{One,...,Eight,} !
  //     So it can instanciate these classes...

  // Save FormulationJFLeeOne for update()
  formulationOne = new FormulationJFLeeOne(domain, *field, k, *proj);

  // Then push it in list
  fList.push_back(formulationOne);
  fList.push_back(new FormulationJFLeeTwo  (domain,  fPhi,*field,  k, *term11));
  fList.push_back(new FormulationJFLeeThree(domain,  fRho,  fPhi, C2, *term01));
  fList.push_back(new FormulationJFLeeFour (domain,  fPhi,         k, *term11));
  fList.push_back(new FormulationJFLeeFive (domain, *field, fPhi,  k, *term11));
  fList.push_back(new FormulationJFLeeSix  (domain, *field, fPhi, C1, *term22));
  fList.push_back(new FormulationJFLeeSeven(domain,  fRho           , *term00));
  fList.push_back(new FormulationJFLeeEight(domain,  fPhi,  fRho    , *term10));
}

FormulationJFLee::~FormulationJFLee(void){
  // Iterate & Delete Formulations //
  list<const FormulationBlock<Complex>*>::iterator end = fList.end();
  list<const FormulationBlock<Complex>*>::iterator it  = fList.begin();

  for(; it !=end; it++)
    delete *it;

  // Delete terms //
  delete proj;
  delete term11;
  delete term01;
  delete term22;
  delete term00;
  delete term10;

  // Delete update stuffs //
  delete jac;
  delete gauss;
}

const list<const FormulationBlock<Complex>*>&
FormulationJFLee::getFormulationBlocks(void) const{
  return fList;
}

bool FormulationJFLee::isBlock(void) const{
  return false;
}

void FormulationJFLee::update(void){
  // Delete RHS
  delete proj;

  // Get DDM Dofs from DDMContext
  const map<Dof, Complex>& ddm = context->getDDMDofs();

  // Pre-evalution
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // New RHS
  proj = new TermProjectionGrad<Complex>(*jac, *basis, *gauss, *field, ddm);

  // Update FormulationJFLeeOne (formulationOne):
  //                                          this FormulationBlock holds RHS
  formulationOne->update(*proj);
}
