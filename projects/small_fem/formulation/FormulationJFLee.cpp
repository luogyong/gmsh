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

  field  = static_cast<const FunctionSpaceVector*>
                                                 (&context.getFunctionSpace());
  fieldG = static_cast<const FunctionSpaceVector*>
                                                 (&context.getFunctionSpaceG());

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationJFLee needs a uniform mesh");

  // Get Basis from primary space //
  basisE             = &field->getBasis(eType); // Saved for update()
  const size_t order = basisE->getOrder();

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
  basisE->preEvaluateFunctions(gC);
  basisE->preEvaluateDerivatives(gC);

  basisPhi.preEvaluateFunctions(gC);
  basisPhi.preEvaluateDerivatives(gC);

  basisRho.preEvaluateFunctions(gC);
  basisRho.preEvaluateDerivatives(gC);

  // Jacobians //
  jac = new GroupOfJacobian(domain, gC, "both"); // Saved for update()

  // Local Terms //
  proj   = new TermProjectionGrad<Complex>(*jac, *basisE, *gauss, *fieldG, ddm);
  termPE = new TermGradGrad<double>       (*jac,  basisPhi, *basisE  ,  *gauss);
  termRP = new TermGradGrad<double>       (*jac,  basisRho,  basisPhi,  *gauss);
  termPP = new TermGradGrad<double>       (*jac,  basisPhi,  basisPhi,  *gauss);
  termEP = new TermGradGrad<double>       (*jac, *basisE  ,  basisPhi,  *gauss);
  term22 = new TermCurlCurl<double>       (*jac, *basisE  ,  basisPhi,  *gauss);
  termRR = new TermFieldField<double>     (*jac,  basisRho           ,  *gauss);
  termPR = new TermGradGrad<double>       (*jac,  basisPhi,  basisRho , *gauss);

  // Formulations //
  // NB: FormulationJFLee is a friend of FormulationJFLee{One,...,Eight,} !
  //     So it can instanciate these classes...

  // Save FormulationJFLeeOne for update()
  formulationOne = new FormulationJFLeeOne(domain, *field, k, *proj);

  // Then push it in list
  fList.push_back(formulationOne);
  fList.push_back(new FormulationJFLeeTwo  (domain,  fPhi,*field,  k, *termPE));
  fList.push_back(new FormulationJFLeeThree(domain,  fRho,  fPhi, C2, *termRP));
  fList.push_back(new FormulationJFLeeFour (domain,  fPhi,         k, *termPP));
  fList.push_back(new FormulationJFLeeFive (domain, *field, fPhi,  k, *termEP));
  fList.push_back(new FormulationJFLeeSix  (domain, *field, fPhi, C1, *term22));
  fList.push_back(new FormulationJFLeeSeven(domain,  fRho           , *termRR));
  fList.push_back(new FormulationJFLeeEight(domain,  fPhi,  fRho    , *termPR));
}

FormulationJFLee::~FormulationJFLee(void){
  // Iterate & Delete Formulations //
  list<const FormulationBlock<Complex>*>::iterator end = fList.end();
  list<const FormulationBlock<Complex>*>::iterator it  = fList.begin();

  for(; it !=end; it++)
    delete *it;

  // Delete terms //
  delete proj;
  delete termPE;
  delete termRP;
  delete termPP;
  delete termEP;
  delete term22;
  delete termRR;
  delete termPR;

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
  basisE->preEvaluateFunctions(gC);

  // New RHS
  proj = new TermProjectionGrad<Complex>(*jac, *basisE, *gauss, *fieldG, ddm);

  // Update FormulationJFLeeOne (formulationOne):
  //                                          this FormulationBlock holds RHS
  formulationOne->update(*proj);
}
