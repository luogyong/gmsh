#include <cmath>

#include "Exception.h"
#include "FormulationOSRCScalarOne.h"
#include "FormulationOSRCScalarTwo.h"
#include "FormulationOSRCScalarThree.h"
#include "FormulationOSRCScalarFour.h"
#include "FormulationOSRCScalar.h"

using namespace std;

FormulationOSRCScalar::FormulationOSRCScalar(DDMContextOSRC& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and auxiliary FunctionSpaces from DDMContext //
  const GroupOfElement&                     domain = context.getDomain();
  const vector<const FunctionSpaceScalar*>& aux    =
    context.getAuxFunctionSpace();

  // Save field FunctionSpace
  field = &context.getFunctionSpace(); // Saved from update()

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRCScalar needs a uniform mesh");

  // Get Basis (same for field and auxiliary function spaces) //
  basis = &field->getBasis(eType); // Saved from update()
  const size_t order = basis->getOrder();

  // k, keps and NPade //
  double k     = context.getWavenumber();
  Complex keps = context.getComplexWavenumber();
  int NPade    = context.getNPade();

  // Gaussian Quadrature //
  gaussFF = new Quadrature(eType, order, 2); // Saved from update()
  Quadrature gaussGG(eType, order - 1, 2);

  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  const fullMatrix<double>& gCGG = gaussGG.getPoints();

  // Local Terms //
  basis->preEvaluateFunctions(gCFF);
  basis->preEvaluateDerivatives(gCGG);

  jacFF = new GroupOfJacobian(domain, gCFF, "jacobian"); // Saved from update()
  GroupOfJacobian jacGG(domain, gCGG, "invert");

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // NB: Since the Formulations share the same basis functions,
  //     the local terms will be the same !
  //     It's the Dof numbering imposed by the function spaces that will differ
  localFF = new TermFieldField<double>(*jacFF, *basis, *gaussFF);
  localGG = new TermGradGrad<double>(jacGG, *basis, gaussGG);
  localPr =
    new TermProjectionField<Complex>(*jacFF, *basis, *gaussFF, *field, ddm);

  // Formulations //
  // NB: FormulationOSRCScalar is a friend
  //     of FormulationOSRCScalar{One,Two,Three,Four,} !
  //     So it can instanciate those classes...

  // Save FormulationOSRCScalarOne for update()
  formulationOne = new FormulationOSRCScalarOne
    (domain, *field, k, NPade, *localFF, *localPr);                // u.u'

  // Then push it in list
  fList.push_back(formulationOne);

  // Loop on phi[j]
  for(int j = 0; j < NPade; j++){
    fList.push_back
      (new FormulationOSRCScalarTwo
       (domain, *aux[j], *field, k, keps, NPade, j+1, *localGG)); // phi[j].u'

    fList.push_back
      (new FormulationOSRCScalarThree
       (domain, *aux[j], keps, NPade, j+1, *localFF, *localGG)); // ph[j].ph[j]'

    fList.push_back
      (new FormulationOSRCScalarFour
       (domain, *field, *aux[j], *localFF));                     // u.phi[j]'
  }
}

FormulationOSRCScalar::~FormulationOSRCScalar(void){
  // Iterate & Delete Formulations //
  list<const FormulationBlock<Complex>*>::iterator end = fList.end();
  list<const FormulationBlock<Complex>*>::iterator it  = fList.begin();

  for(; it !=end; it++)
    delete *it;

  // Delete terms //
  delete localFF;
  delete localGG;
  delete localPr;

  // Delete update stuffs //
  delete jacFF;
  delete gaussFF;
}

const list<const FormulationBlock<Complex>*>&
FormulationOSRCScalar::getFormulationBlocks(void) const{
  return fList;
}

bool FormulationOSRCScalar::isBlock(void) const{
  return false;
}

void FormulationOSRCScalar::update(void){
  // Delete RHS (localPr)
  delete localPr;

  // Get DDM Dofs from DDMContext
  const map<Dof, Complex>& ddm = context->getDDMDofs();

  // Pre-evalution
  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  basis->preEvaluateFunctions(gCFF);

  // New RHS
  localPr =
    new TermProjectionField<Complex>(*jacFF, *basis, *gaussFF, *field, ddm);

  // Update FormulationOSRCScalarOne (formulationOne):
  //                                             this FormulationBlock holds RHS
  formulationOne->update(*localPr);
}
