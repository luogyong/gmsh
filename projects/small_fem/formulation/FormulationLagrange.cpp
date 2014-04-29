#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"

#include "FormulationLagrangeTwo.h"
#include "FormulationLagrangeOne.h"
#include "FormulationLagrange.h"


using namespace std;

FormulationLagrange::FormulationLagrange(const GroupOfElement& domain,
                                         const FunctionSpaceScalar& field,
                                         const FunctionSpaceScalar& lagrange,
                                         double (*f)(fullVector<double>& xyz)){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationLagrange needs a uniform mesh");

  // Get Basis (same for field and Lagrange function spaces) //
  const Basis& basis = lagrange.getBasis(eType);
  const size_t order = basis.getOrder();

  // Gaussian Quadrature //
  Quadrature gauss(eType, order, 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Local Terms //
  basis.preEvaluateFunctions(gC);

  GroupOfJacobian jac(domain, gC, "jacobian");

  // NB: Since the two Formulations (Lagrand . Field & Field . Lagrage)
  //     share the same basis functions, the local terms will be the same !
  //     It's the Dof numbering imposed by the function spaces that will differ
  local = new TermFieldField<double>     (jac, basis, gauss);
  proj  = new TermProjectionField<double>(jac, basis, gauss, f);

  // Formulations //
  // NB: FormulationLagrange is a friend of FormulationLagrange{One,Two,} !
  //     So it can instanciate those classes...

  fList.push_back
    (new FormulationLagrangeOne(domain, field, lagrange, *local, *proj));
  fList.push_back
    (new FormulationLagrangeTwo(domain, lagrange, field, *local));
}

FormulationLagrange::~FormulationLagrange(void){
  // Iterate & Delete Formulations //
  list<const FormulationBlock<Complex>*>::iterator end = fList.end();
  list<const FormulationBlock<Complex>*>::iterator it  = fList.begin();

  for(; it !=end; it++)
    delete *it;

  // Delete terms //
  delete local;
  delete proj;
}

const list<const FormulationBlock<Complex>*>&
FormulationLagrange::getFormulationBlocks(void) const{
  return fList;
}

bool FormulationLagrange::isBlock(void) const{
  return false;
}
