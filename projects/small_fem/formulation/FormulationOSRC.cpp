#include <cmath>

#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"

#include "FormulationOSRCOne.h"
#include "FormulationOSRCTwo.h"
#include "FormulationOSRCThree.h"
#include "FormulationOSRCFour.h"
#include "FormulationOSRC.h"

using namespace std;

FormulationOSRC::FormulationOSRC(const GroupOfElement& domain,
                                 const FunctionSpaceScalar& field,
                                 const vector<const FunctionSpaceScalar*>& aux,
                                 double k,
                                 Complex keps,
                                 int NPade,
                                 const map<Dof, Complex>& ddmDof){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRC needs a uniform mesh");

  // Get Basis (same for field and auxiliary function spaces) //
  const Basis& basis = field.getBasis(eType);
  const size_t order = basis.getOrder();

  // Gaussian Quadrature //
  Quadrature gaussFF(eType, order    , 2);
  Quadrature gaussGG(eType, order - 1, 2);

  const fullMatrix<double>& gCFF = gaussFF.getPoints();
  const fullMatrix<double>& gCGG = gaussGG.getPoints();

  // Local Terms //
  basis.preEvaluateFunctions(gCFF);
  basis.preEvaluateDerivatives(gCGG);

  GroupOfJacobian jacFF(domain, gCFF, "jacobian");
  GroupOfJacobian jacGG(domain, gCGG, "invert");

  // NB: Since the Formulations share the same basis functions,
  //     the local terms will be the same !
  //     It's the Dof numbering imposed by the function spaces that will differ
  localFF = new TermFieldField<double>(jacFF, basis, gaussFF);
  localGG = new TermGradGrad<double>(jacGG, basis, gaussGG);
  localPr =
    new TermProjectionField<Complex>(jacFF, basis, gaussFF, field, ddmDof);

  // Formulations //
  // NB: FormulationOSRC is a friend of FormulationOSRC{One,Two,Three,Four,} !
  //     So it can instanciate those classes...

  fList.push_back
    (new FormulationOSRCOne
     (domain, field, k, NPade, *localFF, *localPr));              // u.u'

  // Loop on phi[j]
  for(int j = 0; j < NPade; j++){
    fList.push_back
      (new FormulationOSRCTwo
       (domain, *aux[j], field, k, keps, NPade, j+1, *localGG)); // phi[j].u'

    fList.push_back
      (new FormulationOSRCThree
       (domain, *aux[j], keps, NPade, j+1, *localFF, *localGG)); // ph[j].ph[j]'

    fList.push_back
      (new FormulationOSRCFour
       (domain, field, *aux[j], *localFF));                      // u.phi[j]'
  }
}

FormulationOSRC::~FormulationOSRC(void){
  // Iterate & Delete Formulations //
  list<const Formulation<Complex>*>::iterator end = fList.end();
  list<const Formulation<Complex>*>::iterator it  = fList.begin();

  for(; it !=end; it++)
    delete *it;

  // Delete terms //
  delete localFF;
  delete localGG;
  delete localPr;
}

const list<const Formulation<Complex>*>&
FormulationOSRC::getFormulations(void) const{
  return fList;
}

double FormulationOSRC::pade_aj(int j, int N){
  double tmp = sin((double)j * M_PI / (2. * N + 1.));

  return 2. / (2. * N + 1.) * tmp * tmp;
}

double FormulationOSRC::pade_bj(int j, int N){
  double tmp = cos((double)j * M_PI / (2. *N + 1.));

  return tmp * tmp;
}

Complex FormulationOSRC::padeC0(int N, double theta){
  Complex sum = Complex(1, 0);
  Complex one = Complex(1, 0);
  Complex z   = Complex(cos(-theta) - 1,  sin(-theta));

  for(int j = 1; j <= N; j++)
    sum += (z * pade_aj(j, N)) / (one + z * pade_bj(j, N));

  z = Complex(cos(theta / 2.), sin(theta / 2.));

  return sum * z;
}

Complex FormulationOSRC::padeAj(int j, int N, double theta){
  Complex one = Complex(1, 0);
  Complex res;
  Complex z;

  z   = Complex(cos(-theta / 2.), sin(-theta / 2.));
  res = z * pade_aj(j, N);

  z   = Complex(cos(-theta) - 1., sin(-theta));
  res = res / ((one + z * pade_bj(j, N)) * (one + z * pade_bj(j, N)));

  return res;
}

Complex FormulationOSRC::padeBj(int j, int N, double theta){
  Complex one = Complex(1, 0);
  Complex res;
  Complex z;

  z   = Complex(cos(-theta), sin(-theta));
  res = z * pade_bj(j, N);

  z   = Complex(cos(-theta) - 1., sin(-theta));
  res = res / (one + z * pade_bj(j, N));

  return res;
}
