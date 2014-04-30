#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"

#include "FormulationOSRC.h"
#include "FormulationUpdateOSRC.h"

using namespace std;

#include <set>
static
void initMap(const FunctionSpace& fs,
             const GroupOfElement& goe,
             map<Dof, Complex>& data){

  set<Dof> dSet;
  fs.getKeys(goe, dSet);

  set<Dof>::iterator it  = dSet.begin();
  set<Dof>::iterator end = dSet.end();

  for(; it != end; it++)
    data.insert(pair<Dof, Complex>(*it, 0));
}

static
void initMap(const vector<const FunctionSpaceScalar*>& fs,
             const GroupOfElement& goe,
             vector<map<Dof, Complex> >& data){

  const size_t size = data.size();

  set<Dof> dSet;
  set<Dof>::iterator end;
  set<Dof>::iterator it;

  if(size != fs.size())
    throw Exception("initMap: vector must have the same size");

  for(size_t i = 0; i < size; i++){
    dSet.clear();
    fs[i]->getKeys(goe, dSet);

    end = dSet.end();

    for(it = dSet.begin(); it != end; it++)
      data[i].insert(pair<Dof, Complex>(*it, 0));
  }
}

FormulationUpdateOSRC::FormulationUpdateOSRC(DDMContext& context){
  // Check if OSRC DDMContext //
  if(context.typeDDM != DDMContext::typeOSRC)
    throw Exception("FormulationUpdateOSRC needs a OSRC DDMContext");

  // Get Domain and FunctionSpace from DDMContext //
  ffspace = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOSRC needs a uniform mesh");

  /*
  // Sanity //
  for(int j = 0; j < NPade; j++)
    if(solPhi[j].size() != solU.size())
      throw Exception
        ("FormulationUpdateOSRC: solPhi[j] and solU should have the same size");

  if(ddm.size() != solU.size())
    throw Exception
      ("FormulationUpdateOSRC: ddm and solU should have the same size");
  */

  // Wavenumber //
  this->k = context.k;

  // Pade //
  int NPade = context.OSRC_NPade;
  A.resize(NPade);
  B.resize(NPade);

  C0 = FormulationOSRC::padeC0(NPade, M_PI / 4.);

  for(int j = 0; j < NPade; j++)
    A[j] = FormulationOSRC::padeAj(j + 1, NPade, M_PI / 4.);
  for(int j = 0; j < NPade; j++)
    B[j] = FormulationOSRC::padeBj(j + 1, NPade, M_PI / 4.);

  // Basis //
  const Basis& basis = ffspace->getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Jacobian //
  GroupOfJacobian jac(*ddomain, gC, "jacobian");

  // Get Volume Solution //
  map<Dof, Complex> solU;
  initMap(*ffspace, *ddomain, solU);
  context.system->getSolution(solU, 0);

  // Get Auxiliary Solutions //
  vector<map<Dof, Complex> > solPhi(NPade);
  initMap(*context.phi, *ddomain, solPhi);
  for(int i = 0; i < NPade; i++)
    context.system->getSolution(solPhi[i], 0);

  // UPhi[d] = sum_j (solU[d] - solPhi[j][d]) * A[j] / B[j] //

  // Init UPhi
  map<Dof, Complex> UPhi = solU;

  map<Dof, Complex>:: iterator UPhiEnd = UPhi.end();
  map<Dof, Complex>:: iterator UPhiIt;

  for(UPhiIt = UPhi.begin(); UPhiIt != UPhiEnd; UPhiIt++)
    UPhiIt->second = Complex(0, 0);

  // Iterator on solU and solPhi[j]
  map<Dof, Complex>::const_iterator solUIt;
  map<Dof, Complex>::const_iterator solPhiJIt;

  // Loop on j (Pade terms) and Degrees of freedom (iterators)
  for(int j = 0; j < NPade; j++)
    for(UPhiIt = UPhi.begin(),
          solUIt = solU.begin(), solPhiJIt = solPhi[j].begin();
        UPhiIt != UPhiEnd; UPhiIt++, solUIt++, solPhiJIt++)

      UPhiIt->second += (solUIt->second - solPhiJIt->second) * A[j] / B[j];

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Local Terms //
  lGout = new TermFieldField<double>(jac, basis, gauss);
  lGin  = new TermProjectionField<Complex>(jac, basis, gauss, *ffspace, ddm);
  lC0   = new TermProjectionField<Complex>(jac, basis, gauss, *ffspace, solU);
  lAB   = new TermProjectionField<Complex>(jac, basis, gauss, *ffspace, UPhi);
}

FormulationUpdateOSRC::~FormulationUpdateOSRC(void){
  delete lGout;
  delete lGin;
  delete lC0;
  delete lAB;
}

Complex FormulationUpdateOSRC::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return
    Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateOSRC::rhs(size_t equationI, size_t elementId) const{
  return
    Complex(-1,  0    ) *      lGin->getTerm(0, equationI, elementId) +
    Complex( 0, -2 * k) * C0 * lC0->getTerm(0, equationI, elementId)  +
    Complex( 0, -2 * k) *      lAB->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationUpdateOSRC::field(void) const{
  return *ffspace;
}

const FunctionSpace& FormulationUpdateOSRC::test(void) const{
  return *ffspace;
}

const GroupOfElement& FormulationUpdateOSRC::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateOSRC::isBlock(void) const{
  return true;
}
