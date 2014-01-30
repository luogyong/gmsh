/////////////////////////////////////////////////
// Templates Implementations for SystemHelper: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include "System.h"
#include "BasisGenerator.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"


template<typename scalar>
SystemHelper<scalar>::SystemHelper(void){
}

template<typename scalar>
SystemHelper<scalar>::~SystemHelper(void){
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          GroupOfElement& goe,
          scalar (*f)(fullVector<double>& xyz)){

  // Check if homogene GoE and get geo type //
  const std::vector<size_t>& gType = goe.getTypeStats();
  const size_t nGType = gType.size();
  size_t eType = (size_t)(-1);

  for(size_t i = 0; i < nGType; i++)
    if((eType == (size_t)(-1)) && (gType[i] != 0))
      eType = i;
    else if((eType != (size_t)(-1)) && (gType[i] != 0))
      throw Exception("SystemHelper::dirichlet needs a uniform mesh");

  // Get this SystemAbstract FunctionSpace & its Basis type //
  const FunctionSpace& fs = sys.getFunctionSpace();

  const std::vector<size_t>& fsGType = fs.getSupport().getTypeStats();
  size_t fsType = (size_t)(-1);

  for(size_t i = 0; i < nGType; i++)
    if((fsType == (size_t)(-1)) && (fsGType[i] != 0))
      fsType = i;

  // Get Function Space for Projection (formFS) //
  Basis* basis = BasisGenerator::generate(eType,
                                          fs.getBasis(fsType).getForm(),
                                          fs.getBasis(fsType).getOrder(),
                                          "hierarchical");

  FunctionSpaceScalar formFS(goe, *basis);

  // Solve Projection //
  FormulationProjectionScalar<scalar> form(f, formFS);

  System<scalar> projection(form);
  projection.assemble();
  projection.solve();

  // Map of Dofs //
  const std::set<Dof>& dof = formFS.getAllDofs();
  std::set<Dof>::iterator it  = dof.begin();
  std::set<Dof>::iterator end = dof.end();

  std::map<Dof, scalar> constr;
  for(; it != end; it++)
    constr.insert(std::pair<Dof, scalar>(*it, 0));

  // Get Solution and Dirichlet Constraint //
  projection.getSolution(constr, 0);
  sys.constraint(constr);

  delete basis;
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          GroupOfElement& goe,
          fullVector<scalar> (*f)(fullVector<double>& xyz)){

  // Check if homogene GoE and get geo type //
  const std::vector<size_t>& gType = goe.getTypeStats();
  const size_t nGType = gType.size();
  size_t eType = (size_t)(-1);

  for(size_t i = 0; i < nGType; i++)
    if((eType == (size_t)(-1)) && (gType[i] != 0))
      eType = i;
    else if((eType != (size_t)(-1)) && (gType[i] != 0))
      throw Exception("SystemHelper::dirichlet needs a uniform mesh");

  // Get this SystemAbstract FunctionSpace & its Basis type //
  const FunctionSpace& fs = sys.getFunctionSpace();

  const std::vector<size_t>& fsGType = fs.getSupport().getTypeStats();
  size_t fsType = (size_t)(-1);

  for(size_t i = 0; i < nGType; i++)
    if((fsType == (size_t)(-1)) && (fsGType[i] != 0))
      fsType = i;

  // Get Function Space for Projection (formFS) //
  Basis* basis = BasisGenerator::generate(eType,
                                          fs.getBasis(fsType).getForm(),
                                          fs.getBasis(fsType).getOrder(),
                                          "hierarchical");

  FunctionSpaceVector formFS(goe, *basis);

  // Solve Projection //
  FormulationProjectionVector<scalar> form(f, formFS);

  System<scalar> projection(form);
  projection.assemble();
  projection.solve();

  // Map of Dofs //
  const std::set<Dof>& dof = formFS.getAllDofs();
  std::set<Dof>::iterator it  = dof.begin();
  std::set<Dof>::iterator end = dof.end();

  std::map<Dof, scalar> constr;
  for(; it != end; it++)
    constr.insert(std::pair<Dof, scalar>(*it, 0));

  // Get Solution and Dirichlet Constraint //
  projection.getSolution(constr, 0);
  sys.constraint(constr);

  delete basis;
}
