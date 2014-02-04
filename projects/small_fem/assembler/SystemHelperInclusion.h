/////////////////////////////////////////////////
// Templates Implementations for SystemHelper: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include "System.h"
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
          const GroupOfElement& goe,
          scalar (*f)(fullVector<double>& xyz)){

  // Check if homogene GoE and get geo type //
  if(!goe.isUniform().first)
    throw Exception("SystemHelper::dirichlet needs a uniform mesh");

  // Get this SystemAbstract FunctionSpace //
  const FunctionSpace& fs = sys.getFunctionSpace();

  // Get Function Space for Projection (formFS) //
  FunctionSpaceScalar formFS(goe, fs.getOrder());

  // Solve Projection //
  FormulationProjectionScalar<scalar> form(goe, formFS, f);

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
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          const FunctionSpaceScalar& fs,
          const GroupOfElement& goe,
          scalar (*f)(fullVector<double>& xyz)){

  // Solve Projection //
  FormulationProjectionScalar<scalar> formulation(goe, fs, f);

  System<scalar> projection(formulation);
  projection.assemble();
  projection.solve();

  // Map of Dofs //
  std::set<Dof> dof;
  std::map<Dof, scalar> constraint;

  fs.getKeys(goe, dof);
  std::set<Dof>::iterator it  = dof.begin();
  std::set<Dof>::iterator end = dof.end();

  for(; it != end; it++)
    constraint.insert(std::pair<Dof, scalar>(*it, 0));

  // Get Solution and Dirichlet Constraint //
  projection.getSolution(constraint, 0);
  sys.constraint(constraint);
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          const GroupOfElement& goe,
          fullVector<scalar> (*f)(fullVector<double>& xyz)){

  // Check if homogene GoE and get geo type //
  if(!goe.isUniform().first)
    throw Exception("SystemHelper::dirichlet needs a uniform mesh");

  // Get this SystemAbstract FunctionSpace //
  const FunctionSpace& fs = sys.getFunctionSpace();

  // Get Function Space for Projection (formFS) //
  FunctionSpaceVector formFS(goe, fs.getOrder());

  // Solve Projection //
  FormulationProjectionVector<scalar> form(goe, formFS, f);

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
}
