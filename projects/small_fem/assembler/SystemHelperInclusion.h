/////////////////////////////////////////////////
// Templates Implementations for SystemHelper: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include "System.h"

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
          const FunctionSpace& fs,
          const GroupOfElement& goe,
          scalar (*f)(fullVector<double>& xyz)){

  // Solve Projection //
  const FunctionSpaceScalar& fsS = static_cast<const FunctionSpaceScalar&>(fs);
  FormulationProjectionScalar<scalar> formulation(goe, fsS, f);

  System<scalar> projection;
  projection.addFormulation(formulation);
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
          const FunctionSpace& fs,
          const GroupOfElement& goe,
          fullVector<scalar> (*f)(fullVector<double>& xyz)){

  // Solve Projection //
  const FunctionSpaceVector& fsV = static_cast<const FunctionSpaceVector&>(fs);
  FormulationProjectionVector<scalar> formulation(goe, fsV, f);

  System<scalar> projection;
  projection.addFormulation(formulation);
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
