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
          const FunctionSpaceVector& fs,
          const GroupOfElement& goe,
          fullVector<scalar> (*f)(fullVector<double>& xyz)){

  // Solve Projection //
  FormulationProjectionVector<scalar> formulation(goe, fs, f);

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
