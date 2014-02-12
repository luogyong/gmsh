/////////////////////////////////////////////////
// Templates Implementations for SystemHelper: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include "System.h"

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
          const GroupOfElement& domain,
          scalar (*f)(fullVector<double>& xyz)){

  // Scalar Formulation //
  FormulationProjection<scalar> formulation(domain, fs, f);

  // Dirichlet
  dirichlet(sys, fs, domain, formulation);
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          const FunctionSpace& fs,
          const GroupOfElement& domain,
          fullVector<scalar> (*f)(fullVector<double>& xyz)){

  // Vector Formulation //
  FormulationProjection<scalar> formulation(domain, fs, f);

  // Dirichlet
  dirichlet(sys, fs, domain, formulation);
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          const FunctionSpace& fs,
          const GroupOfElement& domain,
          FormulationProjection<scalar>& formulation){

  // Solve Projection //
  System<scalar> projection;

  projection.addFormulation(formulation);
  projection.assemble();
  projection.solve();

  // Map of Dofs //
  std::set<Dof> dof;
  std::map<Dof, scalar> constraint;

  fs.getKeys(domain, dof);
  std::set<Dof>::iterator it  = dof.begin();
  std::set<Dof>::iterator end = dof.end();

  for(; it != end; it++)
    constraint.insert(std::pair<Dof, scalar>(*it, 0));

  // Get Solution and Dirichlet Constraint //
  projection.getSolution(constraint, 0);
  sys.constraint(constraint);
}
