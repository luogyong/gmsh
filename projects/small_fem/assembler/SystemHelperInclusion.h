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

  // Get this SystemAbstract FunctionSpace //
  const FunctionSpace& fs = sys.getFunctionSpace();

  // Get Function Space for Projection (formFS) //
  FunctionSpaceScalar formFS(goe, fs.getOrder());

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

  // Get this SystemAbstract FunctionSpace //
  const FunctionSpace& fs = sys.getFunctionSpace();

  // Get Function Space for Projection (formFS) //
  FunctionSpaceVector formFS(goe, fs.getOrder());

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
}
