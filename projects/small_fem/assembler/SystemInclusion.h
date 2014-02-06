///////////////////////////////////////////
// Templates Implementations for System: //
// Inclusion compilation model           //
//                                       //
// Damn you gcc: we want 'export' !      //
///////////////////////////////////////////

#include "SolverMUMPS.h"

template<typename scalar>
System<scalar>::System(const Formulation<scalar>& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fsField());

  // Get Formulation Dofs //
  std::set<Dof> dof;
  this->fs->getKeys(formulation.domain(), dof);

  // Get Dof Manager //
  this->dofM = new DofManager<scalar>();
  this->dofM->addToDofManager(dof);

  // Init //
  A = NULL;
  b = NULL;
  x = NULL;

  // The system is not assembled and not solved //
  this->assembled = false;
  this->solved    = false;
}

template<typename scalar>
System<scalar>::~System(void){
  delete this->dofM;

  if(A)
    delete A;

  if(b)
    delete b;

  if(x)
    delete x;
}

template<typename scalar>
void System<scalar>::addBorderTerm(const Formulation<scalar>& formulation){
  // Get the FunctionSpace of the given formulation
  const FunctionSpace& fs = formulation.fsField();

  // Get All Dofs per Element //
  std::vector<std::vector<Dof> > dof;
  fs.getKeys(formulation.domain(), dof);

  // Get Formulation Term //
  typename SystemAbstract<scalar>::formulationPtr term =
    &Formulation<scalar>::weak;

  // Assemble //
  const size_t E = dof.size();

  #pragma omp parallel for
  for(size_t i = 0; i < E; i++)
    SystemAbstract<scalar>::
      assemble(*A, *b, i, dof[i], term, formulation);
}

template<typename scalar>
void System<scalar>::assemble(void){
  // Enumerate //
  this->dofM->generateGlobalIdSpace();

  // Get All Dofs per Element //
  std::vector<std::vector<Dof> > dof;
  this->fs->getKeys(this->formulation->domain(), dof);

  // Get Formulation Term //
  typename SystemAbstract<scalar>::formulationPtr term =
    &Formulation<scalar>::weak;

  // Alloc //
  const size_t size = this->dofM->getUnfixedDofNumber();

  A = new SolverMatrix<scalar>(size, size);
  b = new SolverVector<scalar>(size);

  // Assemble //
  const size_t E = dof.size();

  #pragma omp parallel for
  for(size_t i = 0; i < E; i++)
    SystemAbstract<scalar>::
      assemble(*A, *b, i, dof[i], term, *this->formulation);

  // The system is assembled //
  this->assembled = true;
}

template<typename scalar>
void System<scalar>::solve(void){
  // Is the System assembled ? //
  if(!this->assembled)
    assemble();

  // Use SolverMUMPS //
  SolverMUMPS<scalar> solver;
  x = new fullVector<scalar>;

  solver.solve(*A, *b, *x);

  // System solved ! //
  this->solved = true;
}

template<typename scalar>
size_t System<scalar>::getNComputedSolution(void) const{
  return 1;
}

template<typename scalar>
void System<scalar>::getSolution(fullVector<scalar>& sol, size_t nSol) const{
  sol.setAsProxy(*x, 0, x->size());
}

template<typename scalar>
void System<scalar>::getSolution(std::map<Dof, scalar>& sol, size_t nSol) const{
  // Get All Dofs
  typename std::map<Dof, scalar>::iterator it  = sol.begin();
  typename std::map<Dof, scalar>::iterator end = sol.end();

  // Loop on Dofs and set Values
  for(; it != end; it++){
    size_t gId = this->dofM->getGlobalId(it->first);

    if(gId == DofManager<scalar>::isFixedId())
      it->second = this->dofM->getValue(it->first);

    else
      it->second = (*x)(gId);
  }
}

template<typename scalar>
void System<scalar>::getSolution(FEMSolution<scalar>& feSol) const{
  // Solved ?
  if(!this->solved)
    throw Exception("System: addSolution -- System not solved");

  // Coefficients //
  // FunctionSpace & Domain
  const FunctionSpace&  fs  = this->formulation->fsField();
  const GroupOfElement& goe = this->formulation->domain();

  // Get Dofs
  std::set<Dof> dof;
  fs.getKeys(goe, dof);

  // Get Coefficient
  const std::set<Dof>::iterator  end = dof.end();
  std::set<Dof>::iterator        it  = dof.begin();
  typename std::map<Dof, scalar> coef;

  for(; it != end; it++)
    coef.insert(std::pair<Dof, scalar>(*it, 0));

  // Populate Map
  getSolution(coef, 0);

  // FEMSolution //
  feSol.addCoefficients(0, 0, goe, fs, coef);
}

template<typename scalar>
void System<scalar>::writeMatrix(std::string fileName,
                                 std::string matrixName) const{
  A->writeToMatlabFile(fileName, matrixName);
}
