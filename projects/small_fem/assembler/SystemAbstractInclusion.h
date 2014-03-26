///////////////////////////////////////////////////
// Templates Implementations for SystemAbstract: //
// Inclusion compilation model                   //
//                                               //
// Damn you gcc: we want 'export' !              //
///////////////////////////////////////////////////

template<typename scalar>
SystemAbstract<scalar>::~SystemAbstract(void){
}

template<typename scalar>
bool SystemAbstract<scalar>::isAssembled(void) const{
  return assembled;
}

template<typename scalar>
bool SystemAbstract<scalar>::isSolved(void) const{
  return solved;
}

template<typename scalar>
size_t SystemAbstract<scalar>::getSize(void) const{
  return dofM.getUnfixedDofNumber();
}

template<typename scalar>
void SystemAbstract<scalar>::
addFormulation(const Formulation<scalar>& formulation){
  // Add formulation in list of formulation //
  this->formulation.push_back(&formulation);

  // Get Formulation Dofs (Field & Test) //
  std::set<Dof> dofField;
  std::set<Dof> dofTest;
  formulation.field().getKeys(formulation.domain(), dofField);
  formulation.test().getKeys(formulation.domain(), dofTest);

  // Add them to DofManager //
  this->dofM.addToDofManager(dofField);
  this->dofM.addToDofManager(dofTest);
}

template<typename scalar>
void SystemAbstract<scalar>::
addFormulation(const CoupledFormulation<scalar>& formulation){
  // Get the list of Formulations //
  const std::list<const Formulation<scalar>*>&
    fList = formulation.getFormulations();

  // Iterate on list and add the pointed Formulations //
  typename std::list<const Formulation<scalar>*>::const_iterator
    end = fList.end();

  typename std::list<const Formulation<scalar>*>::const_iterator
    it  = fList.begin();

  for(; it != end; it++)
    this->addFormulation(**it);
}

template<typename scalar>
void SystemAbstract<scalar>::constraint(const std::map<Dof, scalar>& constr){
  typename std::map<Dof, scalar>::const_iterator it  = constr.begin();
  typename std::map<Dof, scalar>::const_iterator end = constr.end();

  for(; it != end; it++)
    dofM.fixValue(it->first, it->second);
}

template<typename scalar>
void SystemAbstract<scalar>::
assemble(SolverMatrix<scalar>& A,
         SolverVector<scalar>& b,
         size_t elementId,
         const std::vector<Dof>& dofField,
         const std::vector<Dof>& dofTest,
         formulationPtr& term,
         const Formulation<scalar>& formulation){

  const size_t N = dofTest.size();
  const size_t M = dofField.size();

  size_t dofI;
  size_t dofJ;

  for(size_t i = 0; i < N; i++){
    dofI = dofM.getGlobalId(dofTest[i]);

    // If not a fixed Dof line: assemble
    if(dofI != DofManager<scalar>::isFixedId()){
      for(size_t j = 0; j < M; j++){
        dofJ = dofM.getGlobalId(dofField[j]);

        // If not a fixed Dof
        if(dofJ != DofManager<scalar>::isFixedId())
          A.add(dofI, dofJ, (formulation.*term)(i, j, elementId));

        // If fixed Dof (for column 'dofJ'):
        //    add to right hand side (with a minus sign) !
        else
          b.add(dofI,
                minusSign * dofM.getValue(dofField[j]) *
                           (formulation.*term)(i, j, elementId));
      }

      b.add(dofI, formulation.rhs(i, elementId));
    }
  }
}
