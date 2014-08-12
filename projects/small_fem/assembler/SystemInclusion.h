///////////////////////////////////////////
// Templates Implementations for System: //
// Inclusion compilation model           //
//                                       //
// Damn you gcc: we want 'export' !      //
///////////////////////////////////////////

template<typename scalar>
System<scalar>::System(void){
  // Init //
  A = NULL;
  b = NULL;
  x = NULL;

  // Dof Manager //
  this->dofM = new DofManager<scalar>;

  // Per-thread non-zero term //
  #pragma omp parallel
  {
    #pragma omp master
    this->nNZCount.resize(omp_get_num_threads());
  }

  for(size_t i = 0; i < this->nNZCount.size(); i++)
    this->nNZCount[i] = 0;

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
void System<scalar>::assemble(void){
  // Enumerate Dofs in DofManager //
  this->dofM->generateGlobalIdSpace();

  // Formulations Iterators //
  typename std::list<const FormulationBlock<scalar>*>::iterator it;
  typename std::list<const FormulationBlock<scalar>*>::iterator end;

  // Count FE terms //
  it  = this->formulation.begin();
  end = this->formulation.end();

  for(; it != end; it++){
    // Get All Dofs (Field & Test) per Element
    const std::vector<std::vector<Dof> >& dofField =
      (*it)->field().getKeys((*it)->domain());
    const std::vector<std::vector<Dof> >& dofTest  =
      (*it)->test().getKeys((*it)->domain());

    // Count
    const size_t E = dofField.size(); // Should be equal to dofTest.size().?.

    #pragma omp parallel for
    for(size_t i = 0; i < E; i++)
      this->nNZCount[omp_get_thread_num()] =
        SystemAbstract<scalar>::countTerms(this->nNZCount[omp_get_thread_num()],
                                           i, dofField[i], dofTest[i], **it);
  }

  // Alloc //
  const size_t size = this->dofM->getLocalSize();

  A = new SolverMatrix<scalar>(size, size, this->nNZCount);
  b = new SolverVector<scalar>(size);

  // Assemble //
  it  = this->formulation.begin();
  end = this->formulation.end();

  for(; it != end; it++){
    // Get All Dofs (Field & Test) per Element
    const std::vector<std::vector<Dof> >& dofField =
      (*it)->field().getKeys((*it)->domain());
    const std::vector<std::vector<Dof> >& dofTest  =
      (*it)->test().getKeys((*it)->domain());

    // Assemble
    const size_t E = dofField.size(); // Should be equal to dofTest.size().?.

    #pragma omp parallel for
    for(size_t i = 0; i < E; i++)
      SystemAbstract<scalar>::
        assemble(*A, *b, i, dofField[i], dofTest[i], **it);
  }

  // The system is assembled //
  this->assembled = true;
}

template<typename scalar>
void System<scalar>::solve(void){
  // Is the System assembled ? //
  if(!this->assembled)
    assemble();

  // Use SolverMUMPS //
  x = new fullVector<scalar>;     // RHS memory
  this->solver.solve(*A, *b, *x); // Solve

  // System solved ! //
  this->solved = true;
}

template<typename scalar>
void System<scalar>::assembleAgainRHS(void){
  // Is the full system assembled ? //
  if(!this->assembled)
    throw Exception
      ("System::assembleAgainRHS() needs a first call to System::assemble()");

  // Set RHS to zero //
  this->b->reset();

  // Iterate on Formulations //
  typename std::list<const FormulationBlock<scalar>*>::iterator it;
  typename std::list<const FormulationBlock<scalar>*>::iterator end;

  it  = this->formulation.begin();
  end = this->formulation.end();

  for(; it != end; it++){
    // Get All Dofs (Test only) per Element
    const std::vector<std::vector<Dof> >& dofTest =
      (*it)->test().getKeys((*it)->domain());

    // Assemble
    const size_t E = dofTest.size();

    #pragma omp parallel for
    for(size_t i = 0; i < E; i++)
      SystemAbstract<scalar>::assembleRHSOnly(*b, i, dofTest[i], **it);
  }
}

template<typename scalar>
void System<scalar>::solveAgain(void){
  // Is the full system solved ? //
  if(!this->solved)
    throw Exception
      ("System::solveAgain() needs a first call to System::solve()");

  this->solver.setRHS(*b);
  this->solver.solve(*x);
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
void System<scalar>::getSolution(FEMSolution<scalar>& feSol,
                                 const FunctionSpace& fs,
                                 const GroupOfElement& domain) const{
  // Solved ?
  if(!this->solved)
    throw Exception("System: addSolution -- System not solved");

  // Coefficients //
  // Get Dofs
  std::set<Dof> dof;
  fs.getKeys(domain, dof);

  // Get Coefficient
  const std::set<Dof>::iterator  end = dof.end();
  std::set<Dof>::iterator        it  = dof.begin();
  typename std::map<Dof, scalar> coef;

  for(; it != end; it++)
    coef.insert(std::pair<Dof, scalar>(*it, 0));

  // Populate Map
  getSolution(coef, 0);

  // FEMSolution //
  feSol.addCoefficients(0, 0, domain, fs, coef);
}

template<typename scalar>
void System<scalar>::
getSolution(FEMSolution<scalar>& feSol,
            const FunctionSpace& fs,
            const std::vector<const GroupOfElement*>& domain) const{
  // Get size
  const size_t size = domain.size();

  // Get solution for each domain
  for(size_t i = 0; i < size; i++)
    getSolution(feSol, fs, *domain[i]);
}

template<typename scalar>
void System<scalar>::writeMatrix(std::string fileName,
                                 std::string matrixName) const{
  A->writeToMatlabFile(fileName, matrixName);
}
