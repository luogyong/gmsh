///////////////////////////////////////////////////
// Templates Implementations for SystemAbstract: //
// Inclusion compilation model                   //
//                                               //
// Damn you gcc: we want 'export' !              //
///////////////////////////////////////////////////

#include <cmath>

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
  if(isAssembled())
    return dofM->getGlobalSize();

  else
    throw
      Exception("SystemAbstract::getSize() cannot be called before assembly");
}

template<typename scalar>
void SystemAbstract<scalar>::
addFormulation(const Formulation<scalar>& formulation){
  // Add FormulationBlocks to list (this->formulation) //

  // Is this a FormulationBlock ?
  if(formulation.isBlock())
    addFormulationBlock
      (static_cast<const FormulationBlock<scalar>&>(formulation),
       this->formulation);

  else
    addFormulationCoupled
      (static_cast<const FormulationCoupled<scalar>&>(formulation),
       this->formulation);
}

template<typename scalar>
void SystemAbstract<scalar>::
addFormulationBlock(const FormulationBlock<scalar>& formulation,
                    std::list<const FormulationBlock<scalar>*>& fList){
  // Add formulation block in the given list //
  fList.push_back(&formulation);

  // Get Formulation Dofs (Field & Test) //
  const std::vector<std::vector<Dof> >& dofField =
    formulation.field().getKeys(formulation.domain());

  const std::vector<std::vector<Dof> >& dofTest =
    formulation.test().getKeys(formulation.domain());

  // Add them to DofManager //
  this->dofM->addToDofManager(dofField);
  this->dofM->addToDofManager(dofTest);
}

template<typename scalar>
void SystemAbstract<scalar>::
addFormulationCoupled(const FormulationCoupled<scalar>& formulation,
                      std::list<const FormulationBlock<scalar>*>& fList){
  // Get the list of Formulation Blocks //
  const std::list<const FormulationBlock<scalar>*>&
    blockList = formulation.getFormulationBlocks();

  // Iterate on list and add the pointed Formulations //
  typename std::list<const FormulationBlock<scalar>*>::const_iterator
    end = blockList.end();

  typename std::list<const FormulationBlock<scalar>*>::const_iterator
    it  = blockList.begin();

  for(; it != end; it++)
    this->addFormulationBlock(**it, fList);
}

template<typename scalar>
void SystemAbstract<scalar>::constraint(const std::map<Dof, scalar>& constr){
  typename std::map<Dof, scalar>::const_iterator it  = constr.begin();
  typename std::map<Dof, scalar>::const_iterator end = constr.end();

  for(; it != end; it++)
    dofM->fixValue(it->first, it->second);
}

template<typename scalar>
size_t SystemAbstract<scalar>::
countTerms(size_t offset,
           size_t elementId,
           const std::vector<Dof>& dofField,
           const std::vector<Dof>& dofTest,
           const FormulationBlock<scalar>& formulation){

  const size_t N = dofTest.size();
  const size_t M = dofField.size();

  size_t count = offset;
  size_t dofI;
  size_t dofJ;

  for(size_t i = 0; i < N; i++){
    dofI = dofM->getGlobalId(dofTest[i]);

    // If not a fixed Dof line
    if(dofI != DofManager<scalar>::isFixedId()){
      for(size_t j = 0; j < M; j++){
        dofJ = dofM->getGlobalId(dofField[j]);

        // If not a fixed Dof: count!
        if(dofJ != DofManager<scalar>::isFixedId())
          count++;
      }
    }
  }

  return count;
}

template<typename scalar>
void SystemAbstract<scalar>::
assemble(SolverMatrix<scalar>& A,
         SolverVector<scalar>& b,
         size_t elementId,
         const std::vector<Dof>& dofField,
         const std::vector<Dof>& dofTest,
         const FormulationBlock<scalar>& formulation){

  const size_t N = dofTest.size();
  const size_t M = dofField.size();

  size_t dofI;
  size_t dofJ;

  for(size_t i = 0; i < N; i++){
    dofI = dofM->getGlobalId(dofTest[i]);

    // If not a fixed Dof line: assemble
    if(dofI != DofManager<scalar>::isFixedId()){
      for(size_t j = 0; j < M; j++){
        dofJ = dofM->getGlobalId(dofField[j]);

        // If not a fixed Dof
        if(dofJ != DofManager<scalar>::isFixedId())
          A.add(dofI, dofJ, formulation.weak(i, j, elementId));

        // If fixed Dof (for column 'dofJ'):
        //    add to right hand side (with a minus sign) !
        else
          b.add(dofI,
                minusSign * dofM->getValue(dofField[j]) *
                            formulation.weak(i, j, elementId));
      }

      b.add(dofI, formulation.rhs(i, elementId));
    }
  }
}

template<typename scalar>
void SystemAbstract<scalar>::
assembleLHSOnly(SolverMatrix<scalar>& A,
                size_t elementId,
                const std::vector<Dof>& dofField,
                const std::vector<Dof>& dofTest,
                const FormulationBlock<scalar>& formulation){

  const size_t N = dofTest.size();
  const size_t M = dofField.size();

  size_t dofI;
  size_t dofJ;

  for(size_t i = 0; i < N; i++){
    dofI = dofM->getGlobalId(dofTest[i]);

    // If not a fixed Dof line: assemble
    if(dofI != DofManager<scalar>::isFixedId()){
      for(size_t j = 0; j < M; j++){
        dofJ = dofM->getGlobalId(dofField[j]);

        // If not a fixed Dof
        if(dofJ != DofManager<scalar>::isFixedId())
          A.add(dofI, dofJ, formulation.weak(i, j, elementId));
      }
    }
  }
}

template<typename scalar>
void SystemAbstract<scalar>::
assembleRHSOnly(SolverVector<scalar>& b,
                size_t elementId,
                const std::vector<Dof>& dofField,
                const std::vector<Dof>& dofTest,
                const FormulationBlock<scalar>& formulation){

  const size_t N = dofTest.size();
  const size_t M = dofField.size();

  size_t dofI;
  size_t dofJ;

  for(size_t i = 0; i < N; i++){
    dofI = dofM->getGlobalId(dofTest[i]);

    // If not a fixed Dof line: assemble rhs
    if(dofI != DofManager<scalar>::isFixedId()){
      for(size_t j = 0; j < M; j++){
        dofJ = dofM->getGlobalId(dofField[j]);

        // If a fixed Dof column: assemble rhs
        if(dofJ == DofManager<scalar>::isFixedId())
          b.add(dofI,
                minusSign * dofM->getValue(dofField[j]) *
                            formulation.weak(i, j, elementId));
      }

      b.add(dofI, formulation.rhs(i, elementId));
    }
  }
}

template<typename scalar>
void SystemAbstract<scalar>::
getProcSize(size_t nRow, size_t nProc, std::vector<size_t>& size){
  size_t i;
  size_t up  = ceil((double)(nRow) / (double)(nProc));
  size_t rem = nRow - up * nProc;

  size.resize(nProc);

  for(i = 0; i < nProc; i++)
    size[i] = up;

  for(i = 0; rem + size[i] <= 0; i++){
    rem    += size[i];
    size[i] = 0;
  }

  size[i] += rem;
}

template<typename scalar>
void SystemAbstract<scalar>::
getOwnership(const std::vector<size_t>& size, std::vector<size_t>& own){
  // Get full size //
  size_t nProc    = size.size();
  size_t fullSize = 0;

  for(size_t i = 0; i < nProc; i++)
    fullSize += size[i];

  // Alloc //
  own.resize(fullSize);

  // Populate //
  size_t offset = 0;
  size_t owner  = 0;
  for(size_t i = 0; i < nProc; i++){
    for(size_t j = 0; j < size[i]; j++)
      own[j + offset] = owner;

    offset += size[i];
    owner++;
  }
}

template<typename scalar>
void SystemAbstract<scalar>::
getProcMinRange(const std::vector<size_t>& size, std::vector<size_t>& min){
  const size_t nProc = size.size();
  min.resize(nProc);
  min[0] = 0;

  for(size_t i = 1, j = 0; i < nProc; i++, j++)
    min[i] = min[j] + size[j];
}

template<typename scalar>
void SystemAbstract<scalar>::
getProcMaxRange(const std::vector<size_t>& size, std::vector<size_t>& max){
  const size_t nProc = size.size();
  max.resize(nProc);
  max[0] = size[0];

  for(size_t i = 1, j = 0; i < nProc; i++, j++)
    max[i] = max[j] + size[i];
}

template<typename scalar>
void SystemAbstract<scalar>::
petscSparsity(int* nonZero,
              int* row, int* col, size_t size,
              std::vector<size_t>& minRange, std::vector<size_t>& maxRange,
              std::vector<size_t>& owner,
              bool isDiagonal){
  // Init
  size_t iLastRow = row[0] - 1; // row is assumed in Fortran style
  size_t iLastCol = col[0] - 1; // row is assumed in Fortran style
  size_t      own = owner[iLastRow];

  if(iLastRow >= minRange[own] && iLastRow < maxRange[own])
    if((iLastCol >= minRange[own] && iLastCol < maxRange[own]) == isDiagonal)
      nonZero[iLastRow]++;

  // Loop
  size_t iRow;
  size_t iCol;

  for(size_t i = 1; i < size; i++){
    iRow = row[i] - 1;  // row is assumed in Fortran style
    iCol = col[i] - 1;  // col is assumed in Fortran style
    own  = owner[iRow];

    if(iRow != iLastRow || iCol != iLastCol){
      // A new entry is found (row[] and col[] are asumed sorted)
      iLastRow = iRow;
      iLastCol = iCol;

      if(iRow >= minRange[own] && iRow < maxRange[own])
        if((iCol >= minRange[own] && iCol < maxRange[own]) == isDiagonal)
          nonZero[iRow]++;
    }
  }
}

template<typename scalar>
void SystemAbstract<scalar>::
petscSerialize(int* row, int* col, scalar* value, size_t size, Mat& A){
  // Loop
  int iRow;
  int iCol;

  for(size_t i = 0; i < size; i++){
    // Get Row & Col
    iRow = row[i] - 1; // row is assumed in Fortran style
    iCol = col[i] - 1; // col is assumed in Fortran style

    // Add
    MatSetValues(A, 1, &iRow, 1, &iCol, &value[i], ADD_VALUES);
  }
}
