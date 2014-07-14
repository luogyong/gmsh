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
  return dofM.getUnfixedDofNumber();
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
    dofM.fixValue(it->first, it->second);
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
    dofI = dofM.getGlobalId(dofTest[i]);

    // If not a fixed Dof line
    if(dofI != DofManager<scalar>::isFixedId()){
      for(size_t j = 0; j < M; j++){
        dofJ = dofM.getGlobalId(dofField[j]);

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
    dofI = dofM.getGlobalId(dofTest[i]);

    // If not a fixed Dof line: assemble
    if(dofI != DofManager<scalar>::isFixedId()){
      for(size_t j = 0; j < M; j++){
        dofJ = dofM.getGlobalId(dofField[j]);

        // If not a fixed Dof
        if(dofJ != DofManager<scalar>::isFixedId())
          A.add(dofI, dofJ, formulation.weak(i, j, elementId));

        // If fixed Dof (for column 'dofJ'):
        //    add to right hand side (with a minus sign) !
        else
          b.add(dofI,
                minusSign * dofM.getValue(dofField[j]) *
                           formulation.weak(i, j, elementId));
      }

      b.add(dofI, formulation.rhs(i, elementId));
    }
  }
}

template<typename scalar>
void SystemAbstract<scalar>::
assembleRHSOnly(SolverVector<scalar>& b,
                size_t elementId,
                const std::vector<Dof>& dofTest,
                const FormulationBlock<scalar>& formulation){

  const size_t N = dofTest.size();
  size_t dofI;

  for(size_t i = 0; i < N; i++){
    dofI = dofM.getGlobalId(dofTest[i]);

    // If not a fixed Dof line: assemble
    if(dofI != DofManager<scalar>::isFixedId())
      b.add(dofI, formulation.rhs(i, elementId));

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
getProcMinRange(const std::vector<size_t>& size,
                std::vector<size_t>& min){

  const size_t nProc = size.size();

  min.resize(nProc);
  min[0] = 0;

  for(size_t i = 1, j = 0; i < nProc; i++, j++)
    min[i] = min[j] + size[j];
}

template<typename scalar>
void SystemAbstract<scalar>::
getProcMaxRange(const std::vector<size_t>& size,
                std::vector<size_t>& max){

  const size_t nProc = size.size();
  max.resize(nProc);
  max[0] = size[0];

  for(size_t i = 1, j = 0; i < nProc; i++, j++)
    max[i] = max[j] + size[i];
}

template<typename scalar>
void SystemAbstract<scalar>::
petscSparsity(PetscInt* nonZero,
              int* row, int* col, size_t size,
              int iMin, int iMax, bool isDiagonal){
  // Loop
  for(size_t i = 0; i < size; i++){
    int iRow = row[i] - 1;   // row is assumed in Fortran style

    if(iRow >= iMin && iRow < iMax){
      int iCol = col[i] - 1; // col is assumed in Fortran style

      if((iCol >= iMin && iCol < iMax) == isDiagonal)
        nonZero[iRow - iMin]++;
    }
  }
}

template<typename scalar>
void SystemAbstract<scalar>::
petscSerialize(int rowMin, int rowMax,
               int* row, int* col, scalar* value, size_t size, Mat& A){
  // Loop
  int iRow;
  int iCol;

  for(size_t i = 0; i < size; i++){
    // Get Row
    iRow = row[i] - 1;   // row is assumed in Fortran style

    // Is in row range ?
    if(iRow >= rowMin && iRow < rowMax){
      iCol = col[i] - 1; // col is assumed in Fortran style
      MatSetValues(A, 1, &iRow, 1, &iCol, &value[i], ADD_VALUES);
    }
  }
}
