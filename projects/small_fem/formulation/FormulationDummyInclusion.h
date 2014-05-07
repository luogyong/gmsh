/////////////////////////////////////////////////////
// Templates Implementations for FormulationDummy: //
// Inclusion compilation model                     //
//                                                 //
// Damn you gcc: we want 'export' !                //
/////////////////////////////////////////////////////

#include "FormulationDummy.h"

template<typename scalar>
FormulationDummy<scalar>::FormulationDummy(void){
  // Empty List //
  fList.clear();
}

template<typename scalar>
FormulationDummy<scalar>::~FormulationDummy(void){
}

template<typename scalar>
const std::list<const FormulationBlock<scalar>*>&
FormulationDummy<scalar>::getFormulationBlocks(void) const{
  return fList;
}

template<typename scalar>
bool FormulationDummy<scalar>::isBlock(void) const{
  return false;
}
