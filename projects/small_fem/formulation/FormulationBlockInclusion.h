/////////////////////////////////////////////////////
// Templates Implementations for FormulationBlock: //
// Inclusion compilation model                     //
//                                                 //
// Damn you gcc: we want 'export' !                //
/////////////////////////////////////////////////////

template<typename scalar>
FormulationBlock<scalar>::~FormulationBlock(void){
}

template<typename scalar>
bool FormulationBlock<scalar>::isBlock(void) const{
  return true;
}
