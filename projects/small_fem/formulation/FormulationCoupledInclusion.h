///////////////////////////////////////////////////////
// Templates Implementations for FormulationCoupled: //
// Inclusion compilation model                       //
//                                                   //
// Damn you gcc: we want 'export' !                  //
///////////////////////////////////////////////////////

template<typename scalar>
FormulationCoupled<scalar>::~FormulationCoupled(void){
}

template<typename scalar>
bool FormulationCoupled<scalar>::isBlock(void) const{
  return false;
}
