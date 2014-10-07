//////////////////////////////////////////////
// Templates Implementations for TermDummy: //
// Inclusion compilation model              //
//                                          //
// Damn you gcc: we want 'export' !         //
//////////////////////////////////////////////

template<typename scalar>
TermDummy<scalar>::TermDummy(void){
  dummyStat.resize(1);
  dummyStat[0] = 0;

  this->orientationStat = &dummyStat;
  this->nOrientation    = 0;
  this->nFunctionField  = 0;
  this->nFunctionTest   = 0;

  this->allocA(this->nFunctionField * this->nFunctionTest);
}

template<typename scalar>
TermDummy<scalar>::~TermDummy(void){
}
