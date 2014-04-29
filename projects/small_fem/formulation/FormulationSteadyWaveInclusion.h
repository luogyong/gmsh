//////////////////////////////////////////////////////////
// Templates Implementations for FormulationSteadyWave: //
// Inclusion compilation model                          //
//                                                      //
// Damn you gcc: we want 'export' !                     //
//////////////////////////////////////////////////////////

template<typename scalar>
FormulationSteadyWave<scalar>::
FormulationSteadyWave(const GroupOfElement& domain,
                      const FunctionSpace& fs,
                      scalar k){
  // Wave Squared //
  kSquare = k * k;

  // Stiffness and mass formulations //
  stiff = new FormulationStiffness<scalar>(domain, fs, fs);
  mass  = new FormulationMass<scalar>(domain, fs, fs);
}

template<typename scalar>
FormulationSteadyWave<scalar>::~FormulationSteadyWave(void){
  delete stiff;
  delete mass;
}

template<typename scalar>
scalar FormulationSteadyWave<scalar>::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return
    stiff->weak(dofI, dofJ, elementId) -
     mass->weak(dofI, dofJ, elementId) * kSquare;
}

template<typename scalar>
scalar FormulationSteadyWave<scalar>::
rhs(size_t equationI, size_t elementId) const{
  return 0;
}

template<typename scalar>
const FunctionSpace& FormulationSteadyWave<scalar>::field(void) const{
  return stiff->field();
}

template<typename scalar>
const FunctionSpace& FormulationSteadyWave<scalar>::test(void) const{
  return stiff->test();
}

template<typename scalar>
const GroupOfElement& FormulationSteadyWave<scalar>::domain(void) const{
  return stiff->domain();
}

template<typename scalar>
bool FormulationSteadyWave<scalar>::isBlock(void) const{
  return true;
}
