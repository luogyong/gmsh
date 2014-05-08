////////////////////////////////////////////////
// Templates Implementations for SolverMUMPS: //
// Inclusion compilation model                //
//                                            //
// Damn you gcc: we want 'export' !           //
////////////////////////////////////////////////

template<typename scalar>
void SolverMUMPS<scalar>::copy(std::vector<Complex>& data,
                               mumps_double_complex** out){
  // Size
  const size_t size = data.size();

  // Alloc mumps struct
  *out = new mumps_double_complex[size];

  // Copy
  for(size_t i = 0; i < size; i++){
    (*out)[i].r = data[i].real();
    (*out)[i].i = data[i].imag();
  }
}

template<typename scalar>
void SolverMUMPS<scalar>::copy(SolverVector<Complex>& data,
                               mumps_double_complex** out){
  // Size
  const size_t size = data.getSize();

  // Alloc mumps struct
  *out = new mumps_double_complex[size];

  // Copy
  for(size_t i = 0; i < size; i++){
    (*out)[i].r = data.get(i).real();
    (*out)[i].i = data.get(i).imag();
  }
}

template<typename scalar>
void SolverMUMPS<scalar>::copy(SolverVector<double>& data,
                               std::vector<double>& out){
  // Size
  const size_t size = data.getSize();

  // Alloc
  out.resize(size);

  // Copy
  for(size_t i = 0; i < size; i++)
    out[i] = data.get(i);
}

template<typename scalar>
void SolverMUMPS<scalar>::copy(std::vector<double>& data,
                               fullVector<double>& out){
  // Size
  const size_t size = data.size();

  // Alloc
  out.resize(size);

  // Copy
  for(size_t i = 0; i < size; i++)
    out(i) = data[i];
}

template<typename scalar>
void SolverMUMPS<scalar>::copy(mumps_double_complex* data,
                               fullVector<Complex>& out,
                               size_t size){
  // Alloc
  out.resize(size);

  // Copy
  for(size_t i = 0; i < size; i++)
    out(i) = std::complex<double>(data[i].r, data[i].i);
}
