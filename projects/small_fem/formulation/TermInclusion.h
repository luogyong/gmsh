/////////////////////////////////////////
// Templates Implementations for Term: //
// Inclusion compilation model         //
//                                     //
// Damn you gcc: we want 'export' !    //
/////////////////////////////////////////

template<typename scalar>
Term<scalar>::Term(void){
  // One cache per thread
  const size_t nThread = omp_get_max_threads();

  once    = new bool[nThread];
  lastId  = new size_t[nThread];
  lastI   = new size_t[nThread];
  lastCtr = new size_t[nThread];

  // We did not get into the cache
  for(size_t i = 0; i < nThread; i++)
    once[i] = false;
}

template<typename scalar>
Term<scalar>::~Term(void){
  delete[] once;
  delete[] lastId;
  delete[] lastI;
  delete[] lastCtr;

  for(size_t s = 0; s < nOrientation; s++)
    delete aM[s];

  delete[] aM;
}

template<typename scalar>
scalar Term<scalar>::getTermOutCache(size_t dofI, size_t dofJ, size_t elementId,
                                     size_t threadId) const{
  size_t i   = 0;
  size_t ctr = elementId;
  size_t off = (*orientationStat)[0];

  for(; elementId >= off && i < nOrientation; i++){
    off += (*orientationStat)[i + 1];
    ctr -= (*orientationStat)[i];
  }

  once[threadId]    = true;
  lastId[threadId]  = elementId;
  lastI[threadId]   = i;
  lastCtr[threadId] = ctr;

  return (*aM[i])(ctr, dofI * nFunction + dofJ);
}

template<typename scalar>
void Term<scalar>::allocA(size_t nFunction){
  // Alloc //
  aM = new fullMatrix<scalar>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    aM[s] = new fullMatrix<scalar>((*orientationStat)[s], nFunction);
}

template<typename scalar>
void Term<scalar>::computeA(fullMatrix<scalar>**& bM, fullMatrix<scalar>**& cM){
  // Fill //
  for(size_t s = 0; s < nOrientation; s++)
    // GEMM doesn't like matrices with 0 Elements
    if((*orientationStat)[s])
      aM[s]->gemm(*bM[s], *cM[s]);
}

template<typename scalar>
void Term<scalar>::clean(fullMatrix<scalar>**& bM, fullMatrix<scalar>**& cM){

  for(size_t s = 0; s < nOrientation; s++)
    delete cM[s];

  delete[] cM;

  for(size_t s = 0; s < nOrientation; s++)
    delete bM[s];

  delete[] bM;
}
