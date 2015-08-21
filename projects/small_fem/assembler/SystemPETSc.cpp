#include "SmallFem.h"
#include "SystemPETSc.h"

template<>
PetscScalar* SystemPETSc<double>::toPetscScalar(double* a, size_t size){
  PetscScalar* petsc = new PetscScalar[size];
  for(size_t i = 0; i < size; i++)
    petsc[i] = Complex(a[i], 0.0);

  return petsc;
}

template<>
PetscScalar* SystemPETSc<Complex>::toPetscScalar(Complex* a, size_t size){
  PetscScalar* petsc = new PetscScalar[size];
  for(size_t i = 0; i < size; i++)
    petsc[i] = a[i];

  return petsc;
}

template<>
void SystemPETSc<double>::getSolution(void){
  const size_t size = this->dofM->getLocalSize();
  PetscScalar* solution;

  VecGetArray(xPetsc, &solution);
  x = new fullVector<double>(size);

  for(size_t i = 0; i < size; i++)
    (*x)(i) = solution[i].real();
}

template<>
void SystemPETSc<Complex>::getSolution(void){
  const size_t size = this->dofM->getLocalSize();
  PetscScalar* solution;

  VecGetArray(xPetsc, &solution);
  x = new fullVector<Complex>(size);

  for(size_t i = 0; i < size; i++)
    (*x)(i) = solution[i];
}
