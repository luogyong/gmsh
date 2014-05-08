#include "SmallFem.h"

#include "SolverMUMPS.h"
#include "Exception.h"
#include "mpi.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////
// WARNING                                                               //
// -------                                                               //
//                                                                       //
// Implementations are in a weird order...                               //
// Mandatory to avoid : specialization of 'template' after instantiation //
//                                                                       //
// I definitely HATE gcc for its template stuffs...                      //
///////////////////////////////////////////////////////////////////////////

template<>
SolverMUMPS<Complex>::SolverMUMPS(void){
  // State //
  hasMatrix    = false;
  hasRHS       = false;
  isFactorized = false;

  // Clear //
  row.clear();
  col.clear();
  value.clear();
  valueC = NULL;
  rhsR.clear();
  rhsC   = NULL;

  // MPI Self //
  const int FMPICommSelf = MPI_Comm_c2f(MPI_COMM_SELF);

  // Init MUMPS //
  idR = NULL;
  idC = new ZMUMPS_STRUC_C;

  idC->job          =           -1; // Initialize MUMPS instance
  idC->par          =            1; // Host processor participates to the job
  idC->sym          =            0; // Unsymmetric matrix
  idC->comm_fortran = FMPICommSelf; // Use MPI COMM SELF (Fortran)

  zmumps_c(idC);                    // Do what is told in struct 'id'
}

template<>
void SolverMUMPS<Complex>::freeMatrix(void){
  if(valueC)
    delete[] valueC;
}

template<>
void SolverMUMPS<Complex>::freeRHS(void){
  if(rhsC)
    delete[] rhsC;
}

template<>
SolverMUMPS<Complex>::~SolverMUMPS(void){
  // Terminate MUMPS instance //
  idC->job = -2; // Destroy MUMPS Instance
  zmumps_c(idC);

  // Free Matrix and RHS //
  freeMatrix();
  freeRHS();

  // Free MUMPS pointer //
  delete idC;
}

template<>
void SolverMUMPS<Complex>::setMatrix(SolverMatrix<Complex>& A){
  // Size //
  const int size = A.nRows();

  // Is the given matrix square ? //
  if((size_t)(size) != A.nColumns())
    throw Exception("SolverMUMPS -- The given matrix is not square: (%d, %d)",
                    size, A.nColumns());

  // Serialize the given matrix //
  const int nNZ = A.serialize(row, col, value);

  // Convert into MUMPS Complex data //
  copy(value, &valueC);

  // Clear value (not needed any more) //
  value.clear();

  // Define the matrix in MUMPS //
  idC->icntl[4]  = 0;        // Matrix in assembled format
  idC->icntl[17] = 0;        // Matrix is centralized on the host

  idC->n   = size;           // Size of the (square) matrix of unknown
  idC->nz  = nNZ;            // Number of non zero entries in the matrix
  idC->irn = row.data();     // Row vector
  idC->jcn = col.data();     // Column vector
  idC->a   = valueC;         // Value vector

  // State //
  hasMatrix    = true;
  isFactorized = false;
}

template<>
void SolverMUMPS<Complex>::setRHS(SolverVector<Complex>& rhs){
  // Check if coherent with matrix size //
  if(!hasMatrix)
    throw Exception("SolverMUMPS -- Cannot set RHS: need to set Matrix first");

  // Chech if size is coherent //
  if(idC->n != (int)(rhs.getSize()))
    throw Exception("%s -- %s: %d (matrix is %d)",
                    "SolverMUMPS", "The given RHS does not have the right size",
                    rhs.getSize(), idC->n);

  // Copy into MUMPS Complex data //
  copy(rhs, &rhsC);

  // Define the right hand side in MUMPS //
  idC->rhs = rhsC;        // Right hand side

  // State //
  hasRHS = true;
}

template<>
void SolverMUMPS<Complex>::factorize(void){
  // Check State //
  if(!hasMatrix)
    throw Exception("SolverMUMPS -- cannot factorize: no matrix set");

  if(isFactorized)
    return;

  // Output Settings //
  idC->icntl[0] = -1;  // No Output
  idC->icntl[1] = -1;  // ---------
  idC->icntl[2] = -1;  // ---------
  idC->icntl[3] = -1;  // ---------

  // Call MUMPS Analysis //
  idC->job = 1;
  zmumps_c(idC);

  // Call MUMPS Factorization //
  idC->job = 2;
  zmumps_c(idC);

  // State //
  isFactorized = true;
}

template<>
void SolverMUMPS<Complex>::solve(fullVector<Complex>& x){
  // Check State //
  if(!hasMatrix)
    throw Exception("SolverMUMPS -- cannot solve: needs matrix");

  if(!hasRHS)
    throw Exception("SolverMUMPS -- cannot solve: needs right hand side");

  if(!isFactorized)
    throw Exception("SolverMUMPS -- cannot solve: needs factorized matrix");

  // Solve with MUMPS //
  idC->job = 3;
  zmumps_c(idC);

  // The Right hand side is now the solution: copy it into 'x' //
  copy(rhsC, x, idC->n);
}

template<>
void SolverMUMPS<Complex>::solve(SolverMatrix<Complex>& A,
                                 SolverVector<Complex>& rhs,
                                 fullVector<Complex>& x){
  setMatrix(A);
  setRHS(rhs);
  factorize();
  solve(x);
}
