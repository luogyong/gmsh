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
SolverMUMPS<double>::SolverMUMPS(void){
  // State //
  hasMatrix    = false;
  hasRHS       = false;
  isFactorized = false;

  // Clear //
  row    = NULL;
  col    = NULL;
  value  = NULL;
  valueC = NULL;
  rhsR   = NULL;
  rhsC   = NULL;

  // MPI Self //
  const int FMPICommSelf = MPI_Comm_c2f(MPI_COMM_SELF);

  // Init MUMPS //
  idC = NULL;
  idR = new DMUMPS_STRUC_C;

  idR->job          =           -1; // Initialize MUMPS instance
  idR->par          =            1; // Host processor participates to the job
  idR->sym          =            0; // Unsymmetric matrix
  idR->comm_fortran = FMPICommSelf; // Use MPI COMM SELF (Fortran)

  dmumps_c(idR);                    // Do what is told in struct 'id'
}

template<>
void SolverMUMPS<double>::freeMatrix(void){
}

template<>
void SolverMUMPS<double>::freeRHS(void){
}

template<>
SolverMUMPS<double>::~SolverMUMPS(void){
  // Terminate MUMPS instance //
  idR->job = -2; // Destroy MUMPS Instance
  dmumps_c(idR);

  // Free Matrix and RHS //
  freeMatrix();
  freeRHS();

  // Free MUMPS pointer //
  delete idR;
}

template<>
void SolverMUMPS<double>::setMatrix(SolverMatrix<double>& A){
  // Size //
  nUnknown = A.nRows();

  // Is the given matrix square ? //
  if((size_t)(nUnknown) != A.nColumns())
    throw Exception("SolverMUMPS -- The given matrix is not square: (%d, %d)",
                    nUnknown, A.nColumns());

  // Get matrix data //
  const int nNZ = A.get(&row, &col, &value);

  // Define the matrix in MUMPS //
  idR->icntl[4]  = 0;           // Matrix in assembled format
  idR->icntl[17] = 0;           // Matrix is centralized on the host

  idR->n   = nUnknown;          // Size of the (square) matrix of unknown
  idR->nz  = nNZ;               // Number of non zero entries in the matrix
  idR->irn = row;               // Row vector
  idR->jcn = col;               // Column vector
  idR->a   = value;             // Value vector

  // State //
  hasMatrix    = true;
  isFactorized = false;
}

template<>
void SolverMUMPS<double>::setRHS(SolverVector<double>& rhs){
  // Check if coherent with matrix size //
  if(!hasMatrix)
    throw Exception("SolverMUMPS -- Cannot set RHS: need to set Matrix first");

  // Chech if size is coherent //
  if(idR->n != (int)(rhs.getSize()))
    throw Exception("%s -- %s: %d (matrix is %d)",
                    "SolverMUMPS", "The given RHS does not have the right size",
                    rhs.getSize(), idR->n);

  // Get vector data //
  rhsR = rhs.getData();

  // Define the right hand side in MUMPS //
  idR->rhs = rhsR; // Right hand side

  // State //
  hasRHS = true;
}

template<>
void SolverMUMPS<double>::factorize(void){
  // Check State //
  if(!hasMatrix)
    throw Exception("SolverMUMPS -- cannot factorize: no matrix set");

  if(isFactorized)
    return;

  // Output Settings //
  idR->icntl[0] = -1;  // No Output
  idR->icntl[1] = -1;  // ---------
  idR->icntl[2] = -1;  // ---------
  idR->icntl[3] = -1;  // ---------

  // Call MUMPS Analysis //
  idR->job = 1;
  dmumps_c(idR);

  // Call MUMPS Factorization //
  idR->job = 2;
  dmumps_c(idR);

  // State //
  isFactorized = true;
}

template<>
void SolverMUMPS<double>::solve(fullVector<double>& x){
  // Check State //
  if(!hasMatrix)
    throw Exception("SolverMUMPS -- cannot solve: needs matrix");

  if(!hasRHS)
    throw Exception("SolverMUMPS -- cannot solve: needs right hand side");

  if(!isFactorized)
    throw Exception("SolverMUMPS -- cannot solve: needs factorized matrix");

  // Solve with MUMPS //
  idR->job = 3;
  dmumps_c(idR);

  // The Right hand side is now the solution: copy it into 'x' //
  copy(rhsR, x, nUnknown);
}

template<>
void SolverMUMPS<double>::solve(SolverMatrix<double>& A,
                                SolverVector<double>& rhs,
                                fullVector<double>& x){
  setMatrix(A);
  setRHS(rhs);
  factorize();
  solve(x);
}
