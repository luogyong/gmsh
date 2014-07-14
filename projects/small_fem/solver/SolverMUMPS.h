#ifndef _SOLVERMUMPS_H_
#define _SOLVERMUMPS_H_

#include "mumps_c_types.h"
#include "dmumps_c.h"
#include "zmumps_c.h"

#include "Solver.h"
#include "SmallFem.h"

/**
   @class SolverMUMPS
   @brief A Solver using the MUMPS library

   This class implements a Solver using the
   MUltifrontal Massively Parallel sparse direct %Solver (MUMPS) library.

   This library can be download at
    <a href="http://mumps.enseeiht.fr">http://mumps.enseeiht.fr</a> or
    <a href="http://graal.ens-lyon.fr/MUMPS">http://graal.ens-lyon.fr/MUMPS</a>.
*/

template<typename scalar>
class SolverMUMPS: public Solver<scalar>{
 private:
  // State //
  bool hasMatrix;
  bool hasRHS;
  bool isFactorized;

  // Size //
  int nUnknown;

  // Matrix //
  int*                  row;
  int*                  col;
  scalar*               value;  // Complex AND Real case
  mumps_double_complex* valueC; // Complex case only

  // RHS //
  double*               rhsR; // Real    case
  mumps_double_complex* rhsC; // Complex case

  // MUMPS //
  DMUMPS_STRUC_C* idR;
  ZMUMPS_STRUC_C* idC;

 public:
  SolverMUMPS(void);

  virtual ~SolverMUMPS(void);

  virtual void solve(SolverMatrix<scalar>& A,
                     SolverVector<scalar>& rhs,
                     fullVector<scalar>& x);

  void setMatrix(SolverMatrix<scalar>& A);
  void setRHS(SolverVector<scalar>& rhs);

  void factorize(void);
  void solve(fullVector<scalar>& x);

 private:
  void freeMatrix(void);
  void freeRHS(void);

  static void copy(Complex*              data, mumps_double_complex** out,
                   size_t size);
  static void copy(double*               data, fullVector<double>&    out,
                   size_t size);
  static void copy(mumps_double_complex* data, fullVector<Complex>&   out,
                   size_t size);
};

/**
   @fn SolverMUMPS::SolverMUMPS
   Instanciates a new SolverMUMPS
   **

   @fn SolverMUMPS::~SolverMUMPS
   Deletes this SolverMUMPS
   **

   @fn SolverMUMPS::setMatrix
   @param A A SolverMatrix
   The given matrix is now this Solver matrix
   **

   @fn SolverMUMPS::setRHS
   @param rhs A SolverVector
   The given vector is now this Solver right hand side
   **

   @fn SolverMUMPS::factorize
   Computes the LU decomposition of the given matrix
   **

   @fn SolverMUMPS::solve(fullVector<scalar>& x)
   @param x The fullVector where the solution will be stored
   Uses the previous SolverMUMPS::factorize() and SolverMUMPS::setRHS()
   to compute the solution of the linear system
*/


//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SolverMUMPSInclusion.h"

#endif
