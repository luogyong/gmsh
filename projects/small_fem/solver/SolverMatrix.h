#ifndef _SOLVERMATRIX_H_
#define _SOLVERMATRIX_H_

#include <string>
#include <vector>
#include <omp.h>

/**
   @class SolverMatrix
   @brief A class handling a Solver linear system matrix

   This class represents a matrix used by a Solver.

   A SolverMatrix data can be recovered into
   three vectors, r[], c[] and a[], such that a[k] is
   the entry of the matrix at row (r[k] - 1) and column (c[k] - 1).
   Values at the same position should be added.

   During the construction of the matrix, multiple values
   can be added in a thread-safe manner.
*/

template<typename scalar>
class SolverMatrix{
 private:
  static const size_t maxSizeT;

 private:
  // Size //
  size_t nRow;
  size_t nCol;

  // Non zero terms & threads offset //
  std::vector<size_t> nxt;
  std::vector<size_t> min;
  std::vector<size_t> max;

  // Data //
  int*    row;
  int*    col;
  scalar* value;

 public:
   SolverMatrix(size_t nRow, size_t nCol, std::vector<size_t> nonZero);
  ~SolverMatrix(void);

  size_t nRows(void) const;
  size_t nColumns(void) const;

  void   add(size_t row, size_t col, scalar   value);
  size_t get(int**  row, int**  col, scalar** value);
  size_t get(int**  row, int**  col, scalar** value, bool sorted);

  std::string toString(void) const;
  std::string toMatlab(std::string matrixName) const;
  void writeToMatlabFile(std::string fileName, std::string matrixName) const;

 private:
  SolverMatrix(void);

  void   sort(void);
  void   quickSort(size_t start, size_t end);
  size_t partition(size_t start, size_t end);
  void   heapify(size_t size, size_t node);
  void   swap(size_t i, size_t j);

  // void        sortAndReduce(void) const;
  // std::string matlabCommon(std::string matrixName) const;
};

/**
   @fn SolverMatrix::SolverMatrix
   @param nRow The number of row of this SolverMatrix
   @param nCol The number of column of this SolverMatrix
   @param nonZero Per-Thread non zero term number

   Instanciates an new SolverMatrix of zero values with the given number of
   rows and columns.

   The number of non-zero term per-thread should also be given.
   nonZero[i] is the number of non-zero term that will be added by thread i.
   **

   @fn SolverMatrix::~SolverMatrix

   Deletes this SolverMatrix
   **

   @fn SolverMatrix::nRows;
   @return Returns the number of rows of this SolverMatrix
   **

   @fn SolverMatrix::nColumns;
   @return Returns the number of columns of this SolverMatrix
   **

   @fn SolverMatrix::add
   @param row A row index of this SolverMatrix
   @param col A column index of this SolverMatrix
   @param value A real number

   Adds the given value at the given row and column.
   Indices are given in a C style manner (that is starting at 0)
   **

   @fn SolverMatrix::get(int**  row, int**  col, scalar** value)
   @param row A pointer to an int*
   @param col A pointer to an int*
   @param value A pointer to a scalar*

   @return Returns the size of this SolverMatrix non zero entries

   Same as SolverMatrix::get(row, col, value, false).
   **

   @fn SolverMatrix::get(int**  row, int**  col, scalar** value, bool sorted)
   @param row A pointer to an int*
   @param col A pointer to an int*
   @param value A pointer to a scalar*
   @param sorted A boolean value

   @return Returns the size of this SolverMatrix non zero entries

   *row now points the memory array storing the rows.
   *col now points the memory array storing the columns.
   *value now points the memory array storing the values.

   The three memory arrays are such that:
   A[(*row)[k] - 1, (*col)[k] - 1] = *(value)[k],
   where A is this SolverMatrix.

   If sorted is true, the returned vectors will be sorted
   (first by row, and then by columns).
   If sorted is false, the returned vectors will be unsorted.
   **

   @fn SolverMatrix::toString
   @return Returns a string describing this SolverMatrix
   **

   @fn SolverMatrix::toMatlab(std::string matrixName) const
   @param matrixName A string
   @return Returns a string that can be used in Octave/Matlab
   to reproduce this SolverMatrix, whose name will be the given one
   **

   @fn SolverMatrix::writeToMatlabFile
   @param fileName A string
   @param matrixName A string

   Writes this matrix in Octave/Matlab format into the given file,
   and with the given name
 */


//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SolverMatrixInclusion.h"

//////////////////////
// Inline Functions //
//////////////////////

template<typename scalar>
inline void SolverMatrix<scalar>::add(size_t row, size_t col, scalar value){
  // Get Thread num //
  const int threadId = omp_get_thread_num();

  // Next index //
  const size_t myNxt = this->nxt[threadId];

  // Add term (convert to Fortran style) //
  this->row[myNxt]   = row + 1;
  this->col[myNxt]   = col + 1;
  this->value[myNxt] = value;

  // Done //
  this->nxt[threadId]++;
}

template<typename scalar>
void SolverMatrix<scalar>::swap(size_t i, size_t j){
  int tmpRow   = this->row[i];
  this->row[i] = this->row[j];
  this->row[j] = tmpRow;

  int tmpCol   = this->col[i];
  this->col[i] = this->col[j];
  this->col[j] = tmpCol;

  scalar tmpValue = this->value[i];
  this->value[i]  = this->value[j];
  this->value[j]  = tmpValue;
}

#endif
