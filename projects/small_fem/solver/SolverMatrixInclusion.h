/////////////////////////////////////////////////
// Templates Implementations for SolverMatrix: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include "Exception.h"

template<typename scalar>
const size_t SolverMatrix<scalar>::maxSizeT = 0 - 1;

template<typename scalar>
SolverMatrix<scalar>::SolverMatrix(void){
}

template<typename scalar>
SolverMatrix<scalar>::
SolverMatrix(size_t nRow, size_t nCol, std::vector<size_t> nonZero){
  // Matrix size //
  this->nRow = nRow;
  this->nCol = nCol;

  // Get num threads //
  size_t nThreads;
  #pragma omp parallel
  {
    #pragma omp master
    nThreads = omp_get_num_threads();
  }

  // Check thread numbers //
  if(nonZero.size() != nThreads)
    throw Exception("%s %s %s", "SolverMatrix: ", "number of non zero entry ",
                    "doesn't match number of threads");

  // Ranges //
  nxt.resize(nThreads);
  max.resize(nThreads);
  min.resize(nThreads);

  min[0] = 0;
  for(size_t i = 1, j = 0; i < nThreads; i++, j++)
    min[i] = min[j] + nonZero[j];

  max[0] = nonZero[0];
  for(size_t i = 1, j = 0; i < nThreads; i++, j++)
    max[i] = max[j] + nonZero[i];

  for(size_t i = 0; i < nThreads; i++)
    nxt[i] = min[i];

  // Matrix data //
  row   = new int[max[nThreads - 1]];
  col   = new int[max[nThreads - 1]];
  value = new scalar[max[nThreads - 1]];
}

template<typename scalar>
SolverMatrix<scalar>::~SolverMatrix(void){
  delete[] row;
  delete[] col;
  delete[] value;
}

template<typename scalar>
size_t SolverMatrix<scalar>::nRows(void) const{
  return nRow;
}

template<typename scalar>
size_t SolverMatrix<scalar>::nColumns(void) const{
  return nCol;
}

template<typename scalar>
size_t SolverMatrix<scalar>::get(int** row, int** col, scalar** value){
  return get(row, col, value, false);
}

template<typename scalar>
size_t SolverMatrix<scalar>::
get(int** row, int** col, scalar** value, bool sorted){
  if(sorted)
    sort();

  *row   = this->row;
  *col   = this->col;
  *value = this->value;

  return this->max[this->max.size() - 1];
}

template<typename scalar>
std::string SolverMatrix<scalar>::toString(void) const{
  std::stringstream stream;

  return stream.str();
}

template<typename scalar>
void SolverMatrix<scalar>::writeToMatlabFile(std::string fileName,
                                             std::string matrixName) const{
  std::ofstream stream;
  stream.open(fileName.c_str());
  stream << toMatlab(matrixName) << std::endl;
  stream.close();
}

template<typename scalar>
void SolverMatrix<scalar>::sort(void){
  // Array size //
  const size_t size = this->max[this->max.size() - 1];

  // Quick sort //
  //quickSort(0, size - 1);

  // Heap sort //
  // Build heap
  for(size_t i = size / 2; i != 0; i--)
    heapify(size, i);
  heapify(size, 0);

  // Sort
  for(size_t i = size - 1; i != 0; i--){
    swap(i, 0);
    heapify(i, 0);
  }
}

template<typename scalar>
void SolverMatrix<scalar>::quickSort(size_t start, size_t end){
  // If there's nothing to sort, return //
  if(start >= end)
    return;

  // Partition of row, col and value //
  size_t pivotPosition = partition(start, end);

  // Sort left (if it exists) //
  if(pivotPosition > 0)
    quickSort(start, pivotPosition - 1);

  // Sort right (if it exists) //
  if(pivotPosition < maxSizeT)
    quickSort(pivotPosition + 1, end);
}

template<typename scalar>
size_t SolverMatrix<scalar>::partition(size_t start, size_t end){
  // Get pivot at the end of row //
  int pivotRow = this->row[end];
  int pivotCol = this->col[end];
  size_t     i = start - 1;

  // Look for smaller elements than the pivotRow.            //
  // Keep smaller elements to the left,                      //
  // and bigger elements to the right of row, col and value. //

  // If equal element is found, do the same with pivotCol.   //
  for(size_t j = start; j < end; ++j){
    if(row[j] < pivotRow){
      ++i;
      swap(i, j); // Swap row[i & j], col[i & j] and value[i & j]
    }

    if(row[j] == pivotRow && col[j] < pivotCol){
      ++i;
      swap(i, j); // Swap row[i & j], col[i & j] and value[i & j]
    }
  }

  // Swap pivot with the leftmost bigger element //
  ++i;
  swap(i, end);   // Swap row[i & j], col[i & j] and value[i & j]

  return i;
}

template<typename scalar>
void SolverMatrix<scalar>::heapify(size_t size, size_t node){
  const size_t left  = (2 * node) + 1; // Left son
  const size_t right = left + 1;       // Right son

  size_t maxNode = node; // We think that the given node is the biggest

  // If left son exists,
  // and if maxNode is smaller than left son
  if(left < size && row[maxNode] < row[left])
    maxNode = left; // Now, left node is the maxNode

  // If left son is equal to maxNode: look to columns
  if(left < size && row[maxNode] == row[left] && col[maxNode] < col[left])
    maxNode = left; // Now, left node is the maxNode

  // If right son exists,
  // and if maxNode is smaller than right son
  if(right < size && row[maxNode] < row[right])
    maxNode = right; // now, right node is the maxNode

  // If right son is equal to maxNode: look to columns
  if(right < size && row[maxNode] == row[right] && col[maxNode] < col[right])
    maxNode = right; // Now, right node is the maxNode

  // If maxNode is not the given node
  //    --> swap and heapify from new position
  if(maxNode != node){
    swap(node, maxNode);
    heapify(size, maxNode);
  }
}

/*
template<typename scalar>
void SolverMatrix<scalar>::sortAndReduce(void) const{
  // Each thread takes a set of rows
  #pragma omp parallel for
  for(size_t i = 0; i < nRow; i++){
    // Sort row[i]
    data[i].sort(sortPredicate);

    // Remove duplicate in row[i] by adding them
    typename std::list<std::pair<size_t, scalar> >::iterator start =
      data[i].begin();
    typename std::list<std::pair<size_t, scalar> >::iterator stop  =
      data[i].begin();
    typename std::list<std::pair<size_t, scalar> >::iterator end   =
      data[i].end();

    size_t idx;
    scalar val;

    if(start != end){
      idx = start->first;
      val = start->second;
      stop++;
    } // Not needed but compiler cannot see that it is redundant with while loop

    while(start != end){
      // If start == stop, we need to resart
      if(start == stop){
        idx = start->first;
        val = start->second;
        stop++;
      }

      // If we get a element with the same column index:
      // we continue to accumate in 'val'
      if(stop != end && stop->first == idx){
        val += stop->second;
        stop++;
      }

      // Else, we earse the superflous element
      else{
        // Create an iterator to the previous element
        typename std::list<std::pair<size_t, scalar> >::iterator tmp = stop;
        tmp--;

        if(tmp != start){
          // We need to erase superflous elements
          // But first, set 'tmp' to value accumulator
          tmp->second = val;

          // Then erase from start to tmp (excluded)
          data[i].erase(start, tmp);
        }

        // Restart from 'stop'
        start = stop;
      }
    }
  }
}

template<typename scalar>
std::string SolverMatrix<scalar>::matlabCommon(std::string matrixName) const{
  // Stream
  typename std::list<std::pair<size_t, scalar> >::iterator it;
  typename std::list<std::pair<size_t, scalar> >::iterator end;
  std::stringstream stream;

  // Call to 'sparse'
  stream << matrixName << " = sparse(";

  // Rows
  stream << "[";
  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << i << ", ";
  }
  stream << "], ";

  // Columns
  stream << "[";
  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << it->first << ", ";
  }
  stream << "], ";

  // Return
  return stream.str();
}
*/
