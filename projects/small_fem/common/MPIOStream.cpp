#include "mpi.h"
#include "MPIOStream.h"

MPIOStream::MPIOStream(int rank, std::ostream& ostream){
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  this->myStream = (myRank == rank);
  this->ostream  = &ostream;
}

MPIOStream::~MPIOStream(void){
}
