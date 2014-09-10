#ifndef _MPIDOFMAP_H_
#define _MPIDOFMAP_H_

#include <map>
#include "Dof.h"

/**
   @class MPIDofMap
   @brief MPI distributed Dof map

   Helper methods for distributing Dof maps throug MPI
*/

template<typename T>
class MPIDofMap{
 public:
   MPIDofMap(void);
  ~MPIDofMap(void);

 public:
  static void getGlobalMap(const std::map<Dof, T>& local,
                           std::map<Dof, T>& global);

  static void getDofOwners(const std::map<Dof, T>& local,
                           std::multimap<Dof, int>& owners);

 private:
  static int  serialize(const std::map<Dof, T>& in, int** entity, int** type);
  static void unserialize(std::map<Dof, T>& map,
                          int* entity, int* type, int size);

  static void unserialize(std::multimap<Dof, int>& map,
                          int* entity, int* type, int* owner, int size);

  static int* gatherSize(int mySize, int* sum);
  static int* cpteStride(int* size);
  static int* exchange(int* myData,
                       int  mySize, int sizeSum, int* size, int* stride);
};

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "MPIDofMapInclusion.h"

#endif
