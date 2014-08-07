#ifndef _DOFMANAGER_H_
#define _DOFMANAGER_H_

#include <string>
#include <vector>
#include <map>
#include <set>

#include "Dof.h"

/**
   @class DofManager
   @brief This class manages the degrees of freedom (Dof)

   This class numbers the degrees of freedom (Dof).

   It can map a Dof to a unique number, called global ID.

   In addtion, this class allows to assign a Dof to a given value.
   A Dof that has been assigned to a value is called a fixed Dof.
   The global ID of a fixed Dof is not unique and is equal to
   DofManager::isFixedId().

   A DofManager can be distributed across multiple computing node,
   or it can be local to each computing node.
   It the first case, every Dof will have a unique number, no matter the node.
   In the second case, Dof may have different unique number,
   depending on the node.

   Finaly, the global IDs given to the unfixed Dof%s ranges from 0 to
   (total number of Dof - number of fixed Dof - 1).
*/

template<typename scalar>
class DofManager{
 private:
  static const size_t isFixed;
  static const size_t isUndef;

  bool local;

  std::vector<std::vector<size_t> > globalIdV;
  std::map<Dof, size_t>             globalIdM;
  std::map<Dof, scalar>             fixedDof;

  size_t nTotUnfixedLocalDof;
  size_t nTotUnfixedGlobalDof;

  size_t first;
  size_t last;

 public:
  static const size_t isFixedId(void);

 public:
   DofManager(bool isLocal);
   DofManager(void);
  ~DofManager(void);

  bool   isLocal(void)       const;
  size_t getLocalSize(void)  const;
  size_t getGlobalSize(void) const;

  void addToDofManager(const std::set<Dof>& dof);
  void addToDofManager(const std::vector<std::vector<Dof> >& dof);
  void generateGlobalIdSpace(void);

  size_t getGlobalId(const Dof& dof)     const;
  size_t getGlobalIdSafe(const Dof& dof) const;

  bool   fixValue(const Dof& dof, scalar value);
  scalar getValue(const Dof& dof) const;

  std::string toString(void) const;

 private:
  void localSpace (void);
  void globalSpace(void);
  void vectorize(void);

  std::pair<bool, size_t> findSafe(const Dof& dof) const;

  std::string toStringFromMap(void) const;
  std::string toStringFromVec(void) const;

 private:
  template<typename T>
  static void getGlobalMap(std::map<Dof, T>& local, std::map<Dof, T>& global);

  template<typename T>
  static int serialize(const std::map<Dof, T>& in, int** entity, int** type);

  template<typename T>
  static void unserialize(std::map<Dof, T>& map,
                          int* entity, int* type, int size);

  static void count(std::map<Dof, size_t>& dof);
  static void retag(std::map<Dof, size_t>& dof, std::map<Dof, scalar>& fix);
  static int* gatherSize(int mySize, int* sum);
  static int* cpteStride(int* size);
  static int* exchange(int* myData,
                       int  mySize, int sizeSum, int* size, int* stride);
};


/**
   @fn DofManager::isFixedId

   Fixed Dof got a special global ID (which is the highest possible number).

   @return Returns the special ID of a fixed Dof
   **


   @fn DofManager::DofManager(bool)
   @param isLocal A boolean value

   Instantiates a new DofManager.
   If isLocal is true, a local DofManager will be created.
   If isLocal is false, a distributed DofManager will be created.
   **

   @fn DofManager::DofManager(void)

   Same as DofManager::DofManager(true)
   **

   @fn DofManager::~DofManager

   Deletes this DofManager
   **

   @fn DofManager::isLocal(void)
   @return Returns true if this DofManager is local,
   and false if it is distributed.
   **

   @fn DofManager::getLocalSize
   @return Returns the number of local Unfixed Dof%s in this DofManager
   **

   @fn DofManager::getGlobalSize
   @return Returns the number of Unfixed Dof%s in this DofManager
   across all nodes

   If DofManager::isLocal() is true, then DofManager::getGlobalSize()
   is equal to DofManager::getLocalSize()
   **

   @fn DofManager::addToDofManager(const std::set<Dof>& dof);
   @param dof A set of Dof%s

   Adds the given Dof%s in this DofManager.
   The same Dof may be insterd multiple time,
   but it will be given the same unique ID.
   **

   @fn DofManager::addToDofManager(const std::vector<std::vector<Dof> >& dof);
   @param dof A vector of vector of Dof%s

   Adds the given Dof%s in this DofManager.
   The same Dof may be insterd multiple time,
   but it will be given the same unique ID.
   **

   @fn DofManager::generateGlobalIdSpace

   Numbers every non fixed Dof of this DofManager.
   Each Dof will be given a unique global ID.
   **

   @fn DofManager::getGlobalId
   @param dof The Dof from which we want the global ID
   @return Returns the global ID of the given Dof

   If the requested Dof has not been added to this DofManager,
   or if DofManager::generateGlobalIdSpace() has not been called,
   the behaviour of this method is unpredicable.
   Actually, it will most probably lead to a process crash.

   @see DofManager::getGlobalIdSafe
   **

   @fn DofManager::getGlobalIdSafe
   @param dof The Dof from which we want the global ID
   @return Returns the global ID of the given Dof

   If the requested Dof has not been added to this DofManager,
   or if DofManager::generateGlobalIdSpace() has not been called,
   an Exception is thrown.

   This method is safer but slower than DofManager::getGlobalId().

   @see DofManager::getGlobalId
   **

   @fn DofManager::fixValue
   @param dof A Dof
   @param value A real number

   Fixes the given Dof to the given value.

   @return Returns:
   @li true, if the operation is a success
   @li false, otherwise

   Here are three important cases, where DofManager::fixValue() will fail:
   @li The given Dof is not in this DofManager
   @li The given Dof is already fixed
   @li DofManager::generateGlobalIdSpace() has been called
   **

   @fn DofManager::getValue
   @param dof A Dof
   @return Returns the value of the given Dof, if it has been fixed

   This method throws an Exception if the Dof has not been fixed
   **

   @fn  DofManager::toString
   @return Returns the DofManager string
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

template<typename scalar>
inline const size_t DofManager<scalar>::isFixedId(void){
  return isFixed;
}

template<typename scalar>
inline size_t DofManager<scalar>::getGlobalId(const Dof& dof) const{
  return globalIdV[dof.getEntity() - first][dof.getType()];
}

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "DofManagerInclusion.h"

#endif
