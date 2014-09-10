///////////////////////////////////////////////
// Templates Implementations for DofManager: //
// Inclusion compilation model               //
//                                           //
// Damn you gcc: we want 'export' !          //
///////////////////////////////////////////////

#include <mpi.h>
#include <sstream>

#include "MPIDofMap.h"
#include "Exception.h"

template<typename scalar>
const size_t DofManager<scalar>::isFixed = 0 - 1; // Largest size_t

template<typename scalar>
const size_t DofManager<scalar>::isUndef = 0 - 2; // Second Largest size_t

template<typename scalar>
DofManager<scalar>::DofManager(bool isLocal){
  local = isLocal;
}

template<typename scalar>
DofManager<scalar>::DofManager(void){
  local = true;
}

template<typename scalar>
DofManager<scalar>::~DofManager(void){
}

template<typename scalar>
size_t DofManager<scalar>::getLocalSize(void) const{
  if(!globalIdM.empty())
    throw Exception("DofManager::getLocalSize() cannot be used before "
                    "DofManager::generateGlobalIdSpace() is called");
  else
    return nTotUnfixedLocalDof;
}

template<typename scalar>
size_t DofManager<scalar>::getGlobalSize(void) const{
  if(!globalIdM.empty())
    throw Exception("DofManager::getGlobalSize() cannot be used before "
                    "DofManager::generateGlobalIdSpace() is called");
  else
    return nTotUnfixedGlobalDof;
}

template<typename scalar>
void DofManager<scalar>::addToDofManager(const std::set<Dof>& dof){
  // Check if vector has been created //
  if(!globalIdV.empty())
    throw
      Exception
      ("DofManager: global id space generated -> can't add Dof");

  // Iterators //
  std::set<Dof>::const_iterator it  = dof.begin();
  std::set<Dof>::const_iterator end = dof.end();

  // Add to DofManager //
  for(; it != end; it++)
    globalIdM.insert(std::pair<Dof, size_t>(*it, 0));
}

template<typename scalar>
void DofManager<scalar>::
addToDofManager(const std::vector<std::vector<Dof> >& dof){
  // Check if vector has been created //
  if(!globalIdV.empty())
    throw
      Exception
      ("DofManager: global id space generated -> can't add Dof");

  // Add to DofManager //
  const size_t size = dof.size();
        size_t nDof;

  for(size_t i = 0; i < size; i++){
    nDof = dof[i].size();

    for(size_t j = 0; j < nDof; j++)
      globalIdM.insert(std::pair<Dof, size_t>(dof[i][j], 0));
  }
}

template<typename scalar>
void DofManager<scalar>::generateGlobalIdSpace(void){
  int nProc;
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  if(local || nProc == 1) // Even if global DofManager,
    localSpace();         // use local version if only ONE MPI process

  else
    globalSpace();

  vectorize();
  globalIdM.clear();
}

template<typename scalar>
void DofManager<scalar>::localSpace(void){
  // Populate map with global Id
  count(globalIdM);

  // Total Number of local (and global) unfixed Dof
  nTotUnfixedLocalDof  = globalIdM.size() - fixedDof.size();
  nTotUnfixedGlobalDof = globalIdM.size() - fixedDof.size();
}

template<typename scalar>
void DofManager<scalar>::globalSpace(void){
  // Get Global Maps (Global IDs and Fixed Dofs) //
  std::map<Dof, size_t> allGlobalIdM;
  std::map<Dof, scalar> allFixedDof;

  MPIDofMap<size_t>::getGlobalMap(globalIdM, allGlobalIdM);
  MPIDofMap<scalar>::getGlobalMap(fixedDof , allFixedDof);

  // Tag Fixed Dofs into allGlobalIdM //
  retag(allGlobalIdM, allFixedDof);

  // Populate global map with global Id //
  count(allGlobalIdM);

  // Populate local map with global Id from global map //
  // Iterate on LOCAL Dofs
  std::map<Dof, size_t>::iterator it  = globalIdM.begin();
  std::map<Dof, size_t>::iterator end = globalIdM.end();
  std::map<Dof, size_t>::iterator find;

  for(; it != end; it++){
    // Find current Dof counterpart in GLOBAL map
    find = allGlobalIdM.find(it->first);

    // Set current Dof value in LOCAL map
    it->second = find->second;
  }

  // Total Number of local (and global) unfixed Dof
  nTotUnfixedLocalDof  =    globalIdM.size() -    fixedDof.size();
  nTotUnfixedGlobalDof = allGlobalIdM.size() - allFixedDof.size();
}

template<typename scalar>
void DofManager<scalar>::count(std::map<Dof, size_t>& dof){
  // Get Iterators //
  std::map<Dof, size_t>::iterator end = dof.end();
  std::map<Dof, size_t>::iterator it  = dof.begin();

  // Start Counter //
  size_t id = 0;

  // Iterate on Dofs //
  for(; it != end; it++){
    // Check if unknown
    if(it->second != isFixed){
      // If not: count it !
      it->second = id;
      id++;
    }
  }
}

template<typename scalar>
void DofManager<scalar>::
retag(std::map<Dof, size_t>& dof, std::map<Dof, scalar>& fix){
  // Iterate on Dofs //
  std::map<Dof, size_t>::iterator it  = dof.begin();
  std::map<Dof, size_t>::iterator end = dof.end();

  for(; it != end; it++)
    // Look for current Dof into Fixed map
    if(fix.count(it->first) == 1)
      // If Dof is fond, tag it as fixed
      it->second = isFixed;
}

template<typename scalar>
void DofManager<scalar>::vectorize(void){
  // Get Data //
  std::map<Dof, size_t>::iterator end = globalIdM.end();
  std::map<Dof, size_t>::iterator it  = globalIdM.begin();

  // Take the last element *IN* map
  end--;

  first = it->first.getEntity();
  last  = end->first.getEntity();

  // Reset 'end': take the first element *OUTSIDE* map
  end++;


  // Alloc //
  const size_t sizeV = last - first + 1;
  globalIdV.resize(sizeV);

  // Populate //
  size_t nDof;
  size_t max;
  std::map<Dof, size_t>::iterator currentEntity = it;

  // Iterate on vector
  for(size_t i = 0; i < sizeV; i++){
    // No dof found
    nDof = 0;

    // 'currentEntity - first' matches 'i' ?
    if(it != end && currentEntity->first.getEntity() - first == i)
      // Count all elements with same entity in map
      for(; it !=end &&
            currentEntity->first.getEntity() == it->first.getEntity(); it++)
        nDof++; // New Dof found

    // Dof with Biggest type in this 'Same Entity Range'
    it--;
    max = it->first.getType();
    it++;

    //////////////////////////////////////////////////////////////////
    // Here we have the following configuration:                    //
    // ----------------------------------------                     //
    //                                                              //
    // itrators:      currentEntity                      it         //
    // variable:             |      nDof                 |          //
    //                       <----------------->         |          //
    //                       |                 |         |          //
    //                       v                 v         v          //
    // map: ... ; (3, 4) ; (4, 0) ; (4, 2) ; (4, 10) ; (6, 0) ; ... //
    //                                            ^                 //
    //                                            |                 //
    // variable:                                 max                //
    //////////////////////////////////////////////////////////////////

    // Alloc if Dofs were found
    if(nDof)
      // Space for maximum type in this rang of Dof
      // Up to now, values are undefined
      globalIdV[i].resize(max + 1, isUndef);

    // Add globalIds in vector for this range of Dof
    for(; currentEntity != it; currentEntity++)
      globalIdV[i][currentEntity->first.getType()] = currentEntity->second;

    // Now currentEntity is equal to it and we can work on the next entity
  }
}

template<typename scalar>
std::pair<bool, size_t> DofManager<scalar>::findSafe(const Dof& dof) const{
  // Is globabId Vector allocated ?-
  if(globalIdV.empty())
    throw
      Exception
      ("Cannot get Dof %s ID, since global ID space has not been generated",
       dof.toString().c_str());

  // Is 'dof' in globalIdV range ?
  size_t tmpEntity = dof.getEntity();

  if(tmpEntity < first || tmpEntity > last)
    return std::pair<bool, size_t>(false, 42);

  // Offset Entity & Get Type
  const size_t entity = tmpEntity - first;
  const size_t type   = dof.getType();

  // Look for Entity in globalIdV
  const size_t nDof = globalIdV[entity].size();

  size_t globalId;

  if(nDof > 0 && type <= nDof){
    // If we have Dofs associated to this Entity,
    // get the requested Type and fetch globalId
    globalId = globalIdV[entity][type];

    // If globalId is defined return it,
    if(globalId != isUndef)
      return std::pair<bool, size_t>(true, globalId);

    // Else, no Dof and return false
    else
      return std::pair<bool, size_t>(false, 42);
  }

  else
    // If no Dof, return false
    return std::pair<bool, size_t>(false, 42);
}

template<typename scalar>
size_t DofManager<scalar>::getGlobalIdSafe(const Dof& dof) const{
  const std::pair<bool, size_t> search = findSafe(dof);

  if(!search.first)
    throw
      Exception("Their is no Dof %s", dof.toString().c_str());

  else
    return search.second;
}

template<typename scalar>
bool DofManager<scalar>::fixValue(const Dof& dof, scalar value){
  // Check if map is still their
  if(globalIdM.empty())
    return false;

  // Get *REAL* Dof
  const std::map<Dof, size_t>::iterator it = globalIdM.find(dof);

  // Check if 'dof' exists
  if(it == globalIdM.end())
    return false; // 'dof' doesn't exist

  // If 'dof' exists: it becomes fixed
  fixedDof.insert(std::pair<Dof, scalar>(it->first, value));
  it->second = isFixed;
  return true;
}

template<typename scalar>
scalar DofManager<scalar>::getValue(const Dof& dof) const{
  typename std::map<Dof, scalar>::const_iterator end = fixedDof.end();
  typename std::map<Dof, scalar>::const_iterator  it = fixedDof.find(dof);

  if(it == end)
    throw Exception("Dof %s not fixed", dof.toString().c_str());

  return it->second;
}

template<typename scalar>
std::string DofManager<scalar>::toString(void) const{
  if(!globalIdM.empty())
    return toStringFromMap();

  else
    return toStringFromVec();
}

template<typename scalar>
std::string DofManager<scalar>::toStringFromMap(void) const{
  std::stringstream s;

  std::map<Dof, size_t>::const_iterator end = globalIdM.end();
  std::map<Dof, size_t>::const_iterator it  = globalIdM.begin();

  for(; it != end; it++){
    s << "("  << it->first.toString() << ": ";

    if(it->second == isFixed)
      s << fixedDof.find(it->first)->second << " -- Fixed value";

    else
      s << it->second                       << " -- Global ID";

    s << ")"  << std::endl;
  }

  return s.str();
}

template<typename scalar>
std::string DofManager<scalar>::toStringFromVec(void) const{
  const size_t sizeV = globalIdV.size();

  std::stringstream s;
  size_t nDof;
  std::pair<bool, size_t> search;

  for(size_t entity = 0; entity < sizeV; entity++){
    nDof = globalIdV[entity].size();

    for(size_t type = 0; type < nDof; type++){
      Dof dof(entity + first, type);
      search = findSafe(dof);

      if(search.first){
        s << "("  << dof.toString() << ": ";

        if(search.second == isFixed)
          s << fixedDof.find(dof)->second << " -- Fixed value";

        else
          s << search.second              << " -- Global ID";

        s << ")"  << std::endl;
      }
    }
  }

  return s.str();
}
