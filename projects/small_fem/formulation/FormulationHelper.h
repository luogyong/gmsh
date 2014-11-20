#ifndef _FORMULATIONHELPER_H_
#define _FORMULATIONHELPER_H_

#include "FunctionSpace.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "GroupOfElement.h"
#include "Exception.h"
#include "Dof.h"

#include "SmallFem.h"
#include <vector>
#include <map>

/**
   @class FormulationHelper
   @brief A bunch of static method for helping Formulation%s

   A bunch of static method for helping Formulation%s
*/

class FormulationHelper{
 public:
   FormulationHelper(void);
  ~FormulationHelper(void);

  static void initDofMap(const FunctionSpace& fs,
                         const GroupOfElement& goe,
                         std::map<Dof, Complex>& data);

  static void initDofMap(const std::vector<const FunctionSpace*>& fs,
                         const GroupOfElement& goe,
                         std::vector<std::map<Dof, Complex> >& data);

  static void initDofMap(const std::vector<const FunctionSpaceScalar*>& fs,
                         const GroupOfElement& goe,
                         std::vector<std::map<Dof, Complex> >& data);

  static void initDofMap(const std::vector<const FunctionSpaceVector*>& fs,
                         const GroupOfElement& goe,
                         std::vector<std::map<Dof, Complex> >& data);

 private:
  template<typename T>
    static size_t getSize(const std::vector<const T*>& fs,
                          std::vector<std::map<Dof, Complex> >& data);
};

/**
   @fn FormulationHelper::FormulationHelper
   Instanciates a new FormulationHelper (unneeded since this is a static class)
   **

   @fn FormulationHelper::~FormulationHelper
   Deletes this FormulationHelper
   **

   @fn FormulationHelper::initDofMap(const FunctionSpace&,const GroupOfElement&,std::map<Dof, std::complex<double> >&)
   @param fs A FunctionSpace
   @param goe A GroupOfElement
   @param data A Dof -- Value map

   Populates the given map with the Dof%s of the given FunctionSpace
   restricted to the 'goe' GroupOfElement.

   The Dof values are set to zero.
   **

   @fn FormulationHelper::initDofMap(const std::vector<const FunctionSpace*>&,const GroupOfElement&,std::vector<std::map<Dof, std::complex<double> > >&)
   @param fs A vector of FunctionSpace
   @param goe A GroupOfElement
   @param data A vector of Dof -- Value map

   Populates the given map[i] with the Dof%s of the given FunctionSpace[i]
   restricted to the 'goe' GroupOfElement,
   for all i in the given vectors.

   The Dof values are set to zero.

   @fn FormulationHelper::initDofMap(const std::vector<const FunctionSpaceScalar*>&,const GroupOfElement&,std::vector<std::map<Dof, std::complex<double> > >&)
   @param fs A vector of FunctionSpaceScalar
   @param goe A GroupOfElement
   @param data A vector of Dof -- Value map

   Populates the given map[i] with the Dof%s of the given FunctionSpace[i]
   restricted to the 'goe' GroupOfElement,
   for all i in the given vectors.

   The Dof values are set to zero.
   **

   @fn FormulationHelper::initDofMap(const std::vector<const FunctionSpaceVector*>&,const GroupOfElement&,std::vector<std::map<Dof, std::complex<double> > >&)
   @param fs A vector of FunctionSpaceVector
   @param goe A GroupOfElement
   @param data A vector of Dof -- Value map

   Populates the given map[i] with the Dof%s of the given FunctionSpace[i]
   restricted to the 'goe' GroupOfElement,
   for all i in the given vectors.

   The Dof values are set to zero.
 */

///////////////////////
// Template Function //
///////////////////////
template<typename T>
size_t FormulationHelper::getSize(const std::vector<const T*>& fs,
                                  std::vector<std::map<Dof, Complex> >& data){
  const size_t size = data.size();

  if(size != fs.size())
    throw Exception("FormulationHelper::initDofMap: %s",
                    "vectors must have the same size");
  return size;
}

#endif
