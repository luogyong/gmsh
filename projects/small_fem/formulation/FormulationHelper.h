#ifndef _FORMULATIONHELPER_H_
#define _FORMULATIONHELPER_H_

#include "FunctionSpace.h"
#include "FunctionSpaceScalar.h"
#include "GroupOfElement.h"
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

  static void initDofMap(const std::vector<const FunctionSpaceScalar*>& fs,
                         const GroupOfElement& goe,
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
   restricted to the given GroupOfElement.

   The Dof values are set to zero.
   **

   @fn FormulationHelper::initDofMap(const std::vector<const FunctionSpaceScalar*>&,const GroupOfElement&,std::vector<std::map<Dof, std::complex<double> > >&)
   @param fs A vector of FunctionSpaceScalar
   @param goe A GroupOfElement
   @param data A vector of Dof -- Value map

   Populates the given map[i] with the Dof%s of the given FunctionSpace[i]
   restricted to the given GroupOfElement, for all i in the given vectors.

   The Dof values are set to zero.
 */

#endif
