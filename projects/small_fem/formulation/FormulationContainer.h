#ifndef _FORMULATIONCONTAINER_H_
#define _FORMULATIONCONTAINER_H_

#include "FormulationCoupled.h"

/**
   @class FormulationContainer
   @brief A container for FormulationBlock%s

   A container for FormulationBlock%s
 */

template<typename scalar>
class FormulationContainer: public FormulationCoupled<scalar>{
 private:
  std::list<FormulationBlock<scalar>*> fList;

 public:
  FormulationContainer(void);

  virtual ~FormulationContainer(void);

  void addFormulation(FormulationBlock<scalar>&   formulation);
  void addFormulation(FormulationCoupled<scalar>& formulation);

  virtual
    const std::list<FormulationBlock<scalar>*>&
                                               getFormulationBlocks(void) const;

  virtual void update(void);
};

/**
   @fn FormulationContainer::FormulationContainer
   Instantiates a new and empty FormulationContainer
   **

   @fn FormulationContainer::~FormulationContainer
   Deletes this FormulationContainer
   **

   @fn FormulationContainer::addFormulation(FormulationBlock<scalar>&)
   @param formulation A FormulationBlock
   Adds the given FormulationBlock into this FormulationContainer
   **

   @fn FormulationContainer::addFormulation(FormulationCoupled<scalar>&)
   @param formulation A FormulationCoupled
   Adds every FormulationBlock of the given FormulationCoupled
   into this FormulationContainer
   **

   @fn FormulationContainer::update
   Updates every FormulationBlock in this FormulationContainer
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationContainerInclusion.h"

#endif
