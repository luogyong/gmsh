#ifndef _FORMULATIONDUMMY_H_
#define _FORMULATIONDUMMY_H_

#include "FormulationCoupled.h"

/**
   @class FormulationDummy
   @brief Dummy Formulation

   A Formulation that does nothing
 */

template<typename scalar>
class FormulationDummy: public FormulationCoupled<scalar>{
 private:
  std::list<const FormulationBlock<scalar>*> fList;

 public:
  FormulationDummy(void);

  virtual ~FormulationDummy(void);

  virtual
    const std::list<const FormulationBlock<scalar>*>&
                                               getFormulationBlocks(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationDummy::FormulationDummy
   Instantiates a new FormulationDummy
   **

   @fn FormulationDummy::~FormulationDummy
   Deletes this FormulationDummy
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationDummyInclusion.h"

#endif
