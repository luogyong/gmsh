#ifndef _FORMULATIONOSRCSCALAR_H_
#define _FORMULATIONOSRCSCALAR_H_

#include <map>

#include "SmallFem.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "FormulationCoupled.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationOSRCScalarOne.h"
#include "DDMContextOSRCScalar.h"

/**
   @class FormulationOSRCScalar
   @brief Scalar OSRC Formulation for DDM

   Scalar OSRC Formulation for DDM
 */

class FormulationOSRCScalarOne;

class FormulationOSRCScalar: public FormulationCoupled<Complex>{
 private:
  // DDMContext //
  DDMContextOSRCScalar* context;

  // Stuff for updating RHS //
  const Basis*               basis;
  const FunctionSpace*       ffspace;
  const FunctionSpace*       ffspaceG;
  Quadrature*                gaussFF;
  GroupOfJacobian*           jacFF;
  FormulationOSRCScalarOne*  formulationOne;

  // Local Terms //
  TermFieldField<double>*       localFF;
  TermGradGrad<double>*         localGG;
  TermProjectionField<Complex>* localPr;

  // Formulations //
  std::list<const FormulationBlock<Complex>*> fList;

 public:
  FormulationOSRCScalar(DDMContextOSRCScalar& context);

  virtual ~FormulationOSRCScalar(void);

  virtual
    const std::list<const FormulationBlock<Complex>*>&
                                               getFormulationBlocks(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);
};

/**
   @fn FormulationOSRCScalar::FormulationOSRCScalar
   @param context A DDMContextOSRCScalar

   Instantiates a new FormulationOSRCScalar with the given DDMContextOSRCScalar
   **

   @fn FormulationOSRCScalar::~FormulationOSRCScalar
   Deletes this FormulationOSRCScalar
   **

   @fn FormulationOSRCScalar::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
