#ifndef _FORMULATIONUPDATEEMDA_H_
#define _FORMULATIONUPDATEEMDA_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "FormulationBlock.h"

#include "DDMContext.h"

/**
   @class FormulationUpdateEMDA
   @brief Update Formulation for FormulationEMDA

   Update Formulation for FormulationEMDA
 */

class FormulationUpdateEMDA: public FormulationBlock<Complex>{
 private:
  // Wavenumber & Chi //
  double k;
  double chi;

  // Function Space & Domain //
  const FunctionSpace*  fspace;
  const GroupOfElement* ddomain;

  // Local Terms //
  TermFieldField<double>*       lGout;
  TermProjectionField<Complex>* lGin;
  TermProjectionField<Complex>* lU;

 public:
  FormulationUpdateEMDA(DDMContext& context);

  virtual ~FormulationUpdateEMDA(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationUpdateEMDA::FormulationUpdateEMDA
   @param context A DDMContext

   Instantiates a new FormulationUpdateEMDA with the given DDMContext
   **

   @fn FormulationUpdateEMDA::~FormulationUpdateEMDA
   Deletes this FormulationUpdateEMDA
*/

#endif
