#ifndef _FORMULATIONEMDA_H_
#define _FORMULATIONEMDA_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "FormulationBlock.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "DDMContext.h"

/**
   @class FormulationEMDA
   @brief EMDA Formulation for DDM

   EMDA Formulation for DDM
 */

class FormulationEMDA: public FormulationBlock<Complex>{
 private:
  // Wavenumber & Chi //
  double k;
  double chi;

  // DDMContext //
  DDMContext* context;

  // Stuff for updating RHS //
  const Basis*     basis;
  Quadrature*      gauss;
  GroupOfJacobian* jac;

  // Function Space & Domain //
  const FunctionSpace*  fspace;
  const GroupOfElement* ddomain;

  // Local Terms //
  TermFieldField<double>*       localLHS;
  TermProjectionField<Complex>* localRHS;

 public:
  FormulationEMDA(DDMContext& context);

  virtual ~FormulationEMDA(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);
};

/**
   @fn FormulationEMDA::FormulationEMDA
   @param context A DDMContext

   Instantiates a new FormulationEMDA with the given DDMContext
   **

   @fn FormulationEMDA::~FormulationEMDA
   Deletes this FormulationEMDA
   **

   @fn FormulationEMDA::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
