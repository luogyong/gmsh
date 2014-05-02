#ifndef _FORMULATIONOO2_H_
#define _FORMULATIONOO2_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "FormulationBlock.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "DDMContext.h"

/**
   @class FormulationOO2
   @brief OO2 Formulation for DDM

   OO2 Formulation for DDM
 */

class FormulationOO2: public FormulationBlock<Complex>{
 private:
  // a & b //
  Complex a;
  Complex b;

  // DDMContext //
  DDMContext* context;

  // Stuff for updating RHS //
  const Basis*     basis;
  Quadrature*      gaussFF;
  GroupOfJacobian* jacFF;

  // Function Space & Domain //
  const FunctionSpace*  fspace;
  const GroupOfElement* ddomain;

  // Local Terms //
  TermFieldField<double>*       localTermsFF;
  TermGradGrad<double>*         localTermsGG;
  TermProjectionField<Complex>* localTermsPr;

 public:
  FormulationOO2(DDMContext& context);

  virtual ~FormulationOO2(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);
};

/**
   @fn FormulationOO2::FormulationOO2
   @param context A DDMContext

   Instantiates a new FormulationOO2 with the given DDMContext
   **

   @fn FormulationOO2::~FormulationOO2
   Deletes this FormulationOO2
   **

   @fn FormulationOO2::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
