#ifndef _FORMULATIONEMDA_H_
#define _FORMULATIONEMDA_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "FormulationBlock.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Term.h"

#include "DDMContextEMDA.h"

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
  DDMContextEMDA* context;

  // Stuff for updating RHS //
  const Basis*     basis;
  Quadrature*      gauss;
  GroupOfJacobian* jac;

  // Function Space & Domain //
  const FunctionSpace*  fspace;
  const GroupOfElement* ddomain;

  // Local Terms //
  Term<double>*  localLHS;
  Term<Complex>* localRHS;

 public:
  FormulationEMDA(DDMContextEMDA& context);

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
   @param context A DDMContextEMDA

   Instantiates a new FormulationEMDA with the given DDMContextEMDA
   **

   @fn FormulationEMDA::~FormulationEMDA
   Deletes this FormulationEMDA
   **

   @fn FormulationEMDA::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
