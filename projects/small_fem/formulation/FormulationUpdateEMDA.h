#ifndef _FORMULATIONUPDATEEMDA_H_
#define _FORMULATIONUPDATEEMDA_H_

#include <map>

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "FormulationBlock.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "DDMContextEMDA.h"

/**
   @class FormulationUpdateEMDA
   @brief Update Formulation for FormulationEMDA

   Update Formulation for FormulationEMDA.

   Since this Formulation needs the volume solution (restriced at DDM border),
   the FormulationUpdateEMDA::update() must be called @em before
   calling FormulationUpdateEMDA::weak() or FormulationUpdateEMDA::rhs()
   or an @em assembly procedure.
 */

class FormulationUpdateEMDA: public FormulationBlock<Complex>{
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

  // Volume Solution //
  std::map<Dof, Complex> sol;

  // Function Space & Domain //
  const FunctionSpace*  fspace;
  const GroupOfElement* ddomain;

  // Local Terms //
  TermFieldField<double>*       lGout;
  TermProjectionField<Complex>* lGin;
  TermProjectionField<Complex>* lU;

 public:
  FormulationUpdateEMDA(DDMContextEMDA& context);

  virtual ~FormulationUpdateEMDA(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);
};

/**
   @fn FormulationUpdateEMDA::FormulationUpdateEMDA
   @param context A DDMContextEMDA

   Instantiates a new FormulationUpdateEMDA with the given DDMContextEMDA
   **

   @fn FormulationUpdateEMDA::~FormulationUpdateEMDA
   Deletes this FormulationUpdateEMDA
   **

   @fn FormulationUpdateEMDA::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
