#ifndef _FORMULATIONUPDATEOO2_H_
#define _FORMULATIONUPDATEOO2_H_

#include <map>

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermProjectionGrad.h"
#include "TermFieldField.h"
#include "FormulationBlock.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "DDMContextOO2.h"

/**
   @class FormulationUpdateOO2
   @brief Update Formulation for FormulationOO2

   Update Formulation for FormulationOO2.

   Since this Formulation needs the volume solution (restriced at DDM border),
   the FormulationUpdateOO2::update() must be called @em before
   calling FormulationUpdateOO2::weak() or FormulationUpdateOO2::rhs()
   or an @em assembly procedure.
 */

class FormulationUpdateOO2: public FormulationBlock<Complex>{
 private:
  // a & b //
  Complex a;
  Complex b;

  // DDMContext //
  DDMContextOO2* context;

  // Stuff for updating RHS //
  const Basis*     basis;
  Quadrature*      gaussFF;
  Quadrature*      gaussGG;
  GroupOfJacobian* jacFF;
  GroupOfJacobian* jacGG;

  // Volume Solution //
  std::map<Dof, Complex> sol;

  // Function Space & Domain //
  const FunctionSpace*  fspace;
  const FunctionSpace*  fspaceG;
  const GroupOfElement* ddomain;

  // Local Terms //
  TermFieldField<double>*       lGout;
  TermProjectionField<Complex>* lGin;
  TermProjectionField<Complex>* lU;
  TermProjectionGrad<Complex>*  lDU;

 public:
  FormulationUpdateOO2(DDMContextOO2& context);

  virtual ~FormulationUpdateOO2(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);
};

/**
   @fn FormulationUpdateOO2::FormulationUpdateOO2
   @param context A DDMContextOO2

   Instantiates a new FormulationUpdateOO2 with the given DDMContextOO2
   **

   @fn FormulationUpdateOO2::~FormulationUpdateOO2
   Deletes this FormulationUpdateOO2
   **

   @fn FormulationUpdateOO2::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
