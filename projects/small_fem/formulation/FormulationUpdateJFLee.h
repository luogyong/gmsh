#ifndef _FORMULATIONUPDATEJFLEE_H_
#define _FORMULATIONUPDATEJFLEE_H_

#include <map>

#include "SmallFem.h"
#include "FormulationBlock.h"
#include "FunctionSpaceVector.h"

#include "TermProjectionGrad.h"
#include "TermGradGrad.h"

#include "Quadrature.h"
#include "GroupOfJacobian.h"

#include "DDMContextJFLee.h"

/**
   @class FormulationUpdateJFLee
   @brief Update Formulation for FormulationJFLee

   Update Formulation for FormulationJFLee.

   Since this Formulation needs the volume solution (restriced at DDM border),
   the FormulationUpdateJFLee::update() must be called @em before
   calling FormulationUpdateJFLee::weak() or FormulationUpdateJFLee::rhs()
   or an @em assembly procedure.
*/

class FormulationUpdateJFLee: public FormulationBlock<Complex>{
 private:
  // DDMContext //
  DDMContextJFLee* context;

  // Stuff for updating RHS //
  const Basis*     basis;
  Quadrature*      gauss;
  GroupOfJacobian* jac;

  // Volume Solution (field and auxiliary) //
  std::map<Dof, Complex> phi;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  fPhi;
  const GroupOfElement* ddomain;

  // Local Terms //
  TermGradGrad<double>*       lGout;
  TermProjectionGrad<Complex>* lGin;
  TermProjectionGrad<Complex>* lPhi;

 public:
  FormulationUpdateJFLee(DDMContextJFLee& context);

  virtual ~FormulationUpdateJFLee(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);
};

/**
   @fn FormulationUpdateJFLee::FormulationUpdateJFLee
   @param context A DDMContextJFLee

   Instantiates a new FormulationUpdateJFLee with the given DDMContextJFLee
   **

   @fn FormulationUpdateJFLee::~FormulationUpdateJFLee
   Deletes this FormulationUpdateJFLee
   **

   @fn FormulationUpdateJFLee::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
