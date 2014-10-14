#ifndef _FORMULATIONUPDATEOSRCSCALAR_H_
#define _FORMULATIONUPDATEOSRCSCALAR_H_

#include <map>

#include "SmallFem.h"
#include "FormulationBlock.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"

#include "DDMContextOSRC.h"

/**
   @class FormulationUpdateOSRCScalar
   @brief Update Formulation for FormulationOSRCScalar

   Update Formulation for FormulationOSRCScalar.

   Since this Formulation needs the volume solution (restriced at DDM border),
   the FormulationUpdateOSRCScalar::update() must be called @em before
   calling FormulationUpdateOSRCScalar::weak()
   or FormulationUpdateOSRCScalar::rhs() or an @em assembly procedure.
*/

class FormulationUpdateOSRCScalar: public FormulationBlock<Complex>{
 private:
  // Wavenumber //
  double k;

  // DDMContext //
  DDMContextOSRC* context;

  // Stuff for updating RHS //
  int NPade;
  const Basis*     basis;
  Quadrature*      gauss;
  GroupOfJacobian* jac;

  // Volume Solution (field and auxiliary) //
  std::map<Dof, Complex> solU;
  std::map<Dof, Complex> UPhi;
  std::vector<std::map<Dof, Complex> > solPhi;

  // Function Space & Domain //
  const FunctionSpace*  ffspace;
  const GroupOfElement* ddomain;

  // Pade //
  Complex C0;
  std::vector<Complex> A;
  std::vector<Complex> B;

  // Local Terms //
  TermFieldField<double>*       lGout;
  TermProjectionField<Complex>* lGin;
  TermProjectionField<Complex>* lC0;
  TermProjectionField<Complex>* lAB;

 public:
  FormulationUpdateOSRCScalar(DDMContextOSRC& context);

  virtual ~FormulationUpdateOSRCScalar(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);

 private:
  void resetUPhi(void);
  void getUPhi(void);
};

/**
   @fn FormulationUpdateOSRCScalar::FormulationUpdateOSRCScalar
   @param context A DDMContextOSRC

   Instantiates a new FormulationUpdateOSRCScalar with the given DDMContextOSRC
   **

   @fn FormulationUpdateOSRCScalar::~FormulationUpdateOSRCScalar
   Deletes this FormulationUpdateOSRCScalar
   **

   @fn FormulationUpdateOSRCScalar::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
