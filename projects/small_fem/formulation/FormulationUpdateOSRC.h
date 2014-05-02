#ifndef _FORMULATIONUPDATEOSRC_H_
#define _FORMULATIONUPDATEOSRC_H_

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
   @class FormulationUpdateOSRC
   @brief Update Formulation for FormulationOSRC

   Update Formulation for FormulationOSRC.

   Since this Formulation needs the volume solution (restriced at DDM border),
   the FormulationUpdateOSRC::update() must be called @em before
   calling FormulationUpdateOSRC::weak() or FormulationUpdateOSRC::rhs()
   or an @em assembly procedure.
*/

class FormulationUpdateOSRC: public FormulationBlock<Complex>{
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
  FormulationUpdateOSRC(DDMContextOSRC& context);

  virtual ~FormulationUpdateOSRC(void);

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
   @fn FormulationUpdateOSRC::FormulationUpdateOSRC
   @param context A DDMContextOSRC

   Instantiates a new FormulationUpdateOSRC with the given DDMContextOSRC
   **

   @fn FormulationUpdateOSRC::~FormulationUpdateOSRC
   Deletes this FormulationUpdateOSRC
   **

   @fn FormulationUpdateOSRC::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
