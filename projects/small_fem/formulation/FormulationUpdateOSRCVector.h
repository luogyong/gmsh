#ifndef _FORMULATIONUPDATEOSRCVECTOR_H_
#define _FORMULATIONUPDATEOSRCVECTOR_H_

#include <map>

#include "SmallFem.h"
#include "FormulationBlock.h"

#include "TermProjectionGrad.h"
#include "TermGradGrad.h"

#include "Quadrature.h"
#include "GroupOfJacobian.h"

#include "DDMContextOSRCVector.h"

/**
   @class FormulationUpdateOSRCVector
   @brief Update Formulation for FormulationOSRCVector

   Update Formulation for FormulationOSRCVector.

   Since this Formulation needs the volume solution (restriced at DDM border),
   the FormulationUpdateOSRCVector::update() must be called @em before
   calling FormulationUpdateOSRCVector::weak()
   or FormulationUpdateOSRCVector::rhs() or an @em assembly procedure.
*/

class FormulationUpdateOSRCVector: public FormulationBlock<Complex>{
 private:
  // Wavenumber //
  Complex twoJOverK;

  // DDMContext //
  DDMContextOSRCVector* context;

  // Stuff for updating RHS //
  int NPade;
  const Basis*     basis;
  Quadrature*      gauss;
  GroupOfJacobian* jac;

  // Volume Solution (auxiliary) //
  std::map<Dof, Complex>               solR;
  std::vector<std::map<Dof, Complex> > solPhi;

  // Function Space & Domain //
  const FunctionSpace*  ffspace;
  const GroupOfElement* ddomain;

  // Auxiliary FunctionSpace //
  const FunctionSpaceVector*                     fR;
  const std::vector<const FunctionSpaceVector*>* fPhi;

  // Pade //
  Complex R0;
  std::vector<Complex> A;
  std::vector<Complex> B;

  // Local Terms //
  TermGradGrad<double>*        lGout;
  TermProjectionGrad<Complex>* lGin;
  TermProjectionGrad<Complex>* lR;
  TermProjectionGrad<Complex>* lPhi;

 public:
  FormulationUpdateOSRCVector(DDMContextOSRCVector& context);

  virtual ~FormulationUpdateOSRCVector(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);

 private:
  void getAllPhi(void);
};

/**
   @fn FormulationUpdateOSRCVector::FormulationUpdateOSRCVector
   @param context A DDMContextOSRCVector

   Instantiates a new FormulationUpdateOSRCVector
   with the given DDMContextOSRCVector
   **

   @fn FormulationUpdateOSRCVector::~FormulationUpdateOSRCVector
   Deletes this FormulationUpdateOSRCVector
   **

   @fn FormulationUpdateOSRCVector::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
