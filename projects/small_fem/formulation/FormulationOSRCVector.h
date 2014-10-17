#ifndef _FORMULATIONOSRCVECTOR_H_
#define _FORMULATIONOSRCVECTOR_H_

#include <map>

#include "SmallFem.h"
#include "FunctionSpace.h"
#include "GroupOfElement.h"

#include "TermProjectionGrad.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "TermCurlCurl.h"

#include "FormulationCoupled.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationOSRCVectorThree.h"
#include "DDMContextOSRCVector.h"

/**
   @class FormulationOSRCVector
   @brief Vector OSRC Formulation for DDM

   Vector OSRC Formulation for DDM
 */

class FormulationOSRCVectorThree;

class FormulationOSRCVector: public FormulationCoupled<Complex>{
 private:
  // DDMContext //
  DDMContextOSRCVector* context;

  // Stuff for updating RHS //
  const Basis*                basisV;
  const FunctionSpace*        field;
  Quadrature*                 gauss;
  GroupOfJacobian*            jac;
  FormulationOSRCVectorThree* formulationThree;

  // Local Terms //
  TermProjectionGrad<Complex>* RHS;
  TermGradGrad<double>*        GG;
  TermGradGrad<double>*        dFG;
  TermGradGrad<double>*        GdF;
  TermCurlCurl<double>*        CC;
  TermFieldField<double>*      FF;

  // Formulations //
  std::list<const FormulationBlock<Complex>*> fList;

 public:
  FormulationOSRCVector(DDMContextOSRCVector& context);

  virtual ~FormulationOSRCVector(void);

  virtual
    const std::list<const FormulationBlock<Complex>*>&
                                               getFormulationBlocks(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);
};

/**
   @fn FormulationOSRCVector::FormulationOSRCVector
   @param context A DDMContextOSRCVector

   Instantiates a new FormulationOSRCVector with the given DDMContextOSRCVector
   **

   @fn FormulationOSRCVector::~FormulationOSRCVector
   Deletes this FormulationOSRCVector
   **

   @fn FormulationOSRCVector::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
