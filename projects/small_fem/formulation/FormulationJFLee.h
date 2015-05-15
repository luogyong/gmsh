#ifndef _FORMULATIONJFLEE_H_
#define _FORMULATIONJFLEE_H_

#include <map>

#include "SmallFem.h"
#include "GroupOfElement.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "TermCurlCurl.h"
#include "TermGradGrad.h"
#include "TermFieldField.h"
#include "TermProjectionGrad.h"

#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationJFLeeOne.h"
#include "FormulationCoupled.h"
#include "DDMContextJFLee.h"

/**
   @class FormulationJFLee
   @brief JFLee Formulation for DDM

   JFLee Formulation for DDM
 */

class FormulationJFLeeOne;

class FormulationJFLee: public FormulationCoupled<Complex>{
 private:
  // DDMContext //
  DDMContextJFLee* context;

  // Stuff for updating RHS //
  const Basis*               basisE;
  const FunctionSpaceVector* field;
  const FunctionSpaceVector* fieldG;
  Quadrature*                gauss;
  GroupOfJacobian*           jac;
  FormulationJFLeeOne*       formulationOne;

  // Local Terms //
  TermProjectionGrad<Complex>* proj;
  TermGradGrad<double>*        termPE;
  TermGradGrad<double>*        termRP;
  TermGradGrad<double>*        termPP;
  TermGradGrad<double>*        termEP;
  TermCurlCurl<double>*        term22;
  TermFieldField<double>*      termRR;
  TermGradGrad<double>*        termPR;

  // Formulations //
  std::list<FormulationBlock<Complex>*> fList;

 public:
  FormulationJFLee(DDMContextJFLee& context);

  virtual ~FormulationJFLee(void);

  virtual
    const std::list<FormulationBlock<Complex>*>&
                                               getFormulationBlocks(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);
};

/**
   @fn FormulationJFLee::FormulationJFLee
   @param context A DDMContextJFLee

   Instantiates a new FormulationJFLee with the given DDMContextJFLee
   **

   @fn FormulationJFLee::~FormulationJFLee
   Deletes this FormulationJFLee
   **

   @fn FormulationJFLee::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
*/

#endif
