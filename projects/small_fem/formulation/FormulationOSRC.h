#ifndef _FORMULATIONOSRC_H_
#define _FORMULATIONOSRC_H_

#include <map>

#include "SmallFem.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "FormulationCoupled.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationOSRCOne.h"
#include "DDMContext.h"

/**
   @class FormulationOSRC
   @brief OSRC Formulation for DDM

   OSRC Formulation for DDM
 */

class FormulationOSRCOne;

class FormulationOSRC: public FormulationCoupled<Complex>{
 private:
  // DDMContext //
  DDMContext* context;

  // Stuff for updating RHS //
  const Basis*         basis;
  const FunctionSpace* field;
  Quadrature*          gaussFF;
  GroupOfJacobian*     jacFF;
  FormulationOSRCOne*  formulationOne;

  // Local Terms //
  TermFieldField<double>*       localFF;
  TermGradGrad<double>*         localGG;
  TermProjectionField<Complex>* localPr;

  // Formulations //
  std::list<const FormulationBlock<Complex>*> fList;

 public:
  FormulationOSRC(DDMContext& context);

  virtual ~FormulationOSRC(void);

  virtual
    const std::list<const FormulationBlock<Complex>*>&
                                               getFormulationBlocks(void) const;

  virtual bool isBlock(void) const;
  virtual void update(void);

 private:
  static double pade_aj(int j, int N);
  static double pade_bj(int j, int N);

 public:
  static Complex padeC0(int N, double theta);
  static Complex padeAj(int j, int N, double theta);
  static Complex padeBj(int j, int N, double theta);
};

/**
   @fn FormulationOSRC::FormulationOSRC
   @param context A DDMContext

   Instantiates a new FormulationOSRC with the given DDMContext
   **

   @fn FormulationOSRC::~FormulationOSRC
   Deletes this FormulationOSRC
   **

   @fn FormulationOSRC::update
   Updates the DDM Dof%s values from the DDMContext given at construction time
   **

   @fn FormulationOSRC::padeC0
   @param N A natural number
   @param theta A real number
   @return Returns the C0 term of the OSRC Pade approximation for N terms
   and complex rotation of theta.
   **

   @fn FormulationOSRC::padeAj
   @param j A natural number
   @param N A natural number
   @param theta A real number
   @return Returns the jth A term of the OSRC Pade approximation for N terms
   and complex rotation of theta.
   **

   @fn FormulationOSRC::padeBj
   @param j A natural number
   @param N A natural number
   @param theta A real number
   @return Returns the jth B term of the OSRC Pade approximation for N terms
   and complex rotation of theta.
   **
*/

#endif
