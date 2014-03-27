#ifndef _FORMULATIONOSRC_H_
#define _FORMULATIONOSRC_H_

#include <map>

#include "SmallFem.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"

#include "CoupledFormulation.h"

/**
   @class FormulationOSRC
   @brief OSRC Formulation for DDM

   OSRC Formulation for DDM
 */

class FormulationOSRC: public CoupledFormulation<Complex>{
 private:
  // Local Terms //
  TermFieldField<double>*       localFF;
  TermGradGrad<double>*         localGG;
  TermProjectionField<Complex>* localPr;

  // Formulations //
  std::list<const Formulation<Complex>*> fList;

 public:
  FormulationOSRC(const GroupOfElement& domain,
                  const FunctionSpaceScalar& field,
                  const std::vector<const FunctionSpaceScalar*>& aux,
                  double k,
                  Complex keps,
                  int NPade,
                  const std::map<Dof, Complex>& ddmDof);

  virtual ~FormulationOSRC(void);

  virtual
    const std::list<const Formulation<Complex>*>& getFormulations(void) const;

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
   @param domain A GroupOfElement for the domain
   @param field FunctionSpace for the computed field
   @param aux A set of auxiliary FunctionSpace for Pade denominator
   @param k A real number
   @param keps A complex number
   @param NPade Number of approximation terms for Pade
   @param ddmDof A map with the DDM Dof%s and their associated values

   Instantiates a new FormulationOSRC with wavenumber k
   and complexified wavenumber keps.
   The DDM Dof%s are given by ddmDof.

   It is important mentioning that aux.size() and NPade must be equal.
   **

   @fn FormulationOSRC::~FormulationOSRC
   Deletes this FormulationOSRC
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
