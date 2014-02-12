#ifndef _FORMULATIONOO2_H_
#define _FORMULATIONOO2_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "Formulation.h"

/**
   @class FormulationOO2
   @brief OO2 Formulation for DDM

   OO2 Formulation for DDM
 */

class FormulationOO2: public Formulation<Complex>{
 private:
  // a & b //
  Complex a;
  Complex b;

  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField*               localTermsFF;
  TermGradGrad*                 localTermsGG;
  TermProjectionField<Complex>* localTermsPr;

 public:
  FormulationOO2(const GroupOfElement& domain,
                 const FunctionSpaceScalar& fs,
                 Complex a,
                 Complex b,
                 const std::map<Dof, Complex>& ddmDof);

  virtual ~FormulationOO2(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationOO2::FormulationOO2
   @param domain A GroupOfElement for the domain
   @param fs A FunctionSpace for both unknown and test fields
   @param a A complex number
   @param b A complex number
   @param ddmDof A map with the DDM Dof%s and their associated values

   Instantiates a new FormulationEMDA with parameters a and b.
   The DDM Dof%s are given by ddmDof.
   **

   @fn FormulationOO2::~FormulationOO2
   Deletes this FormulationOO2
*/

#endif
