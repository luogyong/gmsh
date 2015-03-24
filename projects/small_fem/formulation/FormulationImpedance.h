#ifndef _FORMULATIONIMPEDANCE_H_
#define _FORMULATIONIMPEDANCE_H_

#include "SmallFem.h"
#include "Term.h"
#include "FunctionSpace.h"
#include "FormulationBlock.h"

/**
   @class FormulationImpedance
   @brief Standard impedance boundary condition (SIBC) formulation

   Weak formulation for the
   standard impedance boundary condition (SIBC) formulation
 */

class FormulationImpedance: public FormulationBlock<Complex>{
 private:
  // Wavenumber //
  double  k;
  Complex epsr;
  Complex mur;

  // Function Space & Domain //
  const FunctionSpace*  fspace;
  const GroupOfElement* goe;

  // Local Terms //
  Term<double>* localTerms;

 public:
  FormulationImpedance(const GroupOfElement& domain,
                       const FunctionSpace& fs,
                       double k, Complex espr, Complex mur);

  virtual ~FormulationImpedance(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationImpedance::FormulationImpedance
   @param domain A GroupOfElement for the domain
   @param fs A FunctionSpace for both the unknown and test field
   @param k The wavenumber
   @param epsr The (complex) relative dielectric permitivity
   @param mur The (complex) relative magnetic permeability

   Instantiates a new FormulationImpedance with the given parameters
   **

   @fn FormulationImpedance::~FormulationImpedance
   Deletes this FormulationImpedance
*/

#endif
