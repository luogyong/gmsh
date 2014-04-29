#ifndef _FORMULATIONEMDA_H_
#define _FORMULATIONEMDA_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "FormulationBlock.h"

/**
   @class FormulationEMDA
   @brief EMDA Formulation for DDM

   EMDA Formulation for DDM
 */

class FormulationEMDA: public FormulationBlock<Complex>{
 private:
  // Wavenumber & Chi //
  double k;
  double chi;

  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField<double>*       localLHS;
  TermProjectionField<Complex>* localRHS;

 public:
  FormulationEMDA(const GroupOfElement& domain,
                  const FunctionSpaceScalar& fs,
                  double k,
                  double chi,
                  const std::map<Dof, Complex>& ddmDof);

  virtual ~FormulationEMDA(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationEMDA::FormulationEMDA
   @param domain A GroupOfElement for the domain
   @param fs A FunctionSpace for both unknown and test fields
   @param k A real number
   @param chi A real number
   @param ddmDof A map with the DDM Dof%s and their associated values

   Instantiates a new FormulationEMDA with wavenumber k and real shift chi
   The DDM Dof%s are given by ddmDof.
   **

   @fn FormulationEMDA::~FormulationEMDA
   Deletes this FormulationEMDA
*/

#endif
