#ifndef _FORMULATIONOSRCVECTORONE_H_
#define _FORMULATIONOSRCVECTORONE_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermProjectionGrad.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorOne
   @brief Helping class for FormulationOSRCVector <r, r> - <g, r>

   Helping class for FormulationOSRCVector <r, r> -<g, r>

   FormulationOSRCVector is a friend of FormulationOSRCVectorOne
 */

class FormulationOSRCVectorOne: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Wavenumber //
  Complex jOverK;
  Complex C0;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>*        localLHS;
  const TermProjectionGrad<Complex>* localRHS;

 private:
  FormulationOSRCVectorOne(void);
  FormulationOSRCVectorOne(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           double  k,
                           Complex C0,
                           const TermGradGrad<double>& localLHS,
                           const TermProjectionGrad<Complex>& localRHS);

 public:
  virtual ~FormulationOSRCVectorOne(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;

 private:
  void update(TermProjectionGrad<Complex>& localRHS);
};

/**
   @fn FormulationOSRCVectorOne::~FormulationOSRCVectorOne
   Deletes this FormulationOSRCVectorOne
*/

#endif
