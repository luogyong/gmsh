#ifndef _FORMULATIONOSRCVECTORTEN_H_
#define _FORMULATIONOSRCVECTORTEN_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorTen
   @brief Helping class for FormulationOSRCVector <grad(rho), r>

   Helping class for FormulationOSRCVector <grad(rho), r>

   FormulationOSRCVector is a friend of FormulationOSRCVectorTen
 */

class FormulationOSRCVectorTen: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Wavenumbers //
  Complex jOverK;
  Complex plusOneOverKEpsSquare;
  Complex Ai;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>* localGG;

 private:
  FormulationOSRCVectorTen(void);
  FormulationOSRCVectorTen(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const FunctionSpace& test,
                           Complex kEps,
                           Complex Ai,
                           double  k,
                           const TermGradGrad<double>& localGG);

 public:
  virtual ~FormulationOSRCVectorTen(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorTen::~FormulationOSRCVectorTen
   Deletes this FormulationOSRCVectorTen
*/

#endif
