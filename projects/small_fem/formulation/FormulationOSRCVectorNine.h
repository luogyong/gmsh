#ifndef _FORMULATIONOSRCVECTORNINE_H_
#define _FORMULATIONOSRCVECTORNINE_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermCurlCurl.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorNine
   @brief Helping class for FormulationOSRCVector <curl(phi), curl(r)>

   Helping class for FormulationOSRCVector <curl(phi), curl(r)>

   FormulationOSRCVector is a friend of FormulationOSRCVectorNine
 */

class FormulationOSRCVectorNine: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Wavenumbers //
  Complex jOverK;
  Complex minusOneOverKEpsSquare;
  Complex Ai;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermCurlCurl<double>* localCC;

 private:
  FormulationOSRCVectorNine(void);
  FormulationOSRCVectorNine(const GroupOfElement& domain,
                            const FunctionSpace& field,
                            const FunctionSpace& test,
                            Complex kEps,
                            Complex Ai,
                            double  k,
                            const TermCurlCurl<double>& localCC);

 public:
  virtual ~FormulationOSRCVectorNine(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorNine::~FormulationOSRCVectorNine
   Deletes this FormulationOSRCVectorNine
*/

#endif
