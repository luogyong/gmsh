#ifndef _FORMULATIONOSRCVECTORTWO_H_
#define _FORMULATIONOSRCVECTORTWO_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermCurlCurl.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorTwo
   @brief Helping class for FormulationOSRCVector <e, r> + <curl(e), curl(r)>

   Helping class for FormulationOSRCVector  <e, r> + <curl(e), curl(r)>

   FormulationOSRCVector is a friend of FormulationOSRCVectorTwo
 */

class FormulationOSRCVectorTwo: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Wavenumber //
  Complex minusOneOverKEpsSquare;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>* localGG;
  const TermCurlCurl<double>* localCC;

 private:
  FormulationOSRCVectorTwo(void);
  FormulationOSRCVectorTwo(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const FunctionSpace& test,
                           Complex kEps,
                           const TermGradGrad<double>& localGG,
                           const TermCurlCurl<double>& localCC);

 public:
  virtual ~FormulationOSRCVectorTwo(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorTwo::~FormulationOSRCVectorTwo
   Deletes this FormulationOSRCVectorTwo
*/

#endif
