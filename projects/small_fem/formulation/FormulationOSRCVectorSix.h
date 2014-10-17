#ifndef _FORMULATIONOSRCVECTORSIX_H_
#define _FORMULATIONOSRCVECTORSIX_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"
#include "TermCurlCurl.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorSix
   @brief Helping class for FormulationOSRCVector -<phi, phi> + <d(phi), d(phi)>

   Helping class for FormulationOSRCVector -<phi, phi> + <curl(phi), curl(phi)>

   FormulationOSRCVector is a friend of FormulationOSRCVectorSix
 */

class FormulationOSRCVectorSix: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Wavenumber & Pade //
  Complex BiOverKEpsSquare;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>* localGG;
  const TermCurlCurl<double>* localCC;

 private:
  FormulationOSRCVectorSix(void);
  FormulationOSRCVectorSix(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const FunctionSpace& test,
                           Complex kEps,
                           Complex Bi,
                           const TermGradGrad<double>& localGG,
                           const TermCurlCurl<double>& localCC);

 public:
  virtual ~FormulationOSRCVectorSix(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorSix::~FormulationOSRCVectorSix
   Deletes this FormulationOSRCVectorSix
*/

#endif
