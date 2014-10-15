#ifndef _FORMULATIONOSRCVECTORFOUR_H_
#define _FORMULATIONOSRCVECTORFOUR_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"
#include "TermCurlCurl.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorFour
   @brief Helping class for FormulationOSRCVector <d(phi), d(phi)> + <phi, phi>

   Helping class for FormulationOSRCVector <curl(phi), curl(phi)> + <phi, phi>

   FormulationOSRCVector is a friend of FormulationOSRCVectorFour
 */

class FormulationOSRCVectorFour: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Wavenumbers //
  Complex minusOneOverKEpsSquare;
  Complex Bi;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>* localGG;
  const TermCurlCurl<double>* localCC;

 private:
  FormulationOSRCVectorFour(void);
  FormulationOSRCVectorFour(const GroupOfElement& domain,
                            const FunctionSpace& field,
                            Complex kEps,
                            Complex Bi,
                            const TermGradGrad<double>& localGG,
                            const TermCurlCurl<double>& localCC);

 public:
  virtual ~FormulationOSRCVectorFour(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorFour::~FormulationOSRCVectorFour
   Deletes this FormulationOSRCVectorFour
*/

#endif
