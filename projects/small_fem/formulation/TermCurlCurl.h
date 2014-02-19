#ifndef _TERMCURLCURL_H_
#define _TERMCURLCURL_H_

#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermCurlCurl
   @brief Term of the Curl Curl type

   Term of the Curl Curl type
 */

template<typename scalar>
class TermCurlCurl: public Term<scalar>{
 private:
  typedef const fullMatrix<double>& (Basis::*bFunction)(size_t s)const;

 public:
  TermCurlCurl(const GroupOfJacobian& goj,
               const Basis& basis,
               const fullVector<double>& integrationWeights);

  virtual ~TermCurlCurl(void);

 private:
  void computeC(const Basis& basis,
                const bFunction& getFunction,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                size_t nG,
                fullMatrix<scalar>**& bM);
};

/**
   @fn TermCurlCurl::TermCurlCurl
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param integrationWeights A set of integration weights

   Instanciates a new Curl-Curl Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The Basis function must be pre-evaluated at the integration points
   (corresponding to the given integration weights)
   **

   @fn TermCurlCurl::~TermCurlCurl
   Deletes this TermCurlCurl
*/


//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermCurlCurlInclusion.h"

#endif
