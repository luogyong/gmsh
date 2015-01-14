#ifndef _TERMCURLCURL_H_
#define _TERMCURLCURL_H_

#include "GroupOfJacobian.h"
#include "Quadrature.h"
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
  typedef const fullMatrix<double>& (Basis::*BFunction)(size_t s)const;

 private:
  std::vector<fullMatrix<scalar> > TJac; //Pre-evaluated tensor * jacobian

 public:
  TermCurlCurl(const GroupOfJacobian& goj,
               const Basis& basis,
               const Quadrature& quadrature);

  TermCurlCurl(const GroupOfJacobian& goj,
               const Basis& field,
               const Basis& test,
               const Quadrature& quadrature);

  TermCurlCurl(const GroupOfJacobian& goj,
               const Basis& basis,
               const Quadrature& quadrature,
               void  (*f)(fullVector<double>& xyz, fullMatrix<scalar>& T));

  virtual ~TermCurlCurl(void);

 private:
  // Init
  void init(const GroupOfJacobian& goj,
            const Basis& field,
            const Basis& test,
            const Quadrature& quadrature);

  // Matrices
  void computeC(const Basis& field,
                const Basis& test,
                const BFunction& getFunctionField,
                const BFunction& getFunctionTest,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                size_t nG,
                fullMatrix<scalar>**& bM);

  // Tensor Pre Evaluator
  void preEvalDummy(const GroupOfJacobian& goj,
                    const Quadrature& quadrature);

  void preEvalT(const GroupOfJacobian& goj,
                const Quadrature& quadrature,
                void  (*f)(fullVector<double>& xyz, fullMatrix<scalar>& T));
};

/**
   @fn TermCurlCurl::TermCurlCurl(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature)
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule

   Instanciates a new Curl-Curl Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the integration points
   **

   @fn TermCurlCurl::TermCurlCurl(const GroupOfJacobian& goj,const Basis& field,const Basis& test,const Quadrature& quadrature)
   @param goj A GroupOfJacobian
   @param field A Basis for the unknwon field
   @param test A Basis for the test functions
   @param quadrature A Quadrature rule

   Instanciates a new Curl-Curl Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the field and test
   @li The given Quadrature is used to compute the Term
   @li The Basis functions must be pre-evaluated at the integration points
   **

   @fn TermCurlCurl::TermCurlCurl(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature,void(*f)(fullVector<double>& xyz,fullMatrix<scalar>& T))
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule
   @param f A multiplicative tensorial function

   Instanciates a new Curl-Curl Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the integration points
   @li The given function is a multiplicative tensor for the whole Term
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
