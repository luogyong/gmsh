#ifndef _TERMFIELDFIELD_H_
#define _TERMFIELDFIELD_H_

#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermFieldField
   @brief A Term of the Field Field type

   A Term of the Field Field type
 */

class TermFieldField: public Term<double>{
 public:
  TermFieldField(const GroupOfJacobian& goj,
                 const Basis& basis,
                 const fullVector<double>& integrationWeights);

  virtual ~TermFieldField(void);

 private:
  void computeC(const Basis& basis,
                const fullVector<double>& gW,
                fullMatrix<double>**& cM);

  void computeB(const GroupOfJacobian& goj,
                size_t nG,
                fullMatrix<double>**& bM);
};

/**
   @fn TermFieldField::TermFieldField
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param integrationWeights A set of integration weights

   Instanciates a new Field-Field Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The Basis function must be pre-evaluated at the integration points
   (corresponding to the given integration weights)
   **

   @fn TermFieldField::~TermFieldField
   Deletes this TermFieldField
*/

#endif
