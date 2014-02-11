#ifndef _TERMPROJECTIONGRAD_H_
#define _TERMPROJECTIONGRAD_H_

#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionGrad
   @brief Term of the Field (in physical space) Grad (in reference space) type

   Term of the Field (in physical space) Grad (in reference space) type
 */

template<typename scalar>
class TermProjectionGrad: public Term<scalar>{
 private:
  typedef const fullMatrix<double>& (Basis::*bFunction)(size_t s)const;

 public:
  TermProjectionGrad(const GroupOfJacobian& goj,
                     const Basis& basis,
                     const fullVector<double>& integrationWeights,
                     const fullMatrix<double>& integrationPoints,
                     fullVector<scalar> (*f)(fullVector<double>& xyz));

  virtual ~TermProjectionGrad(void);


 private:
  void computeC(const Basis& basis,
                const bFunction& getFunction,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                const Basis& basis,
                const fullMatrix<double>& gC,
                fullVector<scalar> (*f)(fullVector<double>& xyz),
                fullMatrix<scalar>**& bM);
};

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermProjectionGradInclusion.h"

#endif
