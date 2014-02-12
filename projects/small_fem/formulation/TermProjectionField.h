#ifndef _TERMPROJECTIONFIELD_H_
#define _TERMPROJECTIONFIELD_H_

#include <map>

#include "FunctionSpaceScalar.h"
#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionField
   @brief Term of the Field (in physical space) Field (in reference space) type

   Term of the Field (in physical space) Field (in reference space) type
 */

template<typename scalar>
class TermProjectionField: public Term<scalar>{
 private:
  // Type def //
  typedef scalar
    (TermProjectionField<scalar>::*Eval)(const MElement& element,
                                         const fullVector<double>& xyz) const;
 private:
  // Function for fContainer //
  scalar (*f)(fullVector<double>& xyz);

  // Data for interpolation //
  const FunctionSpaceScalar*   fsScalar;
  const std::map<Dof, scalar>* dofValue;

 public:
  TermProjectionField(const GroupOfJacobian& goj,
                      const Basis& basis,
                      const fullVector<double>& integrationWeights,
                      const fullMatrix<double>& integrationPoints,
                      scalar (*f)(fullVector<double>& xyz));

  TermProjectionField(const GroupOfJacobian& goj,
                      const Basis& basis,
                      const fullVector<double>& integrationWeights,
                      const fullMatrix<double>& integrationPoints,
                      const FunctionSpaceScalar& fs,
                      const std::map<Dof, scalar>& dof);

  virtual ~TermProjectionField(void);

 private:
  void init(const GroupOfJacobian& goj,
            const Basis& basis,
            const fullVector<double>& integrationWeights,
            const fullMatrix<double>& integrationPoints,
            const Eval& evaluator);

  void computeC(const Basis& basis,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                const Basis& basis,
                const fullMatrix<double>& gC,
                const Eval& evaluator,
                fullMatrix<scalar>**& bM);

  scalar  fContainer(const MElement& elem, const fullVector<double>& xyz) const;
  scalar interpolate(const MElement& elem, const fullVector<double>& xyz) const;
};

/////////////////////
// Inline Function //
/////////////////////
template<typename scalar>
inline scalar TermProjectionField<scalar>::
fContainer(const MElement& element,
           const fullVector<double>& xyz) const{
  return f(const_cast<fullVector<double>&>(xyz));
}

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermProjectionFieldInclusion.h"

#endif
