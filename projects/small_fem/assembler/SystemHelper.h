#ifndef _SYSTEMHELPER_H_
#define _SYSTEMHELPER_H_

/**
   @class SystemHelper
   @brief A bunch of helping class methods for linear systems

   A bunch of helping class methods for linear systems
*/

#include "FormulationProjection.h"
#include "SystemAbstract.h"
#include "GroupOfElement.h"
#include "FunctionSpace.h"

template<typename scalar>
class SystemHelper{
 public:
   SystemHelper(void);
  ~SystemHelper(void);

  static void dirichlet(SystemAbstract<scalar>& sys,
                        const FunctionSpace& fs,
                        const GroupOfElement& domain,
                        scalar (*f)(fullVector<double>& xyz));

  static void dirichlet(SystemAbstract<scalar>& sys,
                        const FunctionSpace& fs,
                        const GroupOfElement& domain,
                        fullVector<scalar> (*f)(fullVector<double>& xyz));
 private:
  static void dirichlet(SystemAbstract<scalar>& sys,
                        const FunctionSpace& fs,
                        const GroupOfElement& domain,
                        FormulationProjection<scalar>& formulation);
};

/**
   @fn SystemHelper::SystemHelper
   Instantiates a new SystemHelper (unneeded since it has only class methods)
   **

   @fn SystemHelper::~SystemHelper
   Deletes this SystemHelper
   **

   @fn SystemHelper::dirichlet(SystemAbstract<scalar>&, const FunctionSpace&, const GroupOfElement&, scalar (*f)(fullVector<double>& xyz))
   @param sys A SystemAbstract
   @param fs A FunctionSpace
   @param domain A GroupOfElement
   @param f A scalar function

   Imposes on the given SystemAbstract a scalar dirichlet condition with:
   @li The Dof%s of the given FunctionSpace (restricted to the given domain)
   @li The given scalar function
   **

   @fn SystemHelper::dirichlet(SystemAbstract<scalar>&, const FunctionSpace&, const GroupOfElement&, fullVector<scalar> (*f)(fullVector<double>& xyz))
   @param sys A SystemAbstract
   @param fs A FunctionSpace
   @param domain A GroupOfElement
   @param f A vectorial function

   Imposes on the given SystemAbstract a vectorial dirichlet condition with:
   @li The Dof%s of the given FunctionSpace (restricted to the given domain)
   @li The given vectorial function
 */

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SystemHelperInclusion.h"

#endif
