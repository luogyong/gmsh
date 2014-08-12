#ifndef _TERMDUMMY_H_
#define _TERMDUMMY_H_

#include "Term.h"

/**
   @class TermDummy
   @brief Term that does nothing

   Term that does nothing
 */

template<typename scalar>
class TermDummy: public Term<scalar>{
 private:
  std::vector<size_t> dummyStat;

 public:
  TermDummy(void);

  virtual ~TermDummy(void);
};

/**
   @fn TermDummy::TermDummy
   Instanciates a new Dummy Term
   **

   @fn TermDummy::~TermDummy
   Deletes this TermDummy
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermDummyInclusion.h"

#endif
