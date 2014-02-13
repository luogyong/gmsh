#ifndef _TERM_H_
#define _TERM_H_

#include <omp.h>
#include "fullMatrix.h"

/**
   @interface Term
   @brief Interface for finite element terms

   Interface for finite element terms
 */

template<typename scalar>
class Term{
 protected:
  size_t nFunction;
  size_t nOrientation;
  const std::vector<size_t>* orientationStat;

  fullMatrix<scalar>** aM;

  mutable bool*   once;
  mutable size_t* lastId;
  mutable size_t* lastI;
  mutable size_t* lastCtr;

 public:
  virtual ~Term(void);

  scalar getTerm(size_t dofI, size_t dofJ, size_t elementId) const;

 private:
  scalar getTermOutCache(size_t dofI, size_t dofJ, size_t elementId,
                         size_t threadId) const;

 protected:
  Term(void);

  void   allocA(size_t nFunction);
  void computeA(fullMatrix<scalar>**& bM, fullMatrix<scalar>**& cM);
  void    clean(fullMatrix<scalar>**& bM, fullMatrix<scalar>**& cM);
};

/**
   @fn Term::~Term
   Deletes this Term
   **

   @fn Term::getTerm
   @param dofI A FunctionSpace function index
   @param dofJ A FunctionSpace function index
   @param elementId The ID of an Element
   @return Returns the finite element term associated to the given values
 */

/////////////////////
// Inline Function //
/////////////////////

template<typename scalar>
inline scalar Term<scalar>::
getTerm(size_t dofI, size_t dofJ, size_t elementId) const{

  const size_t threadId = omp_get_thread_num();

  if(!once[threadId] || elementId != lastId[threadId])
    // If Out Of Cache --> Fetch
    return getTermOutCache(dofI, dofJ, elementId, threadId);

  else
    // Else, rock baby yeah !
    return (*aM[lastI[threadId]])(lastCtr[threadId], dofI * nFunction + dofJ);
}

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermInclusion.h"

#endif
