#ifndef _DDMCONTEXT_H_
#define _DDMCONTEXT_H_

#include <vector>
#include <string>
#include <map>

#include "Dof.h"
#include "FunctionSpace.h"
#include "GroupOfElement.h"
#include "System.h"
#include "SmallFem.h"

/**
   @interface DDMContext
   @brief Common interface for all context in DDM

   This class is the common interface for all context in DDM
   A DDMContext can be:
   @li DDMContextEMDA;
   @li DDMContextOO2;
   @li DDMContextOSRC.

   A DDMContext handles the variables common to all DDM Formulation%s,
   that is:
   @li the DDM border;
   @li the Dirichlet border;
   @li the field FunctionSpace;
   @li the DDM FunctionSpace;
   @li the volume problem System;
   @li the DDM Dof values.

   Child classes must implement method specific variables.
 */

class DDMContext{
 protected:
  const System<Complex>* system;
  const GroupOfElement*  domain;
  const FunctionSpace*   fSpace;
  const FunctionSpace*   fSpaceG;

  std::vector<const GroupOfElement*> dirichlet;

  std::map<Dof, Complex>* ddm;

 public:
  DDMContext(void);
  virtual ~DDMContext(void);

  void setSystem(const System<Complex>& system);
  void setDDMDofs(std::map<Dof, Complex>& ddm);

  std::map<Dof, Complex>& getDDMDofs(void);

  const System<Complex>& getSystem(void)         const;
  const FunctionSpace&   getFunctionSpace(void)  const;
  const FunctionSpace&   getFunctionSpaceG(void) const;
  const GroupOfElement&  getDomain(void)         const;

  const std::vector<const GroupOfElement*>& getDirichletDomain(void) const;
};

/**
   @fn DDMContext::~DDMContext
   Deletes this DDMContext
   **

   @fn DDMContext::setSytem
   @param system A System

   Sets the given System
   as the System solving the volume problem of this DDM problem
   **

   @fn DDMContext::setDDMDofs
   @param ddm A Dof -- Value map

   Sets the given map as the DDM Dof values of this DDM problem
   **

   @fn DDMContext::getSystem
   @return Returns the System solving the volume problem of this DDMContext
   **

   @fn DDMContext::getDDMDofs
   @return Returns the Dof -- Value map (at DDM border) of this DDMContext
   **

   @fn DDMContext::getFunctionSpace
   @return Returns the FunctionSpace used for the unknown field
   of this DDMContext
   **

   @fn DDMContext::getFunctionSpaceG
   @return Returns the FunctionSpace used for the DDM field
   of this DDMContext
   **

   @fn DDMContext::getDomain
   @return Returns the domain defining the DDM border of this DDMContext
   **

   @fn DDMContext::getDirichletDomain
   @return Returns the domain defining the support of the Dirichlet conditions
   of this DDMContext
 */

//////////////////////
// Inline Functions //
//////////////////////

inline void DDMContext::setSystem(const System<Complex>& system){
  this->system = &system;
}

inline void DDMContext::setDDMDofs(std::map<Dof, Complex>& ddm){
  this->ddm = &ddm;
}

inline const System<Complex>& DDMContext::getSystem(void) const{
  return *system;
}

inline std::map<Dof, Complex>& DDMContext::getDDMDofs(void){
  return *ddm;
}

inline const FunctionSpace& DDMContext::getFunctionSpace(void) const{
  return *fSpace;
}

inline const FunctionSpace& DDMContext::getFunctionSpaceG(void) const{
  return *fSpaceG;
}

inline const GroupOfElement& DDMContext::getDomain(void) const{
  return *domain;
}

inline const std::vector<const GroupOfElement*>&
DDMContext::getDirichletDomain(void) const{
  return dirichlet;
}

#endif
