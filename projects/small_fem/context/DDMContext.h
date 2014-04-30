#ifndef _DDMCONTEXT_H_
#define _DDMCONTEXT_H_

#include <vector>
#include <string>
#include <map>

#include "Dof.h"
#include "FunctionSpace.h"
#include "FunctionSpaceScalar.h"
#include "GroupOfElement.h"
#include "System.h"
#include "SmallFem.h"

/**
   @class DDMContext
   @brief Context for DDM

   This class is a DDM context.
   A context can be:
   @li EMDA
   @li OO2
   @li OSRC
 */

class DDMContext{
 private:
  friend class FormulationEMDA;
  friend class FormulationOO2;
  friend class FormulationOSRC;

  friend class FormulationUpdateEMDA;
  friend class FormulationUpdateOO2;
  friend class FormulationUpdateOSRC;

 private:
  static const std::string typeNULL;
  static const std::string typeEMDA;
  static const std::string typeOO2;
  static const std::string typeOSRC;

 private:
  std::string typeDDM;

  const System<Complex>*                         system;
  const GroupOfElement*                          domain;
  const FunctionSpaceScalar*                     fSpace;
  const std::vector<const FunctionSpaceScalar*>* phi;

  const std::map<Dof, Complex>* ddm;

  double  k;
  double  EMDA_Chi;
  Complex OO2_A;
  Complex OO2_B;
  Complex OSRC_keps;
  int     OSRC_NPade;

 public:
   DDMContext(void);
  ~DDMContext(void);

  void setToEMDA(const GroupOfElement& domain,
                 const FunctionSpaceScalar& fSpace,
                 double k, double chil);

  void setToOO2(const GroupOfElement& domain,
                const FunctionSpaceScalar& fSpace,
                Complex a, Complex b);

  void setToOSRC(const GroupOfElement& domain,
                 const FunctionSpaceScalar& fSpace,
                 const std::vector<const FunctionSpaceScalar*>& phi,
                 double k, Complex keps, int NPade);

  void setSystem(const System<Complex>& system);
  void setDDMDofs(const std::map<Dof, Complex>& ddm);

  const System<Complex>&        getSystem(void)        const;
  const std::map<Dof, Complex>& getDDMDofs(void)       const;
  const FunctionSpace&          getFunctionSpace(void) const;
  const GroupOfElement&         getDomain(void)        const;
  std::string                   getType(void)          const;
};

//////////////////////
// Inline Functions //
//////////////////////

inline void DDMContext::setSystem(const System<Complex>& system){
  this->system = &system;
}

inline void DDMContext::setDDMDofs(const std::map<Dof, Complex>& ddm){
  this->ddm = &ddm;
}

inline const System<Complex>& DDMContext::getSystem(void) const{
  return *system;
}

inline const std::map<Dof, Complex>& DDMContext::getDDMDofs(void) const{
  return *ddm;
}

inline const FunctionSpace& DDMContext::getFunctionSpace(void) const{
  return *fSpace;
}

inline const GroupOfElement& DDMContext::getDomain(void) const{
  return *domain;
}

inline std::string DDMContext::getType(void) const{
  return typeDDM;
}


#endif
