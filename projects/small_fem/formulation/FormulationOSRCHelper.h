#ifndef _FORMULATIONOSRCHELPER_H_
#define _FORMULATIONOSRCHELPER_H_

#include "SmallFem.h"

/**
   @class FormulationOSRCHelper
   @brief Helping function for OSRC DDM Formulation

   Helping function for OSRC DDM Formulation (scalar and vectorial)
*/

class FormulationOSRCHelper{
 public:
   FormulationOSRCHelper(void);
  ~FormulationOSRCHelper(void);

 private:
  static double pade_aj(int j, int N);
  static double pade_bj(int j, int N);

 public:
  static Complex padeC0(int N, double theta);
  static Complex padeR0(int N, double theta);
  static Complex padeA(int j, int N, double theta);
  static Complex padeB(int j, int N, double theta);
};

/**
   @fn FormulationOSRCHelper::FormulationOSRCHelper
   Instantiates a new FormulationOSRCHelper.
   This is not required since FormulationOSRCHelper is purely static
   **

   @fn FormulationOSRCHelper::~FormulationOSRCHelper
   Deletes this FormulationOSRCHelper
   **

   @fn FormulationOSRCHelper::padeC0
   @param N A natural number
   @param theta A real number
   @return Returns the C0 term of the OSRC Pade approximation for N terms
   and complex rotation of theta
   **

   @fn FormulationOSRCHelper::padeR0
   @param N A natural number
   @param theta A real number
   @return Returns the R0 term of the OSRC Pade approximation for N terms
   and complex rotation of theta
   **

   @fn FormulationOSRCHelper::padeA
   @param j A natural number
   @param N A natural number
   @param theta A real number
   @return Returns the jth A term of the OSRC Pade approximation for N terms
   and complex rotation of theta
   **

   @fn FormulationOSRCHelper::padeB
   @param j A natural number
   @param N A natural number
   @param theta A real number
   @return Returns the jth B term of the OSRC Pade approximation for N terms
   and complex rotation of theta
*/

#endif
