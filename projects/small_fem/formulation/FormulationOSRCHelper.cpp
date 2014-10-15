#include <cmath>
#include "FormulationOSRCHelper.h"

using namespace std;

FormulationOSRCHelper::FormulationOSRCHelper(void){
}

FormulationOSRCHelper::~FormulationOSRCHelper(void){
}

double FormulationOSRCHelper::pade_aj(int j, int N){
  double tmp = sin((double)j * M_PI / (2. * N + 1.));

  return 2. / (2. * N + 1.) * tmp * tmp;
}

double FormulationOSRCHelper::pade_bj(int j, int N){
  double tmp = cos((double)j * M_PI / (2. * N + 1.));

  return tmp * tmp;
}

Complex FormulationOSRCHelper::padeC0(int N, double theta){
  Complex sum = Complex(1, 0);
  Complex one = Complex(1, 0);
  Complex z   = Complex(cos(-theta) - 1,  sin(-theta));

  for(int j = 1; j <= N; j++)
    sum += (z * pade_aj(j, N)) / (one + z * pade_bj(j, N));

  z = Complex(cos(theta / 2.), sin(theta / 2.));

  return sum * z;
}

Complex FormulationOSRCHelper::padeAj(int j, int N, double theta){
  Complex one = Complex(1, 0);
  Complex res;
  Complex z;

  z   = Complex(cos(-theta / 2.), sin(-theta / 2.));
  res = z * pade_aj(j, N);

  z   = Complex(cos(-theta) - 1., sin(-theta));
  res = res / ((one + z * pade_bj(j, N)) * (one + z * pade_bj(j, N)));

  return res;
}

Complex FormulationOSRCHelper::padeBj(int j, int N, double theta){
  Complex one = Complex(1, 0);
  Complex res;
  Complex z;

  z   = Complex(cos(-theta), sin(-theta));
  res = z * pade_bj(j, N);

  z   = Complex(cos(-theta) - 1., sin(-theta));
  res = res / (one + z * pade_bj(j, N));

  return res;
}
