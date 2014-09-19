#ifndef _HAROCHEHELPER2D_H_
#define _HAROCHEHELPER2D_H_

#include "SmallFem.h"
#include "fullMatrix.h"

/**
   @class Haroche Helper 2D
   @brief Helping stuff for simulation/Haroche2D

   Helping stuff for simulation/Haroche2D
 */

// Material //
class Material{
 public:
  static const double oneOverCSquare;

  class Air{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class XY{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class X{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class Y{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };
};

// PML //
class PML{
 private:
  static const Complex a;
  static const Complex one;
  static const Complex sz;

 public:
  class Air{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class XY{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class X{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class Y{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };
};

#endif
