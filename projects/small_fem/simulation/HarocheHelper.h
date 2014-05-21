#ifndef _HAROCHEHELPER_H_
#define _HAROCHEHELPER_H_

#include "SmallFem.h"
#include "fullMatrix.h"

/**
   @class Haroche Helper
   @brief Helping stuff for simulation/Haroche

   Helping stuff for simulation/Haroche
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

  class XYZ{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class XZ{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class YZ{
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

  class Z{
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

 public:
  class Air{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);
    static Complex sz(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class XYZ{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);
    static Complex sz(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class XZ{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);
    static Complex sz(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class YZ{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);
    static Complex sz(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class XY{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);
    static Complex sz(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class X{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);
    static Complex sz(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class Y{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);
    static Complex sz(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };

  class Z{
  public:
    static Complex sx(fullVector<double>& xyz);
    static Complex sy(fullVector<double>& xyz);
    static Complex sz(fullVector<double>& xyz);

    static Complex Lxx(fullVector<double>& xyz);
    static Complex Lyy(fullVector<double>& xyz);
    static Complex Lzz(fullVector<double>& xyz);
  };
};

#endif
