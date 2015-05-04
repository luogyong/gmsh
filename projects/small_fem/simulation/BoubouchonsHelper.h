#ifndef _BOUBOUCHONSHELPER_H_
#define _BOUBOUCHONSHELPER_H_

#include "SmallFem.h"
#include "fullMatrix.h"

/**
   @class Boubouchons Helper
   @brief Helping stuff for simulation/Boubouchons

   Helping stuff for simulation/Boubouchons
 */

// Material //
class Material{
 public:
  class XYZ{
  public:
    static void Eps(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class XZ{
  public:
    static void Eps(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class YZ{
  public:
    static void Eps(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class XY{
  public:
    static void Eps(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class X{
  public:
    static void Eps(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class Y{
  public:
    static void Eps(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class Z{
  public:
    static void Eps(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };
};

// PML //
class PML{
 private:
  static const Complex a;
  static const Complex one;

 private:
  static Complex dampingX(fullVector<double>& xyz);
  static Complex dampingY(fullVector<double>& xyz);
  static Complex dampingZ(fullVector<double>& xyz);

 public:
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
