#ifndef _DIFFRACTIONHELPER_H_
#define _DIFFRACTIONHELPER_H_

#include "SmallFem.h"
#include "fullMatrix.h"

/**
   @class Diffraction Helper
   @brief Helping stuff for simulation/Diffraction

   Helping stuff for simulation/Diffraction
 */

// Math constant //
class Math{
 public:
  static const double Pi;
  static const double nm;
};

// Wave //
class Wave{
 public:
  static const double lambda0;
  static const double theta0;
  static const double phi0;
  static const double psi0;
};

// Material //
class Material{
 public:
  static const double  mu0;
  static const double  epsilon0;
  static const Complex epsIn;
  static const Complex epsOut;

  class In{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class Out{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class XYZ{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class XZ{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class YZ{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class XY{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class X{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class Y{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };

  class Z{
  public:
    static void Epsilon(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Epsilon1(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Mu(fullVector<double>& xyz, fullMatrix<Complex>& T);
    static void Nu(fullVector<double>& xyz, fullMatrix<Complex>& T);
  };
};

// Constant //
class Constant{
 public:
  static const double cel;
  static const double Freq;
  static const double omega0;
  static const double k0;

  static const double Ae;
  static const double alpha0;
  static const double beta0;
  static const double gamma0;

  static const double Ex0;
  static const double Ey0;
  static const double Ez0;
};

// Signal //
class Signal{
 public:
  static Complex Prop(fullVector<double>& xyz);
  static fullVector<Complex> Einc(fullVector<double>& xyz);

  class In{
  public:
    static fullVector<Complex> source(fullVector<double>& xyz);
  };

  class Out{
  public:
    static fullVector<Complex> source(fullVector<double>& xyz);
  };

  class PML{
  public:
    static fullVector<Complex> source(fullVector<double>& xyz);
  };
};

// PML //
class PML{
 private:
  static const Complex a;
  static const Complex one;

 public:
  class Volume{
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
