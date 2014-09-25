Group{
  /*
  GammaS = Region[5]; // Source
  GammaN = Region[6]; // Neumann
  Omega  = Region[7]; // Omega
  */

  GammaS = Region[1000]; // Source
  GammaN = Region[4000]; // Neumann
  Omega  = Region[100];  // Omega

}

Function{
  k = 5;
  I[] = Complex[0, 1];

  theta_inc = 0;
  XYZdotTheta[] = X[] * Cos[theta_inc] + Y[] * Sin[theta_inc];
  //F[] = Complex[Cos[k*XYZdotTheta[]], Sin[k*XYZdotTheta[]]];

  //F[] = Exp[-((Y[] * 4.2) * (Y[] * 4.2) + (Z[] * 4.2) * (Z[] * 4.2))];
  //F[] = Fabs[Y[]];

  F[] = 1;
}

Constraint{
  { Name Dirichlet ;
    Case {
      { Region GammaS ; Value F[] ; }
    }
  }
}

Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }

  { Name JSur ;
    Case {
      { Region All ; Jacobian Sur ; }
    }
  }
}

Integration {
  { Name I1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Line ; NumberOfPoints 2 ; }
          { GeoElement Triangle ; NumberOfPoints 3 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 4 ; }
        }
      }
    }
  }
}

FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name se; NameOfCoef ee; Function BF_Node; Support Region[{Omega,GammaS,GammaN}] ; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef ee ; EntityType NodesOf ; NameOfConstraint Dirichlet ; }
    }
  }
  /*
  { Name HgradLagrange; Type Form0;
    BasisFunction {
      { Name le; NameOfCoef le; Function BF_Node; Support Region[{GammaS}] ; Entity NodesOf[All]; }
    }
  }
  */
}


Formulation {
  { Name FreeSpace; Type FemEquation;
    Quantity {
      { Name e; Type Local; NameOfSpace Hgrad; }
      //{ Name l; Type Local; NameOfSpace HgradLagrange; }
    }
    Equation {
      // Helmholtz
      Galerkin { [ Dof{d e} , {d e} ];
                 In Omega; Integration I1; Jacobian JVol;  }

      Galerkin { [ -k^2 * Dof{e} , {e} ];
                 In Omega; Integration I1; Jacobian JVol;  }

      // Somerfeld
      Galerkin { [ -1 * I[] * k * Dof{e} , {e} ];
                 In GammaN; Integration I1; Jacobian JSur;  }

      /*
      // Lagrange
      Galerkin { [ Dof{l}, {e} ];
                 In GammaS; Integration I1; Jacobian JSur;  }

      Galerkin { [ Dof{e}, {l} ];
                 In GammaS; Integration I1; Jacobian JSur;  }

      Galerkin { [ -F[], {l} ];
                 In GammaS; Integration I1; Jacobian JSur;  }
      */
    }
  }
}


Resolution {
  { Name FreeSpace ;
    System {
      { Name A ; NameOfFormulation FreeSpace ; Type Complex; }
    }
    Operation {
      Generate[A] ; Solve[A] ; SaveSolution[A] ;
    }
  }
}


PostProcessing {
  { Name FreeSpace ; NameOfFormulation FreeSpace ;
    Quantity {
      { Name e ;
        Value { Local { [ {e} ] ; In Omega; Jacobian JVol ; } } }
      /*
      { Name l ;
        Value { Local { [ {l} ] ; In GammaS; Jacobian JSur ; } } }
      */
    }
  }
}


PostOperation {
  { Name FreeSpace ; NameOfPostProcessing FreeSpace;
    Operation {
      Print[ e, OnElementsOf Omega, File "free.pos"] ;
      //Print[ l, OnElementsOf GammaS, File "lagfree.pos"] ;
    }
  }
}
