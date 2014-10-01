Group{

  GammaS = Region[5]; // Source
  GammaN = Region[6]; // Neumann
  Omega  = Region[7]; // Omega

  /*
  GammaS = Region[1000]; // Source
  GammaN = Region[4000]; // Neumann
  Omega  = Region[100];  // Omega
  */
}

Function{
  k = 5;
  I[] = Complex[0, 1];
  N[] = Normal[];

  //theta_inc = 0;
  //XYZdotTheta[] = X[] * Cos[theta_inc] + Y[] * Sin[theta_inc];
  //F[] = Complex[Cos[k*XYZdotTheta[]], Sin[k*XYZdotTheta[]]];

  //F[] = Exp[-((Y[] * 4.2) * (Y[] * 4.2) + (Z[] * 4.2) * (Z[] * 4.2))];
  //F[] = Fabs[Y[]];

  F[] = Vector[0, 1, 0];
}

Constraint{
  { Name Dirichlet ;
    Case {
      { Region GammaS ;
        Type AssignFromResolution ;
        NameOfResolution Projection ; }
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
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name se; NameOfCoef ee; Function BF_Edge;
        Support Region[{Omega,GammaS,GammaN}] ; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef ee ; EntityType EdgesOf ; NameOfConstraint Dirichlet ; }
    }
  }
}


Formulation {
  { Name FreeSpace; Type FemEquation;
    Quantity {
      { Name e; Type Local; NameOfSpace Hcurl; }
    }
    Equation {
      // Helmholtz
      Galerkin { [ Dof{d e} , {d e} ];
                 In Omega; Integration I1; Jacobian JVol;  }

      Galerkin { [ -k^2 * Dof{e} , {e} ];
                 In Omega; Integration I1; Jacobian JVol;  }

      // Silver-Muller
      // Should be equivalent to [ I[] * k * Dof{e}, {e} ] since we are on Gamma
      Galerkin { [ I[] * k * ( (N[]) /\ (N[] /\ Dof{e}) ) , {e} ];
                 In GammaN; Integration I1; Jacobian JSur; }
    }
  }

  { Name Projection;
    Quantity {
      { Name e; Type Local; NameOfSpace Hcurl; }
    }
    Equation {
      Galerkin { [ Dof{e} , {e} ];
                 In GammaS; Integration I1; Jacobian JSur; }
      Galerkin { [ F[] , {e} ];
                 In GammaS; Integration I1; Jacobian JSur; }
    }
  }
}


Resolution {
  { Name FreeSpace ;
    System {
      { Name A ; NameOfFormulation FreeSpace ;
        Type Complex; }
    }
    Operation {
      Generate[A] ; Solve[A] ; SaveSolution[A] ;
    }
  }

  { Name Projection;
    System {
      { Name B; NameOfFormulation Projection; DestinationSystem A;
        Type Complex; }
    }
    Operation {
      Generate[B]; Solve[B]; TransferSolution[B];
    }
  }
}


PostProcessing {
  { Name FreeSpace ; NameOfFormulation FreeSpace ;
    Quantity {
      { Name e ;
        Value { Local { [ {e} ] ; In Omega; Jacobian JVol ; } } }
    }
  }
}


PostOperation {
  { Name FreeSpace ; NameOfPostProcessing FreeSpace;
    Operation {
      Print[ e, OnElementsOf Omega, File "free.pos"] ;
    }
  }
}
