Include "waveguide.dat";

Group{
  Omega    = Region[1];
  GammaD0  = Region[2];
  GammaD   = Region[3];
  GammaInf = Region[4];
  Domain   = Region[{Omega, GammaD0, GammaD, GammaInf}];
}

Function{
  I[] = Complex[0, 1];

  ky   = MODE_M * Pi / LY;
  kz   = MODE_N * Pi / LZ;
  kc   = Sqrt[ky^2 + kz^2];
  kx[] = (-kc^2 + k^2 >=0 ? Sqrt[-kc^2 + k^2] : -I[] * Sqrt[kc^2 - k^2]);

  If(POLAR == 0)
    // TE Mode
    src[]  = Vector[0,
                    -Cos[ky * Y[]] * Sin[kz * Z[]],
                    +Sin[ky * Y[]] * Cos[kz * Z[]]];
    kInf[] = kx[];
  EndIf

  If(POLAR == 1)
    // TM Mode
    src[]  = Vector[                         Sin[ky * Y[]] * Sin[kz * Z[]],
                    I[] * kx[] * ky / kc^2 * Cos[ky * Y[]] * Sin[kz * Z[]],
                    I[] * kx[] * kz / kc^2 * Sin[ky * Y[]] * Cos[kz * Z[]]];
    kInf[] = k^2 / kx[];
  EndIf
}

Jacobian{
  { Name JVol; Case{ { Region All; Jacobian Vol; } } }
  { Name JSur; Case{ { Region All; Jacobian Sur; } } }
}

Integration{
  { Name I1;
    Case{
      { Type Gauss;
        Case{
          { GeoElement Line       ; NumberOfPoints 2; }
          { GeoElement Triangle   ; NumberOfPoints 3; }
          { GeoElement Tetrahedron; NumberOfPoints 4; }
        }
      }
    }
  }
}

Constraint{
  { Name Dirichlet; Case{
      { Region GammaD0; Type Assign;               Value 0.;              }
      { Region GammaD;  Type AssignFromResolution; NameOfResolution Proj; }
    }
  }
}

FunctionSpace{
  { Name HCurl; Type Form1;
    BasisFunction{
      { Name se1; NameOfCoef ee1;
        Function BF_Edge_1E; Support Region[Domain]; Entity EdgesOf[All]; }
      { Name se2; NameOfCoef ee2;
        Function BF_Edge_2E; Support Region[Domain]; Entity EdgesOf[All]; }
    }

    Constraint{
      { NameOfCoef ee1; EntityType EdgesOf; NameOfConstraint Dirichlet; }
      { NameOfCoef ee2; EntityType EdgesOf; NameOfConstraint Dirichlet; }
    }
  }
}

Formulation{
  { Name Maxwell;
    Quantity{ { Name e; Type Local; NameOfSpace HCurl; } }
    Equation{
      Galerkin{ [Dof{d e}, {d e}];    In Omega; Integration I1; Jacobian JVol; }
      Galerkin{ [-k^2 * Dof{e}, {e}]; In Omega; Integration I1; Jacobian JVol; }
      Galerkin{ [-I[] * kInf[] * Dof{e}, {e}];
        In GammaInf; Integration I1; Jacobian JSur; }
    }
  }

  { Name Proj;
    Quantity{ { Name e; Type Local; NameOfSpace HCurl; } }
    Equation{
      Galerkin{ [Dof{e}, {e}]; In GammaD; Integration I1; Jacobian JSur; }
      Galerkin{ [-src[], {e}]; In GammaD; Integration I1; Jacobian JSur; }
    }
  }
}

Resolution{
  { Name Maxwell;
    System{ { Name A; NameOfFormulation Maxwell; Type Complex; } }
    Operation{ Generate[A]; Solve[A]; }
  }

  { Name Proj;
    System{ { Name B; NameOfFormulation Proj;    Type Complex;
              DestinationSystem A; } }
    Operation{ Generate[B]; Solve[B]; TransferSolution[B]; }
  }
}

PostProcessing{
  { Name Maxwell; NameOfFormulation Maxwell;
    Quantity{ { Name e; Value{ Local{ [{e}]; In Omega; Jacobian JVol; } } } }
  }
}

PostOperation{
  { Name Maxwell ; NameOfPostProcessing Maxwell;
    Operation{ Print[e, OnElementsOf Omega, File "waveguide.pos"]; }
  }
}
