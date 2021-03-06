Include "boubouchons.dat";

Group {
  // Physicals
  PmlXYZ  = Region[1000];
  PmlXZ   = Region[1001];
  PmlYZ   = Region[1002];
  PmlXY   = Region[1003];
  PmlZ    = Region[1004];
  PmlY    = Region[1005];
  PmlX    = Region[1006];
  Air     = Region[1007];
  Rod     = Region[1008];
  Src     = Region[1009];

  // Domains
  Pml    = Region[{ PmlXYZ, PmlXY, PmlXZ, PmlYZ, PmlX, PmlY, PmlZ }];
  Domain = Region[{ Air, Rod, Pml }];
}

Function {
  Omega0 = 2.0 * Pi * Freq;
  K0     = 2.0 * Pi / Lambda;

  a_pml = 1.0;
  b_pml = 1.0;

  sx[Air]    = 1.0;
  sy[Air]    = 1.0;
  sz[Air]    = 1.0;
  sx[Rod]    = 1.0;
  sy[Rod]    = 1.0;
  sz[Rod]    = 1.0;
  sx[PmlXYZ] = Complex[a_pml, -b_pml];
  sy[PmlXYZ] = Complex[a_pml, -b_pml];
  sz[PmlXYZ] = Complex[a_pml, -b_pml];
  sx[PmlXZ]  = Complex[a_pml, -b_pml];
  sy[PmlXZ]  = 1.0;
  sz[PmlXZ]  = Complex[a_pml, -b_pml];
  sx[PmlYZ]  = 1.0;
  sy[PmlYZ]  = Complex[a_pml, -b_pml];
  sz[PmlYZ]  = Complex[a_pml, -b_pml];
  sx[PmlXY]  = Complex[a_pml, -b_pml];
  sy[PmlXY]  = Complex[a_pml, -b_pml];
  sz[PmlXY]  = 1.0;
  sx[PmlX]   = Complex[a_pml, -b_pml];
  sy[PmlX]   = 1.0;
  sz[PmlX]   = 1.0;
  sx[PmlY]   = 1.0;
  sy[PmlY]   = Complex[a_pml, -b_pml];
  sz[PmlY]   = 1.0;
  sx[PmlZ]   = 1.0;
  sy[PmlZ]   = 1.0;
  sz[PmlZ]   = Complex[a_pml, -b_pml];

  Lxx[] = sy[] * sz[] / sx[];
  Lyy[] = sz[] * sx[] / sy[];
  Lzz[] = sx[] * sy[] / sz[];

  EpsilonRAir[] = Complex[EpsRAirRe, EpsRAirIm];
  EpsilonRRod[] = Complex[EpsRRodRe, EpsRRodIm];

  EpsilonR[Rod] = EpsilonRRod[] * TensorDiag[1.0,   1.0,   1.0];
  EpsilonR[Air] = EpsilonRAir[] * TensorDiag[1.0,   1.0,   1.0];
  EpsilonR[Pml] =                 TensorDiag[Lxx[], Lyy[], Lzz[]];

  NuR[Rod] = TensorDiag[1.0,         1.0,         1.0];
  NuR[Air] = TensorDiag[1.0,         1.0,         1.0];
  NuR[Pml] = TensorDiag[1.0 / Lxx[], 1.0 / Lyy[], 1.0 / Lzz[]];

  source[] = Complex[0, 1] * Vector[0, Omega0 * Mu0, 0];
}

Jacobian {
  { Name Vol;
    Case {
      { Region All; Jacobian Vol; }
    }
  }
  { Name Sur;
    Case {
      { Region All; Jacobian Sur; }
    }
  }
  { Name Lin;
    Case {
      { Region All; Jacobian Lin; }
    }
  }
}

Integration {
  { Name I2;
    Case {
      { Type Gauss;
        Case {
          { GeoElement Point      ; NumberOfPoints 1; }
          { GeoElement Line       ; NumberOfPoints 3; }
          { GeoElement Triangle   ; NumberOfPoints 4; }
          { GeoElement Quadrangle ; NumberOfPoints 4; }
          { GeoElement Tetrahedron; NumberOfPoints 4; }
          { GeoElement Hexahedron ; NumberOfPoints 6; }
          { GeoElement Prism      ; NumberOfPoints 6; }
        }
      }
    }
  }
}

Constraint {
  { Name DirichletBC;
    Case {
      { Region Src; Type AssignFromResolution; NameOfResolution Dirichlet; }
    }
  }
}

FunctionSpace {
  { Name HCurl0; Type Form1;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Edge;
        Support Region[{ Domain, Src }]; Entity EdgesOf[All]; }
      // { Name sn2; NameOfCoef un2; Function BF_Edge_2E;
        // Support Region[Domain]; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef un;  EntityType EdgesOf; NameOfConstraint DirichletBC; }
      // { NameOfCoef un2; EntityType EdgesOf; NameOfConstraint DirichletBC; }
    }
  }
}

Formulation {
  { Name Projection; Type FemEquation;
    Quantity {
      { Name u; Type Local; NameOfSpace HCurl0; }
    }
    Equation {
      Galerkin { [ Dof{u},    {u} ]; In Src; Jacobian Sur; Integration I2; }
      Galerkin { [ -source[], {u} ]; In Src; Jacobian Sur; Integration I2; }
    }
  }

  { Name Maxwell; Type FemEquation;
    Quantity {
      { Name u; Type Local; NameOfSpace HCurl0; }
    }
    Equation {
      Galerkin { [-NuR[] * Dof{Curl u}, {Curl u}];
        In Domain; Jacobian Vol; Integration I2; }
      Galerkin { [(Omega0/C0)^2 * EpsilonR[] * Dof{u}, {u}];
        In Domain; Jacobian Vol; Integration I2; }
    }
  }
}

Resolution {
  { Name Boubouchons;
    System {
      { Name S2; NameOfFormulation Maxwell; Type ComplexValue; Frequency Freq; }
    }
    Operation {
      InitSolution S2;
      Generate S2;
      Solve S2;
      SaveSolution S2;
    }
  }

  { Name Dirichlet;
    System {
      { Name S1; NameOfFormulation Projection; Type ComplexValue;
        Frequency Freq; DestinationSystem S2; }
      { Name S2; NameOfFormulation Maxwell;    Type ComplexValue;
        Frequency Freq;}
    }
    Operation {
      Generate S1;
      Solve S1;
      TransferSolution S1;
    }
  }
}

PostProcessing {
  { Name Post; NameOfFormulation Maxwell;
    Quantity {
      { Name E; Value { Local { [{ u }]; In Domain; Jacobian Vol; } } }
    }
  }
}

PostOperation {
  { Name Post; NameOfPostProcessing Post;
    Operation {
      Print[E, OnElementsOf Domain, File "./boubouchons.pos"];
    }
  }
}
