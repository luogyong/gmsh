Group {
  Volume = Region[7];
  Border = Region[5];
  Domain = Region[{ Volume, Border }];
  PrintP  = Region[100];
}

Function {
  c0   = 299792458;
  mu0  = 4 * Pi * 1e-7;
  eps0 = 1.0 / (mu0 * c0 * c0);
  nu0  = 1.0 / mu0;
  mm   = 1e3; // Scaling from milimeters

  nEig = 10;
  //eTarget = 1e-4;  // Wavenumber squared
  eTarget = 6.67e19; // Angular frequency squared
}

Constraint {
  { Name Dirichlet; Type Assign;
    Case {{ Region Border; Value 0.; }}
  }
}

Jacobian {
  { Name JVol; Case {{ Region All; Jacobian Vol; }}}
  { Name JSur; Case {{ Region All; Jacobian Sur; }}}
  { Name JLin; Case {{ Region All; Jacobian Lin; }}}
}

Integration {
  { Name Int;
    Case {
      { Type Gauss;
        Case {
          { GeoElement Point      ; NumberOfPoints   1; }
          { GeoElement Triangle   ; NumberOfPoints   3; }
          { GeoElement Tetrahedron; NumberOfPoints   4; }
        }
      }
    }
  }
}

FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Edge;
        Support Region[{Domain}]; Entity EdgesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_Edge_2E;
        Support Region[{Domain}]; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }
}

Formulation {
  {Name Tesla; Type FemEquation;
    Quantity {
      { Name e; Type Local; NameOfSpace Hcurl;}
    }
    Equation {
      Galerkin {[nu0 * mm * Dof{Curl e} , {Curl e}];
        In Domain; Jacobian JVol; Integration Int; }

      Galerkin {DtDt[eps0 / mm * Dof{e}, {e}];
        In Domain; Jacobian JVol; Integration Int;  }
    }
  }
}

Resolution {
  { Name EVP;
    System {
      { Name M1; NameOfFormulation Tesla; Type ComplexValue; }
    }
    Operation {
      GenerateSeparate[M1];
      EigenSolve[M1, nEig, eTarget, 0];
      SaveSolutions[M1];
    }
  }
}

PostProcessing {
  { Name Tesla; NameOfFormulation Tesla;
    Quantity {
      { Name EigenValuesReal;
        Value { Local{ [$EigenvalueReal]; In PrintP; Jacobian JVol; } } }
      { Name EigenValuesImag;
        Value { Local{ [$EigenvalueImag]; In PrintP; Jacobian JVol; } } }
      { Name e; Value { Local { [ {e} ] ; In Domain; Jacobian JVol; } } }
    }
  }
}

PostOperation {
  { Name Tesla; NameOfPostProcessing Tesla;
    Operation {
      Print[EigenValuesReal, OnElementsOf PrintP, Format TimeTable,
             File "EigenValuesReal.txt"];
      Print[EigenValuesImag, OnElementsOf PrintP, Format TimeTable,
             File "EigenValuesImag.txt"];
      Print[e, OnElementsOf Domain,
             File "eigenVectors.pos", EigenvalueLegend];
    }
  }
}
