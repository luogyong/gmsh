Include "cavity_haroche_2D.dat";

Group{
  // Domains
  pmlX          = Region[1000];
  pmlXY         = Region[2000];
  pmlY          = Region[3000];
  air           = Region[4000];
  Omega         = Region[{air, pmlX, pmlY, pmlXY}];
  pml           = Region[{     pmlX, pmlY, pmlXY}];
  PrintPoint    = Region[10000];

  // Boundaries
  SurfLineOY    = Region[101];
  SurfDirichlet = Region[103];
  SurfPmlExt    = Region[104];
}

Function{
  // PML parameters
  a_pml         = 1.;
  b_pml         = 1.;
  sx[pmlX]      = Complex[a_pml, -b_pml];
  sx[pmlXY]     = Complex[a_pml, -b_pml];
  sx[pmlY]      = 1.;

  sy[pmlX]      = 1.;
  sy[pmlXY]     = Complex[a_pml, -b_pml];
  sy[pmlY]      = Complex[a_pml, -b_pml];

  sz            = 1.;

  epsilonr[air] = TensorDiag[1, 1, 1];
  epsilonr[pml] = TensorDiag[sz   * sy[] / sx[],
                             sx[] * sz   / sy[],
                             sx[] * sy[] / sz];

  mur[air]      = TensorDiag[1, 1, 1];
  mur[pml]      = TensorDiag[sz   * sy[] / sx[],
                             sx[] * sz   / sy[],
                             sx[] * sy[] / sz];
}

Constraint{
  { Name Dirichlet; Type Assign;
    Case{
      { Region SurfDirichlet; Value 0.; }
      { Region SurfLineOY   ; Value 0.; }
    }
  }
}


Jacobian{
  { Name JVol;
    Case{
      { Region All; Jacobian Vol; }
    }
  }
  { Name JSur;
    Case{
      { Region All; Jacobian Sur; }
    }
  }
}

Integration{
  { Name Int_1;
    Case{
      { Type Gauss;
        Case{
          { GeoElement Point   ; NumberOfPoints  1; }
          { GeoElement Triangle; NumberOfPoints  3; }
        }
      }
    }
  }
}

FunctionSpace{
  { Name Hcurl; Type Form1;
    BasisFunction{
      { Name e0; NameOfCoef e0; Function BF_Edge_1E; Support Omega;
        Entity EdgesOf[Omega]; }
      // { Name e1; NameOfCoef e1; Function BF_Edge_2E; Support Omega;
      //   Entity EdgesOf[Omega]; }
     }
    Constraint{
      { NameOfCoef e0; EntityType EdgesOf; NameOfConstraint Dirichlet; }
      // { NameOfCoef e1; EntityType EdgesOf; NameOfConstraint Dirichlet; }
    }
  }
}

Formulation{
  {
    Name Maxwell; Type FemEquation;
    Quantity{
      { Name e; Type Local; NameOfSpace Hcurl;}
    }
    Equation{
      Galerkin{ [ 1 / mur[] * Dof{Curl e}, {Curl e} ];
        In Omega; Jacobian JVol; Integration Int_1; }

      Galerkin{ DtDtDof[ epsilonr[] / cel^2 * Dof{e}, {e} ];
        In Omega; Jacobian JVol; Integration Int_1; }
    }
  }
}

Resolution{
  { Name Haroche;
    System{
      { Name M1; NameOfFormulation Maxwell; Type ComplexValue; }
    }
    Operation{
      GenerateSeparate[M1]; EigenSolve[M1, nEig, target, 0];
      SaveSolutions[M1];
    }
  }
}

PostProcessing{
  { Name Post; NameOfFormulation Maxwell;
    Quantity{
      { Name e              ; Value { Local{ [ {e} ]          ; In Omega;
                                             Jacobian JVol; } } }
      { Name EigenValuesReal; Value { Local{ [$EigenvalueReal]; In PrintPoint;
                                             Jacobian JVol; } } }
      { Name EigenValuesImag; Value { Local{ [$EigenvalueImag]; In PrintPoint;
                                             Jacobian JVol; } } }
    }
  }
}

PostOperation{
  { Name Post; NameOfPostProcessing Post;
    Operation{
      Print[EigenValuesReal, OnElementsOf PrintPoint, Format TimeTable,
            File "EigenValuesReal_Ppol.txt"];
      Print[EigenValuesImag, OnElementsOf PrintPoint, Format TimeTable,
            File "EigenValuesImag_Ppol.txt"];
      Print[e  ,             OnElementsOf Omega,
            File "e.pos", EigenvalueLegend];
    }
  }
}
