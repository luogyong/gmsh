// ======
// GROUPS
// ======
Group{
  //propagation domain
  Omega = Region[7];
  //PML
  PML = Region[8];
  // boundary of the scatterers
  Gama = Region[5];
  // fictitious boundary
  Sigma = Region[4];
}

// =========
// FUNCTIONS
// =========
Function {
  // Source
  fs[Gama] = Vector[0, -1, 0];

  // Wavenumber
  k = 1;

  // Max Value
  SigmaMax  = 1;
  SigmaXmax = SigmaMax;
  SigmaYmax = SigmaMax;

  // PML Size
  SizePML  = 8;
  SizePMLX = SizePML;
  SizePMLY = SizePML;

  // Domain
  Xmax = 10;
  Ymax = 10;

  // Distance between a point (X,Y,Z)
  // and the center of the numerical domain (XF,YF,ZF)
  RF_X[] = Sqrt[X[]*X[]];
  RF_Y[] = Sqrt[Y[]*Y[]];
  // Damping functions of the PML: equal to 0 inside the propagation domain
  // and on the intern boundary of the PML
  // (Boundary in common with the Propagation domain).
  //If(PML_TYPE == PML_LINEAR)
  DampingProfileX[] = SigmaXmax/SizePMLX*(RF_X[] - Xmax);
  DampingProfileY[] = SigmaYmax/SizePMLY*(RF_Y[] - Ymax);
  //EndIf
  //If(PML_TYPE == PML_BERMUDEZ)
  //DampingProfileX[] = 1/(Xmax + SizePMLX - Fabs[X[]]) - 1/(SizePMLX);
  //DampingProfileY[] = 1/(Ymax + SizePMLY - Fabs[Y[]]) - 1/(SizePMLY);
  //EndIf
  //If(PML_TYPE == PML_BERMUDEZ_QUAD)
  //DampingProfileX[] = 1/(Xmax + SizePMLX - Fabs[X[]])^2 - 1/(SizePMLX)^2;
  //DampingProfileY[] = 1/(Ymax + SizePMLY - Fabs[Y[]])^2 - 1/(SizePMLY)^2;
  //EndIf
  //Take Max(0, DampingProfile)
  SigmaX[] = 0.5*(DampingProfileX[] + Fabs[DampingProfileX[]]);
  SigmaY[] = 0.5*(DampingProfileY[] + Fabs[DampingProfileY[]]);

  Kx[] = Complex[1, SigmaX[]/k];
  Ky[] = Complex[1, SigmaY[]/k];
  D[]     = TensorDiag[Ky[]/Kx[], Kx[]/Ky[], 0.];
  OverD[] = TensorDiag[Kx[]/Ky[], Ky[]/Kx[], 0.];
}

// ===========
// CONSTRAINTS
// ===========
Constraint{
  //Dirichlet boundary condition on the ficticious boundary
  // (which truncate the PML)
  //{ Name PMLCondition; Type Assign;
  //Case{ {Region Sigma;  Value 0.; } } }

  //Dirichlet boundary condition on Gama, boundary of the scatterers.
  { Name DirichletCondition; Type AssignFromResolution;
    Case{ {Region Gama; NameOfResolution Projection;} } }
}

// =========
// JACOBIAN
// =========
Jacobian {
  { Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
  { Name JSur ; Case { { Region All ; Jacobian Sur ; } } }
}

// ======================
// INTEGRATION PARAMETERS
// ======================
Integration {
  { Name I1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints  1 ; }
          { GeoElement Line        ; NumberOfPoints  4 ; }
          { GeoElement Triangle    ; NumberOfPoints  6 ; }
          { GeoElement Quadrangle  ; NumberOfPoints  7 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 15 ; }
          { GeoElement Hexahedron  ; NumberOfPoints 34 ; }
        }
      }
    }
  }
}

// ==============
// FUNCTION SPACE
// ==============
FunctionSpace{
  // Dirichlet boundary condition
  { Name H_curl; Type Form1;
    BasisFunction{
      { Name se   ; NameOfCoef ee   ; Function BF_Edge_1E ; Support Region[{Omega,PML, Gama, Sigma}] ; Entity EdgesOf[All]; }

      { Name se2e ; NameOfCoef we2e ; Function BF_Edge_2E ; Support Region[{Omega,PML, Gama, Sigma}] ; Entity EdgesOf[All]; }

      { Name se3fa; NameOfCoef we3fa; Function BF_Edge_3F_a ; Support Region[{Omega,PML, Gama, Sigma}] ; Entity FacetsOf[All]; }
      { Name se3fb; NameOfCoef we3fb; Function BF_Edge_3F_b ; Support Region[{Omega,PML, Gama, Sigma}] ; Entity FacetsOf[All]; }
      { Name se3fc; NameOfCoef we3fc; Function BF_Edge_3F_c ; Support Region[{Omega,PML, Gama, Sigma}] ; Entity FacetsOf[All]; }

      { Name se4e ; NameOfCoef we4e ; Function BF_Edge_4E   ; Support Region[{Omega,PML, Gama, Sigma}] ; Entity EdgesOf[All]; }
    }
    Constraint{
      //Dirichlet boundary condition
      { NameOfCoef ee   ; EntityType EdgesOf  ; NameOfConstraint DirichletCondition; }

      { NameOfCoef we2e ; EntityType EdgesOf  ; NameOfConstraint DirichletCondition; }
      { NameOfCoef we3fa; EntityType FacetsOf ; NameOfConstraint DirichletCondition; }
      { NameOfCoef we3fb; EntityType FacetsOf ; NameOfConstraint DirichletCondition; }
      { NameOfCoef we3fc; EntityType FacetsOf ; NameOfConstraint DirichletCondition; }
      { NameOfCoef we4e ; EntityType EdgesOf  ; NameOfConstraint DirichletCondition; }
      /*
      //PML Constraint
      { NameOfCoef ee   ; EntityType EdgesOf  ; NameOfConstraint PMLCondition; }

      { NameOfCoef we2e ; EntityType EdgesOf  ; NameOfConstraint PMLCondition; }
      { NameOfCoef we3fa; EntityType FacetsOf ; NameOfConstraint PMLCondition; }
      { NameOfCoef we3fb; EntityType FacetsOf ; NameOfConstraint PMLCondition; }
      { NameOfCoef we3fc; EntityType FacetsOf ; NameOfConstraint PMLCondition; }
      { NameOfCoef we4e ; EntityType EdgesOf  ; NameOfConstraint PMLCondition; }
      */
    }
  }
}

// ============
// FORMULATIONS
// ============
Formulation {
  // Formulation for a Dirichlet boundary condition
  { Name Dirichlet; Type FemEquation;
    Quantity{
      { Name u ; Type Local; NameOfSpace H_curl; }
    }
    Equation{
      //Maxwell equation
      Galerkin{[Dof{Curl u}, {Curl u}];
        In Omega; Jacobian JVol; Integration I1; }
      Galerkin{[-k^2*Dof{u}, {u}];
        In Omega; Jacobian JVol; Integration I1; }

      //Modified Helmholtz inside the PML
      Galerkin{[OverD[]* Dof{Curl u}, {Curl u}];
        In PML; Jacobian JVol; Integration I1; }
      Galerkin{[-k^2*D[]*Dof{u}, {u}];
        In PML; Jacobian JVol; Integration I1; }
    }
  }

  // Formulation for Projection
  { Name Projection; Type FemEquation;
    Quantity{
      { Name u ; Type Local; NameOfSpace H_curl; }
    }
    Equation {
      Galerkin { [ Dof{u} , {u} ];
        In Gama; Integration I1; Jacobian JSur; }
      Galerkin { [ fs[] , {u} ];
        In Gama; Integration I1; Jacobian JSur; }
    }
  }
}

// ===========
// RESOLUTIONS
// ===========
Resolution{
  { Name Dirichlet;
    System{
      { Name A; NameOfFormulation Dirichlet; Type Complex; }
    }
    Operation{
      Generate[A]; Solve[A];
    }
  }

  {Name Projection;
    System {
      {Name B; NameOfFormulation Projection; DestinationSystem A; }
    }
    Operation {
      Generate[B]; Solve[B]; TransferSolution[B];
    }
  }
}

// ================
// POST-PROCESSINGS
// ================
PostProcessing{
  { Name Dirichlet; NameOfFormulation Dirichlet;
    Quantity {
      { Name uG; Value {Local { [{u}] ; In Omega; Jacobian JVol; } } }
      { Name uP; Value {Local { [{u}] ; In PML;   Jacobian JVol; } } }
    }
  }
}

// ===============
// POST-OPERATIONS
// ===============
PostOperation{
  { Name uG; NameOfPostProcessing Dirichlet ;
    Operation {
      Print [uG, OnElementsOf Omega, File "uG.pos", Depth 3];
      Print [uP, OnElementsOf PML,   File "uP.pos", Depth 3];
    }
  }
}
