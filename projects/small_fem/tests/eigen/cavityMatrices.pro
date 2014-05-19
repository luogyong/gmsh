Group{
  Border = Region[5]; // Gamma
  Omega  = Region[7]; // Omega
}

Function{
  //Eps[] = TensorDiag[1, 1, 1];
  //Nu[]  = TensorDiag[1, 1, 1];

  //Eps[] = TensorDiag[1, 2, 3];
  //Nu[]  = TensorDiag[1, 2, 3];

  I[] = Complex[0, 1];

  Eps[] = TensorDiag[(X[] + 10) + I[] * (X[] - 10),
                     (Y[] + 10) + I[] * (Y[] - 10),
                     (Z[] + 10) + I[] * (Z[] - 10)];

  Nu[]  = TensorDiag[(X[] + 10) + I[] * (X[] - 10),
                     (Y[] + 10) + I[] * (Y[] - 10),
                     (Z[] + 10) + I[] * (Z[] - 10)];
}

Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
}

Integration {
  { Name IOrder2 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Line ; NumberOfPoints  2 ; }
          { GeoElement Triangle ; NumberOfPoints 3 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 4 ; }
        }
      }
    }
  }

  { Name IOrder4 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Line ; NumberOfPoints  3 ; }
          { GeoElement Triangle ; NumberOfPoints 6 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 15 ; }
        }
      }
    }
  }
}

Constraint {
  { Name Dirichlet_e;
    Case {
      { Region Border ; Value 0 ; }
    }
  }
}

FunctionSpace {
  { Name Hcurl_e ; Type Form1;
    BasisFunction {
      // Ordre 1 Complet //
      { Name se   ; NameOfCoef ee   ; Function BF_Edge    ; Support Region[{Omega,Border}] ; Entity EdgesOf[All]; }
      { Name se2e ; NameOfCoef we2e ; Function BF_Edge_2E ; Support Region[{Omega,Border}] ; Entity EdgesOf[All]; }

      // Ordre 2 Complet //
      { Name se3fa; NameOfCoef we3fa; Function BF_Edge_3F_a ; Support Region[{Omega,Border}] ; Entity FacetsOf[All]; }
      { Name se3fb; NameOfCoef we3fb; Function BF_Edge_3F_b ; Support Region[{Omega,Border}] ; Entity FacetsOf[All]; }
      { Name se3fc; NameOfCoef we3fc; Function BF_Edge_3F_c ; Support Region[{Omega,Border}] ; Entity FacetsOf[All]; }
      { Name se4e ; NameOfCoef we4e ; Function BF_Edge_4E   ; Support Region[{Omega,Border}] ; Entity EdgesOf[All]; }
    }

    Constraint {
      { NameOfCoef ee   ; EntityType EdgesOf  ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we2e ; EntityType EdgesOf  ; NameOfConstraint Dirichlet_e; }

      { NameOfCoef we3fa; EntityType FacetsOf ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we3fb; EntityType FacetsOf ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we3fc; EntityType FacetsOf ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we4e ; EntityType EdgesOf  ; NameOfConstraint Dirichlet_e; }
    }
  }
}

Formulation {
  { Name Maxwell_A; Type FemEquation;
    Quantity {
      { Name e; Type Local;  NameOfSpace Hcurl_e; }
    }
    Equation {
      Galerkin { [ Nu[] * Dof{d e} , {d e} ];
        In Omega; Integration IOrder2; Jacobian JVol;  }
    }
  }

  { Name Maxwell_B; Type FemEquation;
    Quantity {
      { Name e; Type Local;  NameOfSpace Hcurl_e; }
    }
    Equation {
      Galerkin { [ Eps[] * Dof{e} , {e} ];
        In Omega; Integration IOrder4; Jacobian JVol;  }
    }
  }
}

Resolution {
  { Name Maxwell_A ;
    System {
      { Name A ; NameOfFormulation Maxwell_A ; Type Complex ;}
    }
    Operation {
      Generate[A] ; Print[A] ;
    }
  }

  { Name Maxwell_B ;
    System {
      { Name B ; NameOfFormulation Maxwell_B ; Type Complex ;}
    }
    Operation {
      Generate[B] ; Print[B] ;
    }
  }
}
