Include "parameters_gmsh_getdp.dat";

Group {
	
// Physical Volume(138) =  {1  };    // Air
// 	Physical Volume(139) = {138}; // PMLX
// 	Physical Volume(140) = {126}; // PMLXY
// 	Physical Volume(141) = {124}; // PMLY
// 	Physical Volume(142) = {134}; // PMLZ
// 	Physical Volume(143) = {128}; // PMLXYZ
// 	Physical Volume(144) = {136}; // PMLXZ
// 	Physical Volume(145) = {132}; // PMLYZ
// 
// 	Physical Surface(146) = {2, 70, 72, 120};   // XOZ parallel faces (for sym, apply neumann/diri)
// 	Physical Surface(147) = {1, 112, 118, 114}; // YOZ parallel faces (for sym, apply neumann/diri)
//  Physical Surface(148) = {3};   // mirror 	
	// SubDomains
	Air     = Region[138];
	PMLx    = Region[139];
	PMLxy   = Region[140];
	PMLy    = Region[141];
	PMLz    = Region[142];
	PMLxyz  = Region[143];
	PMLxz   = Region[144];
	PMLyz   = Region[145];
  PMLs        = Region[{PMLxyz,PMLxy,PMLxz,PMLyz,PMLx,PMLy,PMLz}];
  All_domains = Region[{Air,PMLxyz,PMLxy,PMLxz,PMLyz,PMLx,PMLy,PMLz}];

	// Boundaries
	Mirror    = Region[148];
	SurfYZ    = Region[147];
	SurfXZ    = Region[146];
	
	
	PrintPoint     = Region[1000000];
}



Function{
  a_pml          = 1.;
  b_pml          = 1.;
  sx[Air]        = 1.;
  sy[Air]        = 1.;
  sz[Air]        = 1.;
  sx[PMLxyz]     = Complex[a_pml,-b_pml];
  sy[PMLxyz]     = Complex[a_pml,-b_pml];
  sz[PMLxyz]     = Complex[a_pml,-b_pml];
  
  sx[PMLxz]      = Complex[a_pml,-b_pml];
  sy[PMLxz]      = 1.0;
  sz[PMLxz]      = Complex[a_pml,-b_pml];
  
  sx[PMLyz]      = 1.0;
  sy[PMLyz]      = Complex[a_pml,-b_pml];
  sz[PMLyz]      = Complex[a_pml,-b_pml];
  
  sx[PMLxy]      = Complex[a_pml,-b_pml];
  sy[PMLxy]      = Complex[a_pml,-b_pml];
  sz[PMLxy]      = 1.0;
  
  sx[PMLx]       = Complex[a_pml,-b_pml];
  sy[PMLx]       = 1.0;
  sz[PMLx]       = 1.0;
  
  sx[PMLy]       = 1.0;
  sy[PMLy]       = Complex[a_pml,-b_pml];
  sz[PMLy]       = 1.0;
  
  sx[PMLz]       = 1.0;
  sy[PMLz]       = 1.0;
  sz[PMLz]       = Complex[a_pml,-b_pml];
            
	Lxx[]	= sy[]*sz[]/sx[];
	Lyy[] = sz[]*sx[]/sy[]; 
	Lzz[] = sx[]*sy[]/sz[];

	epsilon[Air]  = 1. * TensorDiag[1.,1.,1.];
  epsilon[PMLs]	= 1. * TensorDiag[Lxx[],Lyy[],Lzz[]];
	
  nu[Air]       	= TensorDiag[1.,1.,1.];
  nu[PMLs]     		= TensorDiag[1.0/Lxx[],1.0/Lyy[],1.0/Lzz[]];

	// epsilon[PMLbot] = eps_re_L_4*TensorDiag[Complex[1.,-1.],Complex[1.,-1.],Complex[.5,.5]];
	// epsilon[PMLtop] = eps_re_L_1*TensorDiag[Complex[1.,-1.],Complex[1.,-1.],Complex[.5,.5]];
	// mu[PMLbot]      = TensorDiag[Complex[1.,-1.],Complex[1.,-1.],Complex[.5,.5]];
	// mu[PMLtop]      = TensorDiag[Complex[1.,-1.],Complex[1.,-1.],Complex[.5,.5]];
	// nu[PMLbot]      = TensorDiag[Complex[.5,.5],Complex[.5,.5],Complex[1.,-1.]];
	// nu[PMLtop]      = TensorDiag[Complex[.5,.5],Complex[.5,.5],Complex[1.,-1.]];
	
}

Constraint {
	{Name Dirichlet; Type Assign;
		Case {
			{ Region SurfYZ; Value 0.; }
			{ Region Mirror; Value 0.; }
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
  { Name JLin ;
    Case {
      { Region All ; Jacobian Lin ; }
    }
  }
}

Integration {
  { Name Int_1 ;
    Case { 
      { Type Gauss ;
        Case { 
          { GeoElement Point       ; NumberOfPoints   1 ; }
          { GeoElement Tetrahedron ; NumberOfPoints   4 ; }  // 4 pour l'ordre 1 // 15 pour l'ordre 2
        }
      }
    }
  }
}

FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Edge;
        Support Region[{All_domains}]; Entity EdgesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_Edge_2E;
        Support Region[{All_domains}]; Entity EdgesOf[All]; }
     /* { Name sn3; NameOfCoef un3; Function BF_Edge_3F_b;  //!!rajouter BF_Edge_3F_a
        Support Region[All_domains]; Entity FacetsOf[All_domains]; }
      { Name sn4; NameOfCoef un4; Function BF_Edge_3F_c;
        Support Region[All_domains]; Entity FacetsOf[All_domains]; }
      { Name sn5; NameOfCoef un5; Function BF_Edge_4E;
        Support Region[All_domains]; Entity EdgesOf[All_domains]; }*/
    }
    Constraint {
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      
   /*{ NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
     { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
     { NameOfCoef un5; EntityType  EdgesOf ; NameOfConstraint Dirichlet; }*/
      }
  }
}

Formulation {
	{Name modal_helmholtz_vector; Type FemEquation;
		Quantity {
			{ Name u; Type Local; NameOfSpace Hcurl;}
		}
		Equation { 
			Galerkin {[nu[]*Dof{Curl u} , {Curl u} ];
			In All_domains; Jacobian JVol; Integration Int_1;  }
      
			Galerkin { DtDt[epsilon[]/cel^2 * Dof{u},{u} ];
			In All_domains; Jacobian JVol; Integration Int_1;  }
		}
	}
}

Resolution {
	{ Name all;
		System
		{
			{ Name M1; NameOfFormulation modal_helmholtz_vector; Type ComplexValue; }
		}
		Operation
		{
			GenerateSeparate[M1]; EigenSolve[M1,neig,eig_target,0]; SaveSolutions[M1];
		}
	}
}


PostProcessing {   
	{ Name postpro_eigenvectors; NameOfFormulation modal_helmholtz_vector;
		Quantity {
			{ Name EigenValuesReal;  Value { Local{ [$EigenvalueReal]; In PrintPoint; Jacobian JVol; } } }
			{ Name EigenValuesImag;  Value { Local{ [$EigenvalueImag]; In PrintPoint; Jacobian JVol; } } }
			{ Name u; Value { Local { [ {u} ] ; In All_domains; Jacobian JVol; } } }
			{ Name epsilon_XX; Value { Local { [ CompXX[epsilon[]] ] ; In All_domains; Jacobian JVol; } } }
		}
	}
}



PostOperation {
	{ Name postop_eigenvectors; NameOfPostProcessing postpro_eigenvectors ;
		Operation {
			Print [EigenValuesReal, OnElementsOf PrintPoint, Format TimeTable, File "EigenValuesReal.txt"];
			Print [EigenValuesImag, OnElementsOf PrintPoint, Format TimeTable, File "EigenValuesImag.txt"];
			Print [ u  , OnElementsOf All_domains, File "eigenVectors.pos"       , EigenvalueLegend];
			// Print [ epsilon_XX  , OnElementsOf All_domains, File "epsilon_XX.pos"       , EigenvalueLegend];
		}
	}
}





