Include "parameters_gmsh_getdp.dat";
Group {
	// Physical Line(101) = {1, 2, 3};       // neumann OY
	// Physical Line(102) = {9,10};	        // neumann OX
	// Physical Line(103) = {16,4,11};	      // Diri
	// Physical Line(104) = {7,8,14,15};	    // ext pml
	// Physical Surface(1000) = {20};        // PML_X
	// Physical Surface(2000) = {22};        // PML_XY
	// Physical Surface(3000) = {24};        // PML_Y
	// Physical Surface(4000) = {18};        // air
	// Physical Point(10000) = {2};          // printpoint
	
	
	// Domains
	pmlX           = Region[1000];
	pmlXY          = Region[2000];
	pmlY           = Region[3000];
	air             = Region[4000];
	Omega           = Region[{air,pmlX,pmlY,pmlXY}];
	pml             = Region[{pmlX,pmlY,pmlXY}];
  PrintPoint      = Region[10000];
        
	// Boundaries
	SurfNeu   = Region[{101,102}];
	SurfDirichlet   = Region[103];
	SurfPmlExt   = Region[104];
}

Function{
	
	// PML parameters
	a_pml      =   1.;
	b_pml      =   1.;
	sx[pmlX]   =   Complex[a_pml,-b_pml]; 
	sx[pmlXY]  =   Complex[a_pml,-b_pml]; 
	sx[pmlY]   =   1.;

	sy[pmlX]   =   1.;
	sy[pmlXY]  =   Complex[a_pml,-b_pml]; 
	sy[pmlY]   =   Complex[a_pml,-b_pml]; 

	sz         =   1.;

	epsilonr[air]   = Complex[eps_air_re,eps_air_im] * TensorDiag[1,1,1];
	epsilonr[pml]   = eps_air_re*TensorDiag[sz*sy[]/sx[],sx[]*sz/sy[],sx[]*sy[]/sz];

	mur[air]            = TensorDiag[1,1,1];
	mur[pml]            = TensorDiag[sz*sy[]/sx[],sx[]*sz/sy[],sx[]*sy[]/sz];
}

Constraint {
	{Name Dirichlet; Type Assign;
		Case {
			{ Region SurfDirichlet; Value 0.; }
		}
	}
}


Jacobian {
	{Name JVol ;
		Case { 
			{ Region All ; Jacobian Vol ; }
			}
	}
	{ Name JSur ;
		Case{
			{ Region All ; Jacobian Sur ; }
		}
	}
}

Integration {
{ Name Int_1 ;
    Case { 
      { Type Gauss ;
        Case { 
					{ GeoElement Point       ; NumberOfPoints  1 ; }
					{ GeoElement Line        ; NumberOfPoints  4 ; }
					{ GeoElement Triangle    ; NumberOfPoints  6 ; }
				}
			}
		}
	}
}

FunctionSpace {
  { Name Hgrad; Type Form1P;
    BasisFunction {
      { Name un;  NameOfCoef un;  Function BF_PerpendicularEdge_1N;    Support Omega; Entity NodesOf[Omega]; }
      { Name un2; NameOfCoef un2; Function BF_PerpendicularEdge_2E; Support Omega; Entity EdgesOf[Omega]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }
}

Formulation {
	{
		Name modal_helmholtz_Ppol; Type FemEquation;
		Quantity {
			{ Name u; Type Local; NameOfSpace Hgrad;}
		}
		Equation { 
			Galerkin {[-1/mur[] * Dof{Curl u} , {Curl u}];
			In Omega; Jacobian JVol; Integration Int_1;  }

			// // FOR bloch inside weak form - periodic part search
			// Galerkin {[1/mur[] * Complex[0.,1.] * alpha * Dof{u} , {d u}];
			// In Omega; Jacobian JVol; Integration Int_1;  }
			// 
			// Galerkin {[-1/mur[] * Complex[0.,1.] * alpha * Dof{d u} , {u}];
			// In Omega; Jacobian JVol; Integration Int_1;  }
			// 
			// Galerkin {[-1/mur[] * alpha^2 * Dof{u} , {u}];
			// In Omega; Jacobian JVol; Integration Int_1;  }

			Galerkin { DtDtDof[ -epsilonr[]/cel^2 * Dof{u},{u} ];
			In Omega; Jacobian JVol; Integration Int_1;  }


		}
	}
}
			

Resolution {
	{ Name all;
		System
		{
			{ Name M1; NameOfFormulation modal_helmholtz_Ppol    ; Type ComplexValue; }
		}
		Operation
		{
			GenerateSeparate[M1]; EigenSolve[M1,neig,eig_target,0]; SaveSolutions[M1];
		}
	}
	{ Name eig_only;
		System
		{
			{ Name M2; NameOfFormulation modal_helmholtz_Ppol; Type ComplexValue; }
		}
		Operation
		{
			GenerateSeparate[M2]; EigenSolve[M2,neig,eig_target,0]; //SaveSolutions[M1];
		}
	}
}


PostProcessing {
	{ Name postpro_eigenvectors; NameOfFormulation modal_helmholtz_Ppol;
		Quantity {
		  { Name u ; Value { Local { [     CompZ[{u}]  ] ; In Omega; Jacobian JVol; } } }
		  // { Name u ; Value { Local { [     CompZ[{u}]  ] ; In Omega; Jacobian JVol; } } }
			{ Name EigenValuesReal;  Value { Local{ [$EigenvalueReal]; In PrintPoint; Jacobian JVol; } } }
			{ Name EigenValuesImag;  Value { Local{ [$EigenvalueImag]; In PrintPoint; Jacobian JVol; } } }
		}
	}
}

PostOperation {	
	{ Name postop_eigenvectors; NameOfPostProcessing postpro_eigenvectors ;
		Operation {
			Print [EigenValuesReal, OnElementsOf PrintPoint, Format TimeTable, File "EigenValuesReal_Ppol.txt"];
			Print [EigenValuesImag, OnElementsOf PrintPoint, Format TimeTable, File "EigenValuesImag_Ppol.txt"];
			Print [ u  , OnElementsOf Omega,File "u.pos",EigenvalueLegend];
		}
	}
}