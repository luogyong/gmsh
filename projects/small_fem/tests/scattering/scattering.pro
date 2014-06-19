Include "parameters_gmsh_getdp.dat";

Group {
  // SubDomains
  PMLxyz      = Region[1000];
  PMLxz       = Region[1001];
  PMLyz       = Region[1002];
  PMLxy       = Region[1003];
  PMLz        = Region[1004];
  PMLy        = Region[1005];
  PMLx        = Region[1006];
  Scat_In     = Region[1008];
  Scat_Out    = Region[1007];

  // Domains
  Domain      = Region[{Scat_In,Scat_Out}];
  PMLs        = Region[{PMLxyz,PMLxy,PMLxz,PMLyz,PMLx,PMLy,PMLz}];
  All_domains = Region[{Scat_In,Scat_Out,
                        PMLxyz,PMLxy,PMLxz,PMLyz,PMLx,PMLy,PMLz}];
}

Function{
  mu0             = 4 * Pi * 100.0 * nm;
  epsilon0        = 8.854187817e-3 * nm;
  cel             = 1.0 / Sqrt[epsilon0 * mu0];
  Freq            = cel / lambda0;
  omega0          = 2.0 * Pi * Freq;
  k0              = 2.0 * Pi / lambda0;

  alpha0          = k0 * Sin[theta0] * Cos[phi0];
  beta0           = k0 * Sin[theta0] * Sin[phi0];
  gamma0          = k0 * Cos[theta0];

  Ae              = 1.0;
  Ex0             =  Ae * Cos[psi0] * Cos[theta0] * Cos[phi0] -
                     Ae * Sin[psi0] * Sin[phi0];
  Ey0             =  Ae * Cos[psi0] * Cos[theta0] * Sin[phi0] +
                     Ae * Sin[psi0] * Cos[phi0];
  Ez0             = -Ae * Cos[psi0] * Sin[theta0];

  Prop[]          =  Ae * Complex[Cos[alpha0*X[] + beta0*Y[] + gamma0*Z[]] ,
                                  Sin[alpha0*X[] + beta0*Y[] + gamma0*Z[]]];
  Einc[PMLs]      =  Vector[0,0,0];
  Einc[Domain]    =  Vector[Ex0 * Prop[], Ey0 * Prop[], Ez0 * Prop[]];

  a_pml           = 1.;
  b_pml           = 1.;
  sx[Scat_In]     = 1.;
  sy[Scat_In]     = 1.;
  sz[Scat_In]     = 1.;
  sx[Scat_Out]    = 1.;
  sy[Scat_Out]    = 1.;
  sz[Scat_Out]    = 1.;
  sx[PMLxyz]      = Complex[a_pml, -b_pml];
  sy[PMLxyz]      = Complex[a_pml, -b_pml];
  sz[PMLxyz]      = Complex[a_pml, -b_pml];
  sx[PMLxz]       = Complex[a_pml, -b_pml];
  sy[PMLxz]       = 1.0;
  sz[PMLxz]       = Complex[a_pml, -b_pml];
  sx[PMLyz]       = 1.0;
  sy[PMLyz]       = Complex[a_pml, -b_pml];
  sz[PMLyz]       = Complex[a_pml, -b_pml];
  sx[PMLxy]       = Complex[a_pml, -b_pml];
  sy[PMLxy]       = Complex[a_pml, -b_pml];
  sz[PMLxy]       = 1.0;
  sx[PMLx]        = Complex[a_pml, -b_pml];
  sy[PMLx]        = 1.0;
  sz[PMLx]        = 1.0;
  sx[PMLy]        = 1.0;
  sy[PMLy]        = Complex[a_pml, -b_pml];
  sz[PMLy]        = 1.0;
  sx[PMLz]        = 1.0;
  sy[PMLz]        = 1.0;
  sz[PMLz]        = Complex[a_pml, -b_pml];
  Lxx[]           = sy[] * sz[] / sx[];
  Lyy[]           = sz[] * sx[] / sy[];
  Lzz[]           = sx[] * sy[] / sz[];


  epsilon_In[]       = Complex[eps_re_In , eps_im_In];
  epsilon_Out[]      = Complex[eps_re_Out, eps_im_Out];

  epsilon[Scat_In]   = epsilon_In[]  * TensorDiag[1., 1., 1.];
  epsilon[Scat_Out]  = epsilon_Out[] * TensorDiag[1., 1., 1.];
  epsilon[PMLs]      = epsilon_Out[] * TensorDiag[Lxx[], Lyy[], Lzz[]];

  epsilon1[Scat_In]  = epsilon_Out[] * TensorDiag[1., 1., 1.];
  epsilon1[Scat_Out] = epsilon_Out[] * TensorDiag[1., 1., 1.];
  epsilon1[PMLs]     = epsilon_Out[] * TensorDiag[Lxx[], Lyy[], Lzz[]];

  nu[Scat_In]        = TensorDiag[1., 1., 1.];
  nu[Scat_Out]       = TensorDiag[1., 1., 1.];
  nu[PMLs]           = TensorDiag[1.0 / Lxx[], 1.0 / Lyy[], 1.0 / Lzz[]];

  source[]           = (omega0 / cel)^2 * (epsilon[]-epsilon1[]) * Einc[];
}

Jacobian{
  {Name JVol;
    Case{
      {Region All; Jacobian Vol;}
    }
  }
}

Integration{
  {Name Int;
    Case{
      {Type Gauss;
        Case{
          {GeoElement Point      ; NumberOfPoints   1;}
          {GeoElement Line       ; NumberOfPoints   4;}
          {GeoElement Triangle   ; NumberOfPoints   6;}
          {GeoElement Tetrahedron; NumberOfPoints  15;}
        }
      }
    }
  }
}

FunctionSpace{
  {Name Hcurl; Type Form1;
    BasisFunction{
      {Name sn; NameOfCoef un; Function BF_Edge;
       Support Region[All_domains]; Entity EdgesOf[All];}
      // {Name sn2; NameOfCoef un2; Function BF_Edge_2E;
      //  Support Region[All_domains]; Entity EdgesOf[All];}
      // {Name sn3; NameOfCoef un3; Function BF_Edge_3F_b;
      //  Support Region[All_domains]; Entity FacetsOf[All];}
      // {Name sn4; NameOfCoef un4; Function BF_Edge_3F_c;
      //  Support Region[All_domains]; Entity FacetsOf[All];}
      // {Name sn5; NameOfCoef un5; Function BF_Edge_4E;
      //  Support Region[All_domains]; Entity EdgesOf[All];}
    }
  }
}

Formulation{
  {Name helmholtz_vector; Type FemEquation;
    Quantity{
      {Name u; Type Local; NameOfSpace Hcurl;}
    }
    Equation{
      Galerkin{[nu[] * Dof{Curl u}, {Curl u}];
        In All_domains; Jacobian JVol; Integration Int;}
      Galerkin{[-(omega0/cel)^2 * epsilon[] * Dof{u}, {u}];
        In All_domains; Jacobian JVol; Integration Int;}
      Galerkin{[-source[], {u}];
        In All_domains; Jacobian JVol; Integration Int;}
    }
  }
}

Resolution{
  {Name helmholtz_vector;
    System{
      {Name M; NameOfFormulation helmholtz_vector;
       Type ComplexValue; Frequency Freq;}
    }
    Operation{
      Generate[M]; Solve[M]; SaveSolution[M];
    }
  }
}

PostProcessing{
  {Name get_E; NameOfFormulation helmholtz_vector;NameOfSystem M;
    Quantity{
      {Name E_scat; Value{Local{[{u}];        In All_domains; Jacobian JVol;}}}
      {Name E_tot ; Value{Local{[{u}+Einc[]]; In All_domains; Jacobian JVol;}}}
    }
  }
}

PostOperation{
{Name Ed; NameOfPostProcessing get_E;
  Operation{
      Print[E_scat, OnElementsOf All_domains, File "./E_scat.pos"];
      Print[E_tot , OnElementsOf All_domains, File "./E_tot.pos"];
  }
 }
}
