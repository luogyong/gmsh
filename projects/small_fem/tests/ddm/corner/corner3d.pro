Include "corner3d_data.geo";

DefineConstant[
  // Transmission boundary condition
  TC_TYPE = {1, Name "Input/03DDM/01Transmission condition",
    Choices {0="Order 0", 1="Order 2", 2="Pade (OSRC)"}},
  NP_OSRC = 4,

  // Parameters for the DDM iterative solver
  SOLVER  = "gmres", // bcgs, gmsh_pcleft, ...
  TOL     = 1e-9,
  MAXIT   = 250,
  RESTART = MAXIT
];

Function {
  I[] = Complex[0, 1] ;
  N[] = Normal[] ;
  k   = WAVENUMBER ;
  k[] = k ;

  eps0 = 8.854e-12;
  mu0 = 4*Pi*1e-7;
  c = 1 / Sqrt[mu0*eps0];

  omega[] = c*k[] ;
  mu[] = mu0 ;

  // Source
  ky     = MODE_M * Pi / LY;
  kz     = MODE_N * Pi / LZ;
  kc     = Sqrt[ky^2 + kz^2];
  beta[] = (-kc^2 + k[]^2 >=0 ? Sqrt[-kc^2 + k[]^2] : -I[]*Sqrt[kc^2 - k[]^2]);

  einc[] = Vector[ Sin[ky*Y[]]*Sin[kz*Z[]],
                   I[]*beta[]*ky/kc^2*Cos[ky*Y[]]*Sin[kz*Z[]],
                   I[]*beta[]*kz/kc^2*Cos[kz*Z[]]*Sin[ky*Y[]] ];

  // parameter for ABC
  kInf[]    = k;
  alphaBT[] = 0;
  betaBT[]  = 0;

  // parameter for 0th order TC : IBC(0)
  kDtN[] = k;

  // parameters for 2nd order TC Helmholtz: OO2 Gander 2002, pp. 46-47
  xsimin = 0;
  xsimax = Pi / LC;
  deltak[] = Pi;
  alphastar[] = I[] * ((k^2 - xsimin^2) * (k^2 - (k-deltak[])^2))^(1/4);
  betastar[] = ((xsimax^2 - k^2) * ((k+deltak[])^2 - k^2))^(1/4);
  a[] = - (alphastar[] * betastar[] - k^2) / (alphastar[] + betastar[]);
  b[] = - 1 / (alphastar[] + betastar[]);

  // parameters for 2nd order TC Maxwell: J.-F. Lee
  kmax[] = Pi/LC ;
  delt[] = Sqrt[kmax[]^2-k^2]/Sqrt[k^2];
  Coef_Lee1[] = 1/(1 + I[]*delt[]);
  Coef_Lee2[] = -Coef_Lee1[];

  // parameters for Pade-type TC
  keps[] = k;
  theta_branch = Pi/2;

  Printf("N_DOM %g WAVENUMBER %g N_LAMBDA %g TC_TYPE %g NP_OSRC %g",
         N_DOM, WAVENUMBER, N_LAMBDA, TC_TYPE, NP_OSRC);
}

Group{
  For idom In {0:N_DOM-1}
    Omega~{idom}   = Region[(100 + idom)];
    GammaD0~{idom} = Region[(200 + idom)];
    Kappa~{idom}   = Region[{(10 + idom)}];

    Sigma~{idom}~{0} = Region[{(3000 + idom)}];
    Sigma~{idom}~{1} = Region[{(4000 + idom)}];

    If(idom == 0)
      GammaD~{idom}   = Region[{(1000 + idom)}];
      GammaInf~{idom} = Region[{}];
    EndIf

    If(idom == 3)
      GammaD~{idom}   = Region[{(1000 + idom)}];
      GammaInf~{idom} = Region[{}];
    EndIf

    If(idom == 1)
      GammaD~{idom}   = Region[{}];
      GammaInf~{idom} = Region[{(2000 + idom)}];
    EndIf

    If(idom == 2)
      GammaD~{idom}   = Region[{}];
      GammaInf~{idom} = Region[{(2000 + idom)}];
    EndIf

    Sigma~{idom} = Region[{Sigma~{idom}~{0}, Sigma~{idom}~{1}}];

    BndSigma~{idom}~{0} = Region[{}];
    BndSigma~{idom}~{1} = Region[{}];
    BndSigma~{idom}     = Region[{BndSigma~{idom}~{0}, BndSigma~{idom}~{1}}] ;

    BndGammaInf~{idom}~{0} = Region[{}];
    BndGammaInf~{idom}~{1} = Region[{}];
    BndGammaInf~{idom}     = Region[{BndGammaInf~{idom}~{0},
                                     BndGammaInf~{idom}~{1}}] ;
  EndFor
}

Function{
  // definitions for parallel (MPI) runs:
  ListOfDom = {} ; // the domains that I'm in charge of
  ListOfField = {}; // my fields
  ListOfNeighborField = {}; // my neighbors (left blank: getdp will handle that)

  For idom In {0:N_DOM-1}
    If (idom % MPI_Size == MPI_Rank)
      If(idom == 0)
        myFieldLeft  = {7};
        myFieldRight = {0};

        exchangeFieldLeft  = {6};
        exchangeFieldRight = {1};
      EndIf
      If(idom == 1)
        myFieldLeft  = {1};
        myFieldRight = {2};

        exchangeFieldLeft  = {0};
        exchangeFieldRight = {3};
      EndIf
      If(idom == 2)
        myFieldLeft  = {3};
        myFieldRight = {4};

        exchangeFieldLeft  = {2};
        exchangeFieldRight = {5};
      EndIf
      If(idom == 3)
        myFieldLeft = {5};
        myFieldRight = {6};

        exchangeFieldLeft = {4};
        exchangeFieldRight = {7};
      EndIf
      ListOfDom += idom;
      ListOfField += {myFieldLeft(), myFieldRight()};

      If(ANALYSIS == 0)
        g_in~{idom}~{0}[Sigma~{idom}~{0}] = ComplexScalarField[XYZ[]]{exchangeFieldLeft()};
        g_in~{idom}~{1}[Sigma~{idom}~{1}] = ComplexScalarField[XYZ[]]{exchangeFieldRight()};
      EndIf
      If(ANALYSIS == 1)
        g_in~{idom}~{0}[Sigma~{idom}~{0}] = ComplexVectorField[XYZ[]]{exchangeFieldLeft()};
        g_in~{idom}~{1}[Sigma~{idom}~{1}] = ComplexVectorField[XYZ[]]{exchangeFieldRight()};
      EndIf
    EndIf
  EndFor

  // MPI_Printf["ListOfDom = ", ListOfDom()];
  // MPI_Printf["ListOfField = ", ListOfField()];
  // MPI_Printf["ListOfNeighborField = ", ListOfNeighborField()];
}

If(ANALYSIS == 0)
  Include "Helmholtz.pro" ;
EndIf
If(ANALYSIS == 1)
  Include "Maxwell.pro" ;
EndIf

DefineConstant[
  // default getdp parameters for onelab
  R_ = {"DDM", Name "GetDP/1ResolutionChoices", Visible 0},
  C_ = {"-solve -v 3 -bin -ksp_monitor", Name "GetDP/9ComputeCommand", Visible 0},
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 0}
];
