Include "parameters_gmsh_getdp.dat";

 /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!L4_hh      = L3_hh;*/


Group {
        // SubDomains
                PMLxyz         = Region[1000];
                PMLxz          = Region[1001];
                PMLyz          = Region[1002];
                PMLxy          = Region[1003];
                PMLz           = Region[1004];
                PMLy           = Region[1005];
                PMLx           = Region[1006];
                Scat_In        = Region[1008];
                Scat_Out       = Region[1007];

        // Domains
                Domain         = Region[{Scat_In,Scat_Out}];
                PMLs           = Region[{PMLxyz,PMLxy,PMLxz,PMLyz,PMLx,PMLy,PMLz}];
                All_domains    = Region[{Scat_In,Scat_Out,PMLxyz,PMLxy,PMLxz,PMLyz,PMLx,PMLy,PMLz}];
                
//                SurfNeumann    = Region[{SurfBlochXm,SurfBlochXp,SurfBlochYm,SurfBlochYp}];
}



Function{       
		mu0              = 4*Pi*100.0*nm;
		epsilon0         = 8.854187817e-3*nm;
                cel              = 1.0/Sqrt[epsilon0 * mu0];
                Freq             = cel/lambda0;
                omega0           = 2.0*Pi*Freq;
		k0               = 2.0*Pi/lambda0;
                Ae               = 1.0;
                Ah               = Ae*Sqrt[epsilon0/mu0];
                alpha0           = k0*Sin[theta0]*Cos[phi0];
                beta0            = k0*Sin[theta0]*Sin[phi0];
                gamma0           = k0*Cos[theta0];
	        Ex0              =  Ae * Cos[psi0]*Cos[theta0]*Cos[phi0] - Ae* Sin[psi0]*Sin[phi0];
	        Ey0              =  Ae * Cos[psi0]*Cos[theta0]*Sin[phi0] + Ae* Sin[psi0]*Cos[phi0];
        	Ez0              = -Ae * Cos[psi0]*Sin[theta0];
        	Hx0              = -1/(omega0*mu0)*(beta0  * Ez0 - gamma0 * Ey0);
        	Hy0              = -1/(omega0*mu0)*(gamma0 * Ex0 - alpha0 * Ez0);
        	Hz0              = -1/(omega0*mu0)*(alpha0 * Ey0 - beta0  * Ex0);
        	Prop[]           =  Ae * Complex[ Cos[alpha0*X[]+beta0*Y[]+gamma0*Z[]] , Sin[alpha0*X[]+beta0*Y[]+gamma0*Z[]] ];
        	Einc[PMLs]       =  Vector[0,0,0];
        	Einc[Domain]     =  Vector[Ex0*Prop[],Ey0*Prop[],Ez0*Prop[]];
                Pinc             =  0.5*Ae*Ae*Sqrt[epsilon0/mu0] * Cos[theta0]; 
        	
                a_pml           = 1.;
                b_pml           = 1.;
                sx[Scat_In]          = 1.;
                sy[Scat_In]          = 1.;
                sz[Scat_In]          = 1.;
                sx[Scat_Out]         = 1.;
                sy[Scat_Out]         = 1.;
                sz[Scat_Out]         = 1.;
                sx[PMLxyz]      = Complex[a_pml,-b_pml];
                sy[PMLxyz]      = Complex[a_pml,-b_pml];
                sz[PMLxyz]      = Complex[a_pml,-b_pml];
                
                sx[PMLxz]       = Complex[a_pml,-b_pml];
                sy[PMLxz]       = 1.0;
                sz[PMLxz]       = Complex[a_pml,-b_pml];
                
                sx[PMLyz]       = 1.0;
                sy[PMLyz]       = Complex[a_pml,-b_pml];
                sz[PMLyz]       = Complex[a_pml,-b_pml];
                
                sx[PMLxy]       = Complex[a_pml,-b_pml];
                sy[PMLxy]       = Complex[a_pml,-b_pml];
                sz[PMLxy]       = 1.0;
                
                sx[PMLx]        = Complex[a_pml,-b_pml];
                sy[PMLx]        = 1.0;
                sz[PMLx]        = 1.0;
                
                sx[PMLy]        = 1.0;
                sy[PMLy]        = Complex[a_pml,-b_pml];
                sz[PMLy]        = 1.0;
                
                sx[PMLz]        = 1.0;
                sy[PMLz]        = 1.0;
                sz[PMLz]        = Complex[a_pml,-b_pml];
                
		Lxx[]           = sy[]*sz[]/sx[];
		Lyy[]           = sz[]*sx[]/sy[]; 
		Lzz[] 		= sx[]*sy[]/sz[];

                
                epsilon_In[]    = Complex[eps_re_In  , eps_im_In];
                epsilon_Out[]   = Complex[eps_re_Out , eps_im_Out];
                
                epsilon[Scat_In]   	= epsilon_In[]    * TensorDiag[1.,1.,1.];
                epsilon[Scat_Out]  	= epsilon_Out[]  * TensorDiag[1.,1.,1.];
                epsilon[PMLs]    	= epsilon_Out[]  * TensorDiag[Lxx[],Lyy[],Lzz[]];

                epsilon1[Scat_In]  	= epsilon_Out[]   * TensorDiag[1.,1.,1.];
                epsilon1[Scat_Out] 	= epsilon_Out[]  * TensorDiag[1.,1.,1.];                
                epsilon1[PMLs] 		= epsilon_Out[]  * TensorDiag[Lxx[],Lyy[],Lzz[]];

                mu[Scat_In]        	= TensorDiag[1.,1.,1.];
                mu[Scat_Out]       	= TensorDiag[1.,1.,1.];
                mu[PMLs]      		= TensorDiag[Lxx[],Lyy[],Lzz[]];

                nu[Scat_In]        	= TensorDiag[1.,1.,1.];
                nu[Scat_Out]       	= TensorDiag[1.,1.,1.];
                nu[PMLs]      		= TensorDiag[1.0/Lxx[],1.0/Lyy[],1.0/Lzz[]];

                source[] 	        = (omega0/cel)^2*(epsilon[]-epsilon1[])*Einc[];

}

// Constraint {
//         {Name Dirichlet; Type Assign;
//                 Case {
//                         { Region SurfDirichlet; Value 0.; }
//                 }
//         }
// }

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
          { GeoElement Point       ; NumberOfPoints   4 ; }
          { GeoElement Line        ; NumberOfPoints  32 ; }
          { GeoElement Triangle    ; NumberOfPoints  16 ; } //1, 3, 4, 6, 7, 12, 13, 16
          { GeoElement Tetrahedron ; NumberOfPoints  29 ; }
          { GeoElement Prism       ; NumberOfPoints  51 ; } 
        }
      }
    }
  }
}

FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Edge;
        Support Region[All_domains]; Entity EdgesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_Edge_2E;
        Support Region[All_domains]; Entity EdgesOf[All]; }
      { Name sn3; NameOfCoef un3; Function BF_Edge_3F_b;
        Support Region[All_domains]; Entity FacetsOf[All_domains]; }
      { Name sn4; NameOfCoef un4; Function BF_Edge_3F_c;
        Support Region[All_domains]; Entity FacetsOf[All_domains]; }
      { Name sn5; NameOfCoef un5; Function BF_Edge_4E;
        Support Region[All_domains]; Entity EdgesOf[All_domains]; }
    }
//    Constraint {
//      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint Dirichlet; }
//      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      
/*     { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint BlochX; }*/
/*     { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint BlochY; }*/
/*     { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint Dirichlet; }*/
/*     { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint BlochX; }*/
/*     { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint BlochY; }*/
/*     { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint Dirichlet; }*/
/*     { NameOfCoef un5; EntityType EdgesOf ; NameOfConstraint BlochX; }*/
/*     { NameOfCoef un5; EntityType EdgesOf ; NameOfConstraint BlochY; }*/
/*     { NameOfCoef un5; EntityType EdgesOf ; NameOfConstraint Dirichlet; }*/
//      }
  }
}

Formulation {{Name helmholtz_vector; Type FemEquation;
    		Quantity {{ Name u; Type Local; NameOfSpace Hcurl;}}    
		Equation { Galerkin {[-nu[]*Dof{Curl u} , {Curl u}];
                 		In All_domains; Jacobian JVol; Integration Int_1;  }
                   	   Galerkin { [(omega0/cel)^2*epsilon[]*Dof{u} , {u}];
                 		In All_domains; Jacobian JVol; Integration Int_1;  }
                   	   Galerkin { [source[] , {u}];
                 		In All_domains; Jacobian JVol; Integration Int_1;  }
                }
            }
        }

Resolution {
  { Name helmholtz_vector;
    System {
      { Name M; NameOfFormulation helmholtz_vector; Type ComplexValue; Frequency Freq; }
    }
    Operation { 
      Generate[M]; Solve[M]; SaveSolution[M];  
    }
  }
}



        // Hinc[] : Complex[0,1] * 1/omega0 * 1/mu0 * Curl Einc[];
        // H_d    : Complex[0,1] * 1/omega0 * 1/mu0 * {Curl u};

PostProcessing {   
    { Name get_E; NameOfFormulation helmholtz_vector;NameOfSystem M;
            Quantity {
//             E diffracted 3 components, Im and Re parts

                { Name ex_re_scat; Value { Local { [Re[  CompX[{u}] ]]; In All_domains; Jacobian JVol; } } }
                { Name ey_re_scat; Value { Local { [Re[  CompY[{u}] ]]; In All_domains; Jacobian JVol; } } }
                { Name ez_re_scat; Value { Local { [Re[  CompZ[{u}] ]]; In All_domains; Jacobian JVol; } } }
                { Name ex_im_scat; Value { Local { [Im[  CompX[{u}] ]]; In All_domains; Jacobian JVol; } } }
                { Name ey_im_scat; Value { Local { [Im[  CompY[{u}] ]]; In All_domains; Jacobian JVol; } } }
                { Name ez_im_scat; Value { Local { [Im[  CompZ[{u}] ]]; In All_domains; Jacobian JVol; } } }


                { Name normE2; Value { Local { [ SquNorm[{u}+Einc[]] ]; In All_domains; Jacobian JVol; } } }


/*                { Name sourcex_re; Value { Local { [Re[  CompX[source[]/((omega0/cel)^2)] ]]; In Domain; Jacobian JVol; } } }*/
/*                { Name sourcey_re; Value { Local { [Re[  CompY[source[]/((omega0/cel)^2)] ]]; In Domain; Jacobian JVol; } } }*/
/*                { Name sourcez_re; Value { Local { [Re[  CompZ[source[]/((omega0/cel)^2)] ]]; In Domain; Jacobian JVol; } } }*/
/*                { Name sourcex_im; Value { Local { [Im[  CompX[source[]/((omega0/cel)^2)] ]]; In Domain; Jacobian JVol; } } }*/
/*                { Name sourcey_im; Value { Local { [Im[  CompY[source[]/((omega0/cel)^2)] ]]; In Domain; Jacobian JVol; } } }*/
/*                { Name sourcez_im; Value { Local { [Im[  CompZ[source[]/((omega0/cel)^2)] ]]; In Domain; Jacobian JVol; } } }
*/

/*                { Name epsXX_re; Value { Local { [Re[  CompXX[ epsilon[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name epsXX_im; Value { Local { [Im[  CompXX[ epsilon[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name epsYY_re; Value { Local { [Re[  CompYY[ epsilon[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name epsYY_im; Value { Local { [Im[  CompYY[ epsilon[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name epsZZ_re; Value { Local { [Re[  CompZZ[ epsilon[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name epsZZ_im; Value { Local { [Im[  CompZZ[ epsilon[] ] ]]; In All_domains; Jacobian JVol; } } }*/

/*                { Name nuXX_re; Value { Local { [Re[  CompXX[ nu[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name nuXX_im; Value { Local { [Im[  CompXX[ nu[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name nuYY_re; Value { Local { [Re[  CompYY[ nu[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name nuYY_im; Value { Local { [Im[  CompYY[ nu[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name nuZZ_re; Value { Local { [Re[  CompZZ[ nu[] ] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name nuZZ_im; Value { Local { [Im[  CompZZ[ nu[] ] ]]; In All_domains; Jacobian JVol; } } }
*/
		
		{ Name   dummy1; Value { Local { [Re[  CompXX[epsilon[]] ]]; In All_domains; Jacobian JVol; } } }
                
/*                { Name test; Value {  Pinc; In Domain; Jacobian JVol; } } }*/

                { Name normalized_losses ; Value { Integral { [ 0.5*omega0*epsilon0*Fabs[Im[epsilon_In[]]]*(SquNorm[{u}+Einc[]])  ] ; In Scat_In ; Integration Int_1 ; Jacobian JVol ; } } }
                { Name inc ; Value { Local { [ (Pinc*Pi*ro*ro) ] ; In Scat_In ; Integration Int_1 ; Jacobian JVol ; } } }
/*                / (Pinc*Pi*ro*ro)*/


/*sigma_rode              = omega*abs(imag(epsilon0*eps_rode));*/

/*    int_sigmaE2(i_R)  = postint(fem,'1/2*sigma_rode*(abs(Ex+Ex0)^2+abs(Ey+Ey0)^2+abs(Ez+Ez0)^2)','unit','','dl',[19]);*/
/*    P_inc_enW(i_R)    = k0/(omega*mu0)*(pi*(R_particule)^2)*cos(theta_0);*/
/*    AAQ(i_R)          = int_sigmaE2(i_R)/P_inc_enW(i_R);*/
/*    absorptivity(i_R) = AAQ(i_R)/(pi*(R_particule)^2);*/



/*                { Name ex_cm_re; Value { Local { [Re[  CompX[Einc[]] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name ex_cm_im; Value { Local { [Im[  CompX[Einc[]] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name ey_cm_re; Value { Local { [Re[  CompY[Einc[]] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name ey_cm_im; Value { Local { [Im[  CompY[Einc[]] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name ez_cm_re; Value { Local { [Re[  CompZ[Einc[]] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name ez_cm_im; Value { Local { [Im[  CompZ[Einc[]] ]]; In All_domains; Jacobian JVol; } } }*/

/*                { Name dummy1; Value { Local { [Re[  CompXX[epsilon[]] ]]; In All_domains; Jacobian JVol; } } }*/
/*                { Name dummy2; Value { Local { [Re[  CompXX[epsilon1[]] ]]; In All_domains; Jacobian JVol; } } }*/

                
                
/*//             H diffracted 3 components, Im and Re parts                 */
/*                { Name hx_d_re; Value { Local { [ -Im[CompX[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name hx_d_im; Value { Local { [  Re[CompX[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name hy_d_re; Value { Local { [ -Im[CompY[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name hy_d_im; Value { Local { [  Re[CompY[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name hz_d_re; Value { Local { [ -Im[CompZ[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name hz_d_im; Value { Local { [  Re[CompZ[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*//             Poynting diffracted field => used in reflexion coefficient (real quantity)*/
/*                { Name pox; Value { Local { [ 0.5*CompX[ Im[ 1/omega0 * 1/mu0*{u}/\Conj[{Curl u}] ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name poy; Value { Local { [ 0.5*CompY[ Im[ 1/omega0 * 1/mu0*{u}/\Conj[{Curl u}] ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name poz; Value { Local { [ 0.5*CompZ[ Im[ 1/omega0 * 1/mu0*{u}/\Conj[{Curl u}] ] ] ]; In All_domains; Jacobian JVol; } } }*/
/*//             To be integrated via GMSH Integrate Plugin (use the cutplane view) */
/*                { Name poz_r_int           ; Value {    Local { [ 1/period_x * 1/period_y * CompZ[ 0.5*Re[ Complex[0,1]* 1/omega0 * 1/mu0 * {u}/\Conj[{Curl u}] ] ] / Pinc ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name poz_t_int           ; Value {    Local { [ 1/period_x * 1/period_y * CompZ[ 0.5*Re[ (Einc[]+{u})/\Conj[ Complex[0,1]/(omega0*mu0)*{Curl u} + Hinc[]] ] ] / Pinc ]; In All_domains; Jacobian JVol; } } }*/
/*                { Name joule_losses        ; Value {    Local { [  epsilon0*omega0 * 0.5*Fabs[Im[epsilon_rode[]]]*(SquNorm[{u}+Einc[]]) / (Pinc*period_x*period_y)  ] ; In All_domains; Jacobian JVol; } } }*/
/*//             R (using E_d) then T (using E_d+E_tot). Makes use of the fact that we know a priori that the efficiencies are constant along the substrate */
/*                {  Name int_vol_div_P_Super   ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*1/(h_supb-h_supa)*CompZ[ 0.5*Re[ Complex[0,1]* 1/omega0 * 1/mu0 * {u}/\Conj[{Curl u}] ] ]   ]             ; In Super ; Integration Int_1 ; Jacobian JVol ; } } }*/
/*                {  Name int_vol_div_P_Subs    ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*1/(h_subb-h_suba)*CompZ[ 0.5*Re[ (Einc[]+{u})/\Conj[ Complex[0,1]/(omega0*mu0)*{Curl u} + Hinc[]] ] ]   ] ; In Subs  ; Integration Int_1 ; Jacobian JVol ; } } }*/
/*//             R (using E_d) then T (using E_d+E_tot) with 2 plane cuts in the superstrate and substrate*/
/*//                 {  Name int_sur_div_P_Sup1 ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*CompZ[ 0.5*Re[ Complex[0,1]* 1/omega0 * 1/mu0 * {u}/\Conj[{Curl u}] ] ]   ]             ; In SurfIntegSup1 ; Integration Int_1 ; Jacobian JSur ; } } }*/
/*//                 {  Name int_sur_div_P_Sup2 ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*CompZ[ 0.5*Re[ Complex[0,1]* 1/omega0 * 1/mu0 * {u}/\Conj[{Curl u}] ] ]   ]             ; In SurfIntegSup2 ; Integration Int_1 ; Jacobian JSur ; } } }*/
/*//                 {  Name int_sur_div_P_Sub1 ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*CompZ[ 0.5*Re[ (Einc[]+{u})/\Conj[ Complex[0,1]/(omega0*mu0)*{Curl u} + Hinc[]] ] ]   ] ; In SurfIntegSub1 ; Integration Int_1 ; Jacobian JSur ; } } }*/
/*//                 {  Name int_sur_div_P_Sub2 ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*CompZ[ 0.5*Re[ (Einc[]+{u})/\Conj[ Complex[0,1]/(omega0*mu0)*{Curl u} + Hinc[]] ] ]   ] ; In SurfIntegSub2 ; Integration Int_1 ; Jacobian JSur ; } } }*/
/*//             Computes the volume of a region */
/*                {  Name blabla              ; Value { Integral { [1] ; In Rode ; Integration Int_1 ; Jacobian JVol ; } } }*/
                
}}}



PostOperation {
{ Name Ed; NameOfPostProcessing get_E ;
    Operation {
        Print [ ex_re_scat  , OnElementsOf All_domains, File "./Views/Ex_re_scat.pos", Smoothing ];
        Print [ ey_re_scat  , OnElementsOf All_domains, File "./Views/Ey_re_scat.pos", Smoothing ];
        Print [ ez_re_scat  , OnElementsOf All_domains, File "./Views/Ez_re_scat.pos", Smoothing ];
        Print [ ex_im_scat  , OnElementsOf All_domains, File "./Views/Ex_im_scat.pos", Smoothing ];
        Print [ ey_im_scat  , OnElementsOf All_domains, File "./Views/Ey_im_scat.pos", Smoothing ];
        Print [ ez_im_scat  , OnElementsOf All_domains, File "./Views/Ez_im_scat.pos", Smoothing ];
        Print [ inc         , OnElementsOf All_domains, File "./Views/inc.pos", Smoothing ];
        Print [ normE2      , OnElementsOf All_domains, File "./Views/normE2.pos", Smoothing ];
        
/*        */
/*        */
/*        Print [ sourcex_re  , OnElementsOf All_domains, File "./Views/sourcex_re.pos"];*/
/*        Print [ sourcey_re  , OnElementsOf All_domains, File "./Views/sourcey_re.pos"];*/
/*        Print [ sourcez_re  , OnElementsOf All_domains, File "./Views/sourcez_re.pos"];*/
/*        Print [ sourcex_im  , OnElementsOf All_domains, File "./Views/sourcex_im.pos"];*/
/*        Print [ sourcey_im  , OnElementsOf All_domains, File "./Views/sourcey_im.pos"];*/
/*        Print [ sourcez_im  , OnElementsOf All_domains, File "./Views/sourcez_im.pos"];*/

/*        Print [ epsXX_re  , OnElementsOf PMLs, File "./Views/epsXX_re.pos"];*/
/*        Print [ epsXX_im  , OnElementsOf PMLs, File "./Views/epsXX_im.pos"];*/
/*        Print [ epsYY_re  , OnElementsOf PMLs, File "./Views/epsYY_re.pos"];*/
/*        Print [ epsYY_im  , OnElementsOf PMLs, File "./Views/epsYY_im.pos"];*/
/*        Print [ epsZZ_re  , OnElementsOf PMLs, File "./Views/epsZZ_re.pos"];        */
/*        Print [ epsZZ_im  , OnElementsOf PMLs, File "./Views/epsZZ_im.pos"];*/

/*        Print [ nuXX_re  , OnElementsOf PMLs, File "./Views/nuXX_re.pos"];*/
/*        Print [ nuXX_im  , OnElementsOf PMLs, File "./Views/nuXX_im.pos"];*/
/*        Print [ nuYY_re  , OnElementsOf PMLs, File "./Views/nuYY_re.pos"];*/
/*        Print [ nuYY_im  , OnElementsOf PMLs, File "./Views/nuYY_im.pos"];*/
/*        Print [ nuZZ_re  , OnElementsOf PMLs, File "./Views/nuZZ_re.pos"];        */
/*        Print [ nuZZ_im  , OnElementsOf PMLs, File "./Views/nuZZ_im.pos"];*/

/*        Print [ test  , OnElementsOf PMLs, File "./Views/test.pos"];*/
/*        Print [ dummy1  , OnElementsOf PMLs, File "./Views/dummy1.pos"];*/
        
        Print[ normalized_losses[Scat_In]           , OnGlobal, File "./Views/temp-Q.txt", Format Table ];





        //Print[ int_vol_div_P_Subs[Subs]          , OnGlobal, File "temp.txt"  , Format Table ];
        //Print[ int_vol_div_P_Super[Super]        , OnGlobal, File > "temp.txt", Format Table ];
/*        Print [ poz_r_int, OnElementsOf All_domains, File "poz_r_int.pos"]; //, Smoothing ];*/
/*        Print [ poz_t_int, OnElementsOf All_domains, File "poz_t_int.pos"]; //, Smoothing ];*/
/*        Print [ dummy1  , OnElementsOf All_domains, File "dummy1.pos"];*/
/*        Print [ dummy2  , OnElementsOf All_domains, File "dummy2.pos"];*/
        
/*        Print [ ex_d_im  , OnElementsOf All_domains, File "map_imEX_diffacted.pos"];*/
/*        Print [ ey_d_im  , OnElementsOf All_domains, File "map_imEY_diffacted.pos"];*/
/*        Print [ ez_d_im  , OnElementsOf All_domains, File "map_imEZ_diffacted.pos"];*/
/*        Print [ ex_cm_re  , OnElementsOf All_domains, File "map_reEX_cm.pos"];*/
/*        Print [ ey_cm_re  , OnElementsOf All_domains, File "map_reEY_cm.pos"];*/
/*        Print [ ez_cm_re  , OnElementsOf All_domains, File "map_reEZ_cm.pos"];*/
/*        Print [ ex_cm_im  , OnElementsOf All_domains, File "map_imEX_cm.pos"];*/
/*        Print [ ey_cm_im  , OnElementsOf All_domains, File "map_imEY_cm.pos"];*/
/*        Print [ ez_cm_im  , OnElementsOf All_domains, File "map_imEZ_cm.pos"];*/
/*        Print [ hy_d_re  , OnElementsOf All_domains, File "map_reHY_diffacted.pos"]; //, Smoothing ];*/
        }
}
}



