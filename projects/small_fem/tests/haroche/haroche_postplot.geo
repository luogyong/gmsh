Merge "eigenVectors.pos" ;
View[0].Visible = 0;

Plugin(MathEval).View = 0;
Plugin(MathEval).Expression0 = 'v0' ;
Plugin(MathEval).Run ;
View[1].Visible = 0;

Plugin(CutPlane).A = 1; 
Plugin(CutPlane).B = 0; 
Plugin(CutPlane).C = 0; 
Plugin(CutPlane).D = 0; 
Plugin(CutPlane).View = 1 ;
Plugin(CutPlane).Run ; 

Plugin(CutPlane).A = 0; 
Plugin(CutPlane).B = 1; 
Plugin(CutPlane).C = 0; 
Plugin(CutPlane).D = 0; 
Plugin(CutPlane).View = 1 ;
Plugin(CutPlane).Run ; 

Plugin(CutPlane).A = 0; 
Plugin(CutPlane).B = 0; 
Plugin(CutPlane).C = 1; 
Plugin(CutPlane).D = 0; 
Plugin(CutPlane).View = 1 ;
Plugin(CutPlane).Run ; 
