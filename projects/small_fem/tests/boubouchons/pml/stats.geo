// Displays stats //
////////////////////

DefineConstant[
  SCelN = {CellN,    Name "Output/00Stats/00Number of Cells",       ReadOnly 1},
  SRodN = {RodN,     Name "Output/00Stats/01Number of Rod cells",   ReadOnly 1},
  SEndN = {EndN * 2, Name "Output/00Stats/02Number of End cells",   ReadOnly 1},
  SPmlN = {PmlN * 2, Name "Output/00Stats/03Number of Pml cells",   ReadOnly 1},
  SDdmN = {DdmN,     Name "Output/00Stats/04Sub-domain size [cell]",ReadOnly 1}
];

SCelN = {};
SRodN = {};
SEndN = {};
SPmlN = {};
SDdmN = {};