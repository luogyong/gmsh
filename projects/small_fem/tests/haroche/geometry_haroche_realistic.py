#!/usr/bin/env python
from gmshpy import *
import os

# 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
# 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D)

GmshSetOption('Mesh', 'CharacteristicLengthFactor', 0.9)
GmshSetOption('Mesh', 'Algorithm', 1.0)
GmshSetOption('Mesh', 'Algorithm3D',4.0)

mm = 1.e-3

r=39.4*mm
R=40.6*mm
L_cav=27.57*mm
deltaz=5.*mm
radius_mirror = 25.*mm
dist2PML_xy   =  0.5*mm
dist2PML_z    =  3.7*mm
thick_mirror_atcenter = 1.415*mm

freq     = 51.099e9
cel      = 299792458.
lambdaeg = cel/freq


box_x=radius_mirror+dist2PML_xy
box_y=radius_mirror+dist2PML_xy
box_z=L_cav/2.+thick_mirror_atcenter+dist2PML_z
pml_x=lambdaeg
pml_y=lambdaeg
pml_z=lambdaeg

print R
print (R-r)
print (-R+L_cav/2.)

paramaille = 0.08
myModela = GModel()
myModelb = GModel()
myModel1 = GModel()
myModel2 = GModel()
myModel3 = GModel()
myModel4 = GModel()

myModel1.addCylinder([0,0,-R+L_cav/2.],[0,0,L_cav/2.+thick_mirror_atcenter],radius_mirror)
# Torus - almost spherical - centered in 0,0,0
myModel2.addTorus([0,0,-R+L_cav/2.],[0,1,-R+L_cav/2.],R-r,r)
myModel3.addSphere(0,0,-R+L_cav/2.,r+r/1000)
myModel4.addBlock([0,0,0],[box_x,box_y,box_z])
 
myModel1.computeBooleanDifference(myModel2);
myModel1.computeBooleanDifference(myModel3);
myModel4.computeBooleanDifference(myModel1);



paramaille_pmlx=2
paramaille_pmly=2
paramaille_pmlz=2
myModel4.addVertex(box_x+pml_x,0          ,0,paramaille_pmlx)
myModel4.addVertex(box_x+pml_x,box_y      ,0,paramaille_pmlx)
myModel4.addVertex(box_x+pml_x,box_y+pml_y,0,paramaille_pmlx)
myModel4.addVertex(0          ,box_y+pml_y,0,paramaille_pmly)
myModel4.addVertex(box_x      ,box_y+pml_y,0,paramaille_pmly)

myModel4.addVertex(box_x+pml_x,0          ,box_z,paramaille_pmlx)
myModel4.addVertex(box_x+pml_x,box_y      ,box_z,paramaille_pmlx)
myModel4.addVertex(box_x+pml_x,box_y+pml_y,box_z,paramaille_pmlx)
myModel4.addVertex(0          ,box_y+pml_y,box_z,paramaille_pmly)
myModel4.addVertex(box_x      ,box_y+pml_y,box_z,paramaille_pmly)

myModel4.addVertex(box_x+pml_x,0          ,box_z+pml_z,paramaille_pmlz)
myModel4.addVertex(box_x+pml_x,box_y      ,box_z+pml_z,paramaille_pmlz)
myModel4.addVertex(box_x+pml_x,box_y+pml_y,box_z+pml_z,paramaille_pmlz)
myModel4.addVertex(0          ,box_y+pml_y,box_z+pml_z,paramaille_pmlz)
myModel4.addVertex(box_x      ,box_y+pml_y,box_z+pml_z,paramaille_pmlz)

myModel4.addVertex(0          ,0          ,box_z+pml_z,paramaille_pmlz)
myModel4.addVertex(box_x      ,box_y      ,box_z+pml_z,paramaille_pmlz)
myModel4.addVertex(box_x      ,0          ,box_z+pml_z,paramaille_pmlz)
myModel4.addVertex(0          ,box_y      ,box_z+pml_z,paramaille_pmlz)

myModel4.setAsCurrent();
myModel4.setVisibility(1);
FlGui.instance()
FlGui.run()
FlGui.close()

# 
# 
# ## doesnt work :<(
# # myModel4.addBlock([box_x,0,0],[box_x+pml_x,box_y+pml_y,box_z+pml_z])
# # myModel4.salomeconnect()
# # myModel4.addBlock([0,box_y,0],[box_x+pml_x,box_y+pml_y,box_z+pml_z])
# # myModel4.salomeconnect()
# # myModel4.addBlock([0,0,box_z],[box_x+pml_x,box_y+pml_y,box_z+pml_z])
# # myModel4.salomeconnect()
# # myModel4.addBlock([box_x,box_y,0],[box_x+pml_x,box_y+pml_y,box_z])
# # myModel4.salomeconnect()
# 
# # myModel4.addBlock([0,box_y,0],[box_x,box_y+pml_y,box_z])
# # myModel4.salomeconnect()
# # myModel4.addBlock([box_x,0,0],[box_x+pml_x,box_y,box_z])
# # myModel4.salomeconnect()
# 
# 
# 
# 
# # # e11=myModel1.getEdgeByTag(11)
# # # f7=myModel1.getFaceByTag(7)
# # # v7=myModel1.getVertexByTag(7)
# # # myModel1.remove(f7)
# # # myModel1.remove(e11)
# # # myModel1.remove(v7)
# # 
# # myModel1.computeBooleanDifference(myModel3);
# # myModela.computeBooleanIntersection(myModel1)
# # myModelb.computeBooleanDifference(myModel1)
# # 
# # 
# 
# 
# # myModel1.setVisibility(0);
# # myModel2.setVisibility(0);
# # myModel3.setVisibility(0);
# # myModel4.setVisibility(0);
# # myModelb.setAsCurrent();
# # myModelb.setVisibility(1);
# 
