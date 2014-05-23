#!/usr/bin/env python
from gmshpy import *
import numpy

# 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
# 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D)

GmshSetOption('Mesh', 'CharacteristicLengthFactor', 0.9)
GmshSetOption('Mesh', 'Algorithm', 1.0)
GmshSetOption('Mesh', 'Algorithm3D',4.0)

## Pyhsics ##
#############
## Scaling
mm = 1.e-3
nm = 1.e-9

## Speed of light
epsilon0 = 8.854187817e-3 * nm
mu0      = 400. * numpy.pi * nm
c        = 1.0 / (numpy.sqrt(epsilon0 * mu0))

## Haroche Frequency & Wavelength
freq_haroche   = 51.099e9
lambda_haroche = c / freq_haroche


## Geomtrical Parameters ##
###########################
## Mirror
r                     = 39.4   * mm
R                     = 40.6   * mm
L_cav                 = 27.57  * mm
deltaz                =  5.    * mm
radius_mirror         = 25.    * mm
dist2PML_xy           =  0.5   * mm
dist2PML_z            =  3.7   * mm
thick_mirror_atcenter =  1.415 * mm

## Box
box_x =                      radius_mirror + dist2PML_xy
box_y =                      radius_mirror + dist2PML_xy
box_z = L_cav / 2. + thick_mirror_atcenter + dist2PML_z

## PML
pml_x = lambda_haroche
pml_y = lambda_haroche
pml_z = lambda_haroche


## Mesh Parameters ##
#####################
## Mesh Size
refinement      = 3.

paramaille_air  = lambda_haroche / refinement
paramaille_pml  = lambda_haroche / refinement
paramaille_mir  = lambda_haroche / refinement


## Geomtrical Model ##
######################
myModel1 = GModel()
myModel2 = GModel()
myModel3 = GModel()
myModel4 = GModel()

myModel1.addCylinder([0, 0, -R + L_cav / 2.], \
                     [0, 0,  L_cav / 2. + thick_mirror_atcenter], \
                     radius_mirror)

## Torus - almost spherical - centered in 0,0,0
myModel2.addTorus([0, 0, -R + L_cav / 2.], [0, 1, -R + L_cav / 2.], R-r, r)
myModel3.addSphere(0, 0, -R + L_cav / 2., r+r / 1000)
myModel4.addBlock([0, 0, 0], [box_x, box_y, box_z])

myModel1.computeBooleanDifference(myModel2);
myModel1.computeBooleanDifference(myModel3);
myModel4.computeBooleanDifference(myModel1);

## Extrusion of boundary surface to get PML region
face = myModel4.bindingsGetFaces()
myModel4.extrude(face[8] , [0, 0, 0], [pml_x, 0, 0])

face = myModel4.bindingsGetFaces()
myModel4.extrude(face[11], [0, 0, 0], [0, 0, pml_z])

face = myModel4.bindingsGetFaces()
myModel4.extrude(face[7] , [0, 0, 0], [0, 0, pml_z])

face = myModel4.bindingsGetFaces()
myModel4.extrude(face[22], [0, 0, 0], [0, pml_y, 0])

face = myModel4.bindingsGetFaces()
myModel4.extrude(face[13], [0, 0, 0], [0, pml_y, 0])

face = myModel4.bindingsGetFaces()
myModel4.extrude(face[24], [0, 0, 0], [0, pml_y, 0])

face = myModel4.bindingsGetFaces()
myModel4.extrude(face[18], [0, 0, 0], [0, pml_y, 0])

## Mesh size at Vertices
vertex = myModel4.bindingsGetVertices()

for i in range(1, len(vertex)):
    myModel4.getVertexByTag(i).setPrescribedMeshSizeAtVertex(paramaille_pml)

vertex[36].setPrescribedMeshSizeAtVertex(paramaille_mir)
vertex[37].setPrescribedMeshSizeAtVertex(paramaille_mir)
vertex[41].setPrescribedMeshSizeAtVertex(paramaille_mir)

vertex[0].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[10].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[17].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[22].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[26].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[30].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[33].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[35].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[38].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[39].setPrescribedMeshSizeAtVertex(paramaille_air)
vertex[40].setPrescribedMeshSizeAtVertex(paramaille_air)

## Physical Entites
myModel4.getRegionByTag(8).addPhysicalEntity(138)  # Air
myModel4.getRegionByTag(7).addPhysicalEntity(139)  # PML X
myModel4.getRegionByTag(4).addPhysicalEntity(141)  # PML Y
myModel4.getRegionByTag(6).addPhysicalEntity(142)  # PML Z
myModel4.getRegionByTag(2).addPhysicalEntity(140)  # PML XY
myModel4.getRegionByTag(3).addPhysicalEntity(145)  # PML XZ
myModel4.getRegionByTag(5).addPhysicalEntity(144)  # PML YZ
myModel4.getRegionByTag(1).addPhysicalEntity(143)  # PML XYZ

myModel4.getFaceByTag(40).addPhysicalEntity(146)   # XOZ
myModel4.getFaceByTag(32).addPhysicalEntity(146)   # XOZ
myModel4.getFaceByTag(25).addPhysicalEntity(146)   # XOZ
myModel4.getFaceByTag(35).addPhysicalEntity(146)   # XOZ

myModel4.getFaceByTag(39).addPhysicalEntity(147)   # YOZ
myModel4.getFaceByTag(22).addPhysicalEntity(147)   # YOZ
myModel4.getFaceByTag(14).addPhysicalEntity(147)   # YOZ
myModel4.getFaceByTag(30).addPhysicalEntity(147)   # YOZ

myModel4.getFaceByTag(8).addPhysicalEntity(149)    # XOY
myModel4.getFaceByTag(44).addPhysicalEntity(149)   # XOY
myModel4.getFaceByTag(36).addPhysicalEntity(149)   # XOY
myModel4.getFaceByTag(19).addPhysicalEntity(149)   # XOY

myModel4.getFaceByTag(41).addPhysicalEntity(148)   # Mirror

myModel4.getVertexByTag(1).addPhysicalEntity(1000000) # Dummy point for GetDP


## Mesh ##
##########
myModel4.mesh(3)
myModel4.save("haroche_mesh.msh")


## Display ##
#############
myModel4.setAsCurrent();
myModel4.setVisibility(1);
FlGui.instance()
FlGui.run()
FlGui.close()
