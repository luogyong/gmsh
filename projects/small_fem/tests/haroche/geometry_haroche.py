#!/usr/bin/env python
from gmshpy import *
import numpy
import sys


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
if len(sys.argv) == 1:
    refinement_air  = 5.
    refinement_pml  = 5.
    refinement_mir  = 5.
    order           = 2.

    print('No mesh size given: using 5 5 5 2')

elif len(sys.argv) != 5:
    refinement_air = float(sys.argv[1])
    refinement_pml = float(sys.argv[1])
    refinement_mir = float(sys.argv[1])
    order          = 2.

    print('Using first mesh data for air, pml and mirror (order 2)')

else:
    refinement_air = float(sys.argv[1])
    refinement_pml = float(sys.argv[2])
    refinement_mir = float(sys.argv[3])
    order          = float(sys.argv[4])

paramaille_air  = lambda_haroche / refinement_air
paramaille_pml  = lambda_haroche / refinement_pml
paramaille_mir  = lambda_haroche / refinement_mir


## Options ##
#############
# 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
# 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D)

GmshSetOption('Mesh', 'Algorithm',    1.0)
GmshSetOption('Mesh', 'Algorithm3D',  1.0)
GmshSetOption('General', 'Verbosity', 4.0)
GmshSetOption('Mesh', 'ElementOrder', order)


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
myModel4.extrude(myModel4.getFaceByTag(9) , [0, 0, 0], [pml_x, 0, 0])
myModel4.extrude(myModel4.getFaceByTag(12), [0, 0, 0], [0, 0, pml_z])
myModel4.extrude(myModel4.getFaceByTag(8) , [0, 0, 0], [0, 0, pml_z])
myModel4.extrude(myModel4.getFaceByTag(23), [0, 0, 0], [0, pml_y, 0])
myModel4.extrude(myModel4.getFaceByTag(14), [0, 0, 0], [0, pml_y, 0])
myModel4.extrude(myModel4.getFaceByTag(25), [0, 0, 0], [0, pml_y, 0])
myModel4.extrude(myModel4.getFaceByTag(19), [0, 0, 0], [0, pml_y, 0])

myModel4.occconnect()

## Mesh size at Vertices
for i in range(1, myModel4.getNumVertices()):
    myModel4.getVertexByTag(i).setPrescribedMeshSizeAtVertex(paramaille_pml)

myModel4.getVertexByTag(28).setPrescribedMeshSizeAtVertex(paramaille_mir)
myModel4.getVertexByTag(29).setPrescribedMeshSizeAtVertex(paramaille_mir)
myModel4.getVertexByTag(33).setPrescribedMeshSizeAtVertex(paramaille_mir)

myModel4.getVertexByTag(1).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(3).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(5).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(7).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(9).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(11).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(18).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(27).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(30).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(31).setPrescribedMeshSizeAtVertex(paramaille_air)
myModel4.getVertexByTag(32).setPrescribedMeshSizeAtVertex(paramaille_air)


## Physical Entites ##
######################
myModel4.getRegionByTag(8).addPhysicalEntity(138)  # Air
myModel4.getRegionByTag(1).addPhysicalEntity(139)  # PML X
myModel4.getRegionByTag(4).addPhysicalEntity(141)  # PML Y
myModel4.getRegionByTag(2).addPhysicalEntity(142)  # PML Z
myModel4.getRegionByTag(6).addPhysicalEntity(140)  # PML XY
myModel4.getRegionByTag(3).addPhysicalEntity(145)  # PML XZ
myModel4.getRegionByTag(5).addPhysicalEntity(144)  # PML YZ
myModel4.getRegionByTag(7).addPhysicalEntity(143)  # PML XYZ

myModel4.getFaceByTag(1).addPhysicalEntity(146)    # XOZ
myModel4.getFaceByTag(10).addPhysicalEntity(146)   # XOZ
myModel4.getFaceByTag(14).addPhysicalEntity(146)   # XOZ
myModel4.getFaceByTag(35).addPhysicalEntity(146)   # XOZ

myModel4.getFaceByTag(7).addPhysicalEntity(147)    # YOZ
myModel4.getFaceByTag(20).addPhysicalEntity(147)   # YOZ
myModel4.getFaceByTag(24).addPhysicalEntity(147)   # YOZ
myModel4.getFaceByTag(34).addPhysicalEntity(147)   # YOZ

myModel4.getFaceByTag(4).addPhysicalEntity(149)    # XOY
myModel4.getFaceByTag(17).addPhysicalEntity(149)   # XOY
myModel4.getFaceByTag(28).addPhysicalEntity(149)   # XOY
myModel4.getFaceByTag(39).addPhysicalEntity(149)   # XOY

myModel4.getFaceByTag(36).addPhysicalEntity(148)   # Mirror
myModel4.getFaceByTag(37).addPhysicalEntity(148)   # Mirror
myModel4.getFaceByTag(38).addPhysicalEntity(148)   # Mirror

myModel4.getVertexByTag(1).addPhysicalEntity(1000000) # Dummy point for GetDP


## Mesh & Save ##
#################
brepName = 'haroche.brep'
meshName = 'haroche'                         + \
           '_geoorder_%d'    %order          + \
           '_air_%d'         %refinement_air + \
           '_pml_%d'         %refinement_pml + \
           '_mir_%d'         %refinement_mir + '.msh'

myModel4.mesh(3)
myModel4.save(meshName)
myModel4.save(brepName)


## Display ##
#############
print('Data used: ')
print('  ** Air   : ' + str(refinement_air))
print('  ** PML   : ' + str(refinement_pml))
print('  ** Mirror: ' + str(refinement_mir))
print('  ** Order : ' + str(order))

# myModel4.setAsCurrent();
# myModel4.setVisibility(1);
# FlGui.instance()
# FlGui.run()
# FlGui.close()
