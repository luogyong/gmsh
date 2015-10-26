#!/usr/bin/env python
from gmshpy import *
import numpy
import sys
import os


## Pyhsics ##
#############
## Scaling
mm = 1.e-3
nm = 1.e-9

## Speed of light
epsilon0 = 8.854187817e-3  * nm
mu0      = 400. * numpy.pi * nm
c        = 1.0 / (numpy.sqrt(epsilon0 * mu0))

## Haroche Frequency, Wavelength & Wavenumber
freq_haroche   = 51.099e9
lambda_haroche = c / freq_haroche
k_haroche      = freq_haroche * 2 * numpy.pi / c


## Geomtrical Parameters ##
###########################
r                     = 39.4   * mm
R                     = 40.6   * mm
L_cav                 = 27.57  * mm
deltaz                =  5.    * mm
radius_mirror         = 25.    * mm
thick_mirror_atcenter =  1.415 * mm

## Geomtrical Model ##
######################
myModel11 = GModel()
myModel21 = GModel()
myModel31 = GModel()

myModel12 = GModel()
myModel22 = GModel()
myModel32 = GModel()

myModel11.addCylinder([0, 0, -R + L_cav / 2.], \
                      [0, 0, +L_cav / 2. + thick_mirror_atcenter], \
                      radius_mirror)

myModel12.addCylinder([0, 0, +R - L_cav / 2.], \
                      [0, 0, -L_cav / 2. - thick_mirror_atcenter], \
                      radius_mirror)

## Torus - almost spherical - centered in 0,0,0
myModel21.addTorus([0, 0, -R + L_cav / 2.], [0, 1, -R + L_cav / 2.], R-r, r)
myModel22.addTorus([0, 0, +R - L_cav / 2.], [0, 1, +R - L_cav / 2.], R-r, r)

myModel31.addSphere(0, 0, -R + L_cav / 2., r+r / 1000)
myModel32.addSphere(0, 0, +R - L_cav / 2., r+r / 1000)

myModel11.computeBooleanDifference(myModel21);
myModel11.computeBooleanDifference(myModel31);
myModel12.computeBooleanDifference(myModel22);
myModel12.computeBooleanDifference(myModel32);

myModel11.computeBooleanUnion(myModel12);

## Save ##
##########
brepName = 'haroche.brep'
myModel11.save(brepName)
