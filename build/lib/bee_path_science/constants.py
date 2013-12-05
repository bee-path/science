# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 09:49:05 2013

@author: Oleguer Sagarra <osagarra@ub.edu>


This module defines constants and variables for the package
"""
# External modules
import pyproj as pr

# Internal modules



## Geodesic data ###
zone=31 # zone utm 31T
origin=(431769,4582211) # map origin #Coordenades geodèsiqueS: N41º23.312 E02º11.034
ellips='WGS84'
proj='utm'


def utm():
    return pr.Proj(ellips=ellips, proj=proj, zone=zone, preserve_units=False) # Not preserving units (meters)



## Other ##
eps=0.0000001 # very small value (to avoid problems)




