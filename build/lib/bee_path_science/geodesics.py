# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 09:59:43 2013

@author: Oleguer Sagarra <osagarra@ub.edu>

This module takes care of geodesic transforms and projections, mainly uses UTM projection (2-D)
"""

# External modules
import pyproj as pr
import numpy as np

# Internal modules
import constants as cnt

utm=cnt.utm()

##### Geodesic functions ######

def UTM_2_latlon(points,origin,zone):
    """ Converts (x-y) input to lat-lon files on given UTM zone with given origin
    input:
        - origin must be a 2 index tupple (x,y)
        - np.array with at least 2 columns x-y
        - Scalar UTM zone
    output:
        - array with lat,lon
    """
    points=np.array(points)
    coords=utm(points.T[0]+origin[0],points.T[1]+origin[1],inverse=True)
    return np.array(coords)[::-1].T # transformacio geo (lat-lon)

def latlon_2_UTM(points,origin,zone,time=False):
    """ Converts lon-lat-time input to x-y-files on given UTM zone with given origin
    input:
        - origin must be a 2 index tupple (x,y)
        - np.array with at least 2 columns lat-lon
        - Scalar UTM zone
    output:
        - array with lat,lon
    """
    points=np.array(points)
    xy=np.array(utm(points.T[0],points.T[1])).T-origin # transformacio geo (lon-lat)
    if time:
        return np.array([xy.T[0],xy.T[1],points.T[-1]]).T # transformacio geo (lon-lat)
    else:
        return np.array(xy) # transformacio geo (lon-lat)
        