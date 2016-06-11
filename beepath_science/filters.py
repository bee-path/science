# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:11:08 2013

@author: Oleguer Sagarra <osagarra@ub.edu>

This module includes filtering functions for sets of points (also known as Traces)

"""
## INternal modules ##
import constants as cnt


## External modules ##
import pyproj as pr
import numpy as np


##### Flight selections functions #####


def filter_v_updates(points,v_max=50/3.6,zone=cnt.zone,UTM=True):
    """ Splits trace into sets of points according to sconsecutive updates with impossible velocity jumps
        Given N updates, computes de N-1 inst. velocity updates V_n=(X_n+1-X_n)/(T_n+1-T_n), 
        then if V_n>V_max then the dataset gets separated 
        input: Points (3-col numpy array lat,lon,time)
               v_max (maximum velocity permitted)
               zone UTM zone definition
               mode: If UTM, then no conversion
        output: list of sets of points valids
    """
    if not UTM:
    # Convert points
        utm=cnt.utm()
        xy=np.array(utm(points.T[0],points.T[1])).T # transformacio geo (lon-lat)
        xy=xy-xy[0] # refer to arbitrary origin (to ease calculations)
    else:
        xy=points.T[:-1].T
    #deltar=np.array(map(np.linalg.norm,np.diff(xy.T).T),dtype=float)
    #deltat=np.array(np.diff(points.T[-1]),dtype=int) # time diff
    #v=(deltar)/(deltat+1) # avoid problems
    #group_points=[e for e in np.split(points,np.where(v>v_max)[0]) if len(e)>1] # only sets of points bigger than one (if not, means isolated bad point)
    lastp=points[0]
    gpoints=[]
    for p in points[1:]:
        deltar=np.linalg.norm(p[:-1]-lastp[:-1])
        deltat=p[-1]-lastp[-1] # time diff
        if deltat >0:
            v=deltar/deltat
            if v>v_max:
                #print "Warning, bad point! v=%f" % (v)
                pass
            else:
                lastp=p
                gpoints.append(p)
        else:
            #print "Warning, repeated time! deltat = %f \t deltar = %f" %(deltat,deltar)
            pass
    return [np.array(gpoints)]
