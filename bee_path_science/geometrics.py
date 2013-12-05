# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:04:07 2013

@author: Oleguer Sagarra <osagarra@ub.edu>

This modules takes care of geometric transforms using projected points on x,y space
"""

# External modules #
import numpy as np



def angle_diff(points):
    """ given a set of vectors defined by delta_x,delta_y, get their angles """
    angles=np.arctan2(points.T[1],points.T[0])
    norms=np.hypot(points.T[1],points.T[0])
    ang_dif=np.empty(len(points)-1)
    for i in range(len(points)-1):
        dot=np.sum(points[i]*points[i+1])/(norms[i]*norms[i+1])
        ang_dif[i]=np.arccos(dot)
    return ang_dif



def distancepointtoline(p1,p2,p3,eps=0.000001):
    """Calculates the ortogonal distance from point p3=(x3,y3) to the line across
    the points p1=(x1,y1) and p2=(x2,y2)
    intput:
        - p1
        - p2
        - p3
        - minimum error tolerance
    output:
        - distance (scalar)
    """
    if np.linalg.norm(p1-p2)<eps:
        d=0
    else:
        a=(p2[1]-p1[1])/(p2[0]-p1[0])
        b=(p1[1]*p2[0]-p2[1]*p1[0])/(p2[0]-p1[0])
        sq=np.sqrt(a*a+1)
        d=abs(a*p3[0]-p3[1]+b)/sq
    return d
