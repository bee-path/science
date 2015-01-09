# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:04:07 2013

@author: Oleguer Sagarra <osagarra@ub.edu> and Mario Gutierrez

This modules takes care of geometric transforms using projected points on x,y space
"""

# External modules #
import numpy as np
import shapely.geometry as sh
import shapely.affinity as aff


############# Box things #################
def pointsinbox(Points,R): #works fine
    """ From a group of points, check if they fit all in a box ""
    input: 
        - Group of points (x,y,s) [list]
        - box_coordinates (minx,miny,maxx,maxy) [4 item list]
    output:
        True or False
    Uses shapely
    """
    box=generate_box(Points[0],Points[-1],R)
    return all(box.contains(sh.Point(e)) or box.touches(sh.Point(e)) for e in Points[1:-1]) # exclude boundaries
    
def generate_box(p0,p1,R): # works fine
    """ Generates a box through the line joining p1 and p2, with width R """
    # Compute turning angle
    p0=np.array(p0)
    p1=np.array(p1)
    a=angle_vect(p0,p1) # in radians
    #print a
    r=np.linalg.norm(p1-p0) # distance
    p11=[p0[0],p0[1]-R]
    p12=[p0[0],p0[1]+R]
    p13=[p0[0]+r,p0[1]-R]
    p14=[p0[0]+r,p0[1]+R]
    box = sh.Polygon([p11,p12,p14,p13,p11])
    return aff.rotate(box,a,origin=p0,use_radians=True)
    
def positions_box(b,h,angle,origin=None):
    """ Calculates p0 and p1 to fit real data to simulation box 
    b= Base of the box
    h= Heigh of the box
    angle= Rotation of the box (in counter-clockwise is positive and Radians)
    origin= Origin of the UTM map, if set to none, use default
    """
    if not origin:
        import constants as cnt
        origin = cnt.origin
    origin=np.array(origin)
    R = h/2.
    p0 = origin + R*np.array([-np.sin(angle),np.cos(angle)])
    p1 = p0 + b*np.array([np.cos(angle),np.sin(angle)])
    return p0,p1

def crop_points_box(Points, b,h,angle,origin=None):
    """ Returns only points in the box given by origin, base, heigh and angle 
        Points: Array with x-y coordinates
        Angle: In radians
        b,h in arbitrary units
    """
    ps = positions_box(b,h,angle,origin)
    bb = generate_box(ps[0],ps[1],h/2.) # last param is dummy
    return np.array([e for e in Points if bb.contains(sh.Point(e[:2])) or bb.touches(sh.Point(e[:2])) ])
    
################# Geometric things #################


def radius2(p):
    """ Returns radius^2 of points (x-y) """
    return p.T[0]*p.T[0]+p.T[1]*p.T[1]



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
        a=abs(p3[1]*p1[0]-p3[1]*p2[0]-p1[1]*p3[0]+p2[1]*p3[0]-p2[1]*p1[0]+p1[1]*p2[0])
        b=np.sqrt((p1[1]-p2[1])*(p1[1]-p2[1])+(p1[0]-p2[0])*(p1[0]-p2[0])) #ok
        d=a/b
    return d


def rotate_xyt(points,angle):
    """ Rotates  input to x-y-t points by certain angle (in radians)
    input:
        - points x-y-t in utm
        - angle
    output:
        - rotated x-y-t
    """
    points = np.array(points)
    x = points[:,0] * np.cos(angle) - points[:,1] * np.sin(angle)
    y = points[:,0] * np.sin(angle) + points[:,1] * np.cos(angle)
    t = points[:,-1]
    
    return np.vstack((x,y,t)).T
    
def angle_vect(p0,p1):
    """ Given two points origin p0 and des p1, get their angle (uniquely defined) """
    p=np.array(p1)-np.array(p0)
    return np.arctan2(p[1],p[0])
