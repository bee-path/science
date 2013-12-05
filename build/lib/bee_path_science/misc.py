# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:04:22 2013

@author: Oleguer Sagarra <osagarra@ub.edu>

Misc module with handles for code analysis
"""
import numpy as np


def unique(a):
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1) 
    return a[ui]
    

def type_move(Trace,params={'v':(4,5),'r_gyr':(75)},mov_types={(1,1):'shark',(1,-1):'fox',(0,1):'dollar',(0,-1):'bee',(-1,1):'pollen',(-1,-1):'cell'},rand=False):
    """ Gets a type of animanl movement
    input:
        - params: set threshold for calssif (velocity and r_gyr) in (km/h and m)
        - mov_types: types of movement
        - rand: If true, make up the type (in case something fails)
    output:
        - type_mov: a string selected according to pre-set parametters dic
        - status: If rand True, else false
    Warning: This was only used for a given event, results are just for divulgation purposes
    """

    if rand:
        return mov_types.values().pop(np.random.randint(0,high=len(mov_types))),False
    else:
        v=Trace.v(kmh=True)[0]
        rgyr=Trace.rgyr()[0][0]
        if len(params)==2: #then v and r_gyr
            if v<params['v'][0]:
                v_stat=-1
            elif v>=params['v'][0] and v<params['v'][1]:
                v_stat=0
            else:
                v_stat=1
            if rgyr < params['r_gyr']:
                r_stat=-1
            else:
                r_stat=1
            return mov_types[v_stat,r_stat],False
        else:
            print "Not enough params or too many params to classify, returning rand"
            return mov_types.values().pop(np.random.randint(0,high=len(mov_types))),True