# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:32:00 2013

@author: Oleguer Sagarra <osagarra@ub.edu>


Stats module:
    
    Takes care of statistics of flights and stops
"""
# External Modules
import numpy as np
#try:
#    import pandas as pd
#except ImportError:
#    print "Install pandas or module stats will not work entirely"


#Internal Modules
import geometrics as geom
##################################
##################################
#### Point based statistics ######
##################################
##################################

def gyrad(points,time=False): #unchecked
    """ Calculates center of masses and radius of gyration of the given points in UTM format (x,y)
    As defined in [1].
    input:
        - array with (x,y) GPS updates
        - time: if time, plots evolution in time
    output:
        - 2 tuple (x,y) value
        - gyration radius (float)
        - standard derivation of gyration radius distr. (aprox statistical error)
    if time: Output is t,cm[x],cm[y],radius,std_radius
    
    [1] M. C. Gonzalez, C. a C. A. Hidalgo, A. L. A.-L. Barabási, M. C. González, C. A. H. 
    & A.-L. B. Marta C. González, and M. C. Gonz, Nature 453, 779 (2008).
    """
    if len(points[0])== 3 and not time:
        points=points.T[:-1].T
    if time:
        t=points.T[-1]
        points=points.T[:-1].T
        cm=np.array([np.mean(points[:e],axis=0) for e in xrange(1,len(points)+1)])
        scaledpoints=points-cm[-1] # last CM
        rr=geom.radius2(scaledpoints)
        r2=np.cumsum(rr)/np.arange(1,len(rr)+1)
        return np.array([t,cm.T[0],cm.T[1],r2]).T
    else:
        cm=np.mean(points,axis=0) #mean position
        scaledpoints=points-cm
        return cm,np.sqrt(np.mean(geom.radius2(scaledpoints)))


#### Mean sq dsipl and so on ### -> To compute over aggregated traces!

def r_delta_iter(p,min_lag=1, max_lag=10000): #unchecked
    """ Calculates position increments with moving lag
        Input: 
            Set of points (x-y-t)
            Minimum lag: Minimium point window (in updates, not time!)
            Max lag: Maximum point window
        Output: 
            Unordered iterator (lag,deltax,deltay)
    Important: Computes lags as:
            Delta_r (tau)  = r(t+tau)- r(t)
    Over all possible ranges of t.
    """
    i=0
    lag=min_lag
    lags=np.zeros(3)
    while True:
        q=len(p[lag:])
        if q<2 or lag>max_lag:
            break
        dp = p[lag:]-p[:-lag]
        lags=np.vstack([lags,[(e[-1],e[0],e[1]) for e in dp]])
        #print r[i]
        lag+=1
        i+=1
    return lags




def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]




def tortuosity(): #not implemented
    """ Computes tortuosity of the trace
    output:
        Tortuosity= < 1 - ( direct distance between stops / real distance between stops ) >
    averaged over all stops
    input:
        - radius around stop points (tolerance)
        - 'UTM' or 'latlon' (depend on the unit of the geodesic data)
    output:
        - float between 0 and 1
    """
    raise NotImplementedError("Not yet implemented")
    return 0
    #mode={'UTM':['StartEnd','mean_point'],'latlon':['StartEnd','mean_point']}
    #from collections import defaultdict
    #spots_st=defaultdict(list)
    #spots_end=defaultdict(list)
    #for nst,stop in self.iter_stops(att=mode[unit][-1],count=True):
    #    for nfli,fli in self.iter_flights(att=mode[unit][0],count=True):
    #        if abs(np.linalg.norm(fli[0]-stop))<r:spots_st[nst].append(nfli)
    #    for nfli,fli in self.iter_flights(att=mode[unit][0],count=True):
    #        if abs(np.linalg.norm(fli[-1]-stop))<r:spots_end[nst].append(nfli)
    #return spots_st,spots_end
    
    
    
##################################
##################################
#### Trace based statistics ######
##################################
##################################

#### Statistical Functions  -> han de cridar funcions totals ####

def data_PDF(x, hbin=10, hmin=0, hmax=100):
    """ Returns an histogram 
    input:
        - Number of bins
        - Minimum of the range
        - Maximum of the range
    output:
        - Array with counts
        - Array with bins
    """
    hist, bins = np.histogram(x, hbin, [hmin,hmax])
    bins = bins[:-1] + (bins[1:]-bins[:-1])/2.
    pdf = hist/float(max(np.cumsum(hist)))
    err = np.sqrt(pdf*(1.-pdf)/float(max(np.cumsum(hist))))
    return np.vstack((bins, pdf, err)).T

    
def data_CCDF(x, hbin=10, hmin=0, hmax=100):
    """ Returns Cumulative Complementary Distribution Function 
     input:
        - Number of bins
        - Minimum of the range
        - Maximum of the range
    output:
        - Array with bins, pdf, error
    """
    hist, bins = np.histogram(x, hbin - 1, [hmin,hmax])
    hist = np.insert(hist, 0, 0)
    pdf = hist/float(max(np.cumsum(hist)))
    ccdf = 1. - np.cumsum(pdf)
    err = np.sqrt(ccdf*(1. - ccdf)/float(max(np.cumsum(hist))))
    return np.vstack((bins, ccdf, err)).T

