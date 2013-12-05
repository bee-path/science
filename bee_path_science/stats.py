# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:32:00 2013

@author: Oleguer Sagarra <osagarra@ub.edu>


Stats module:
    
    Takes care of statistics of flights and stops
"""

import numpy as np

def gyrad(points):
    """ Calculates center of masses and radius of gyration of the given points in UTM format (x,y)
    input:
        - array with (x,y) GPS updates
    output:
        - 2 tuple (x,y) value
        - gyration radius (float)
        - standard derivation of gyration radius distr. (aprox statistical error)
    """
    cm=np.mean(points,axis=0) #mean position
    scaledpoints=points-cm
    return cm,np.mean(map(np.linalg.norm,scaledpoints)),np.std(map(np.linalg.norm,scaledpoints))



#### Mean sq dsipl and so on ###




def r_sigma():#not implemented
    """ Calculates exponent gamma of sigma(r)~t^{gamma} of the given points in UTM format (x,y,t)
    input:
        - array with (x,y,t), which alternatively are the start and the end of flights & stops sorted in time
    output:
        - gamma value,goodness of fit value,real points,fitted func
    """
    raise NotImplementedError("Not yet implemented")
    return None
#    if not self.N_points()[-2]: # not flights
#        gamma,p_value,rt,x,y=0,0,np.array([0,0]),0,0
#    else:
#        try:
#            from collections import defaultdict
#        except ImportError:
#            print "Needs collections module "
#            return
#        rt=defaultdict(float)
#        rcount=defaultdict(float)
#        points=self.start(stops=False)
#        N=len(points)
#        for d in iter(range(1,N)):
#            for i in iter(range(0,N-d)):
#                r0=points[i][0:2]
#                r1=points[i+d][0:2]
#                dr2=np.linalg.norm(r1-r0)*np.linalg.norm(r1-r0)
#                dt=points[i+d][2]-points[i][2]
#                rt[dt]+=dr2
#                rcount[dt]+=1
#        for key,val in rt.items():
#            if val==0:
#                del rt[key]
#            else:
#                rt[key]=rt[key]/float(rcount[key])
#        rt=np.array([ (k,rt[k]) for k in sorted(rt.keys())])
#        rt_log=np.log10(rt[:int(len(rt)/div)])
#        if N>1:
#            gamma, intercept, r_value, p_value, std_err = stats.linregress(rt_log.T[0],rt_log.T[1])
#        else: # It did not work
#            print "Not enough data"
#            return 0,0,0,0,0
#        #for a,b in rt:
#        #    print a,b
#        #print gamma, intercept, r_value, p_value, std_err
#        if add:
#            self.add_atts('gamma',(gamma,p_value))
#        x=np.linspace(min(rt.T[0]),max(rt.T[0]))
#        y=x**gamma*10.**intercept
#    return gamma,p_value,rt,np.array([x,y])







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
    #mode={'UTM':['StartEnd_UTM','mean_point_UTM'],'latlon':['StartEnd_latlon','mean_point_latlon']}
    #from collections import defaultdict
    #spots_st=defaultdict(list)
    #spots_end=defaultdict(list)
    #for nst,stop in self.iter_stops(att=mode[unit][-1],count=True):
    #    for nfli,fli in self.iter_flights(att=mode[unit][0],count=True):
    #        if abs(np.linalg.norm(fli[0]-stop))<r:spots_st[nst].append(nfli)
    #    for nfli,fli in self.iter_flights(att=mode[unit][0],count=True):
    #        if abs(np.linalg.norm(fli[-1]-stop))<r:spots_end[nst].append(nfli)
    #return spots_st,spots_end
    #print "Feature not yet implemented"