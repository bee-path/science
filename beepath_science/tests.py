# -*- coding: utf-8 -*-
"""

@author: Oleguer Sagarra <osagarra@ub.edu>

MOdule with tools to tests the models
"""
import numpy as np
import models
import classes


def generate_RW_trace(T,beta,Delta_t=0.1,tau=15,Rs=8,Rf=8):
    """
    Genrates a TRACE of a RW
        Inputs:
            T: Total simulation time
            Beta: Temperature of the walk (see generate_RW_points docs)
            Delta_t: Simulation time-step [default=0.1s]
            tau: Sampling ratio        [default = 15s]
            Rs: Stop criteria
            Rf: Flight criteria
        Output:
            Class Trace
    """
    # generate points
    Npoints = np.floor(T/Delta_t)
    xyt = generate_RW_points(beta,Npoints,Delta_t)
    # sampple
    xytt = filter_RW_points(xyt,tau)
    # Filter
    p = models.stopormove(xytt,Rs)
    Tr = models.rectangularmodel(p,Rf,Rs)
    return Tr
    
    
def generate_RW_points(beta,T = 3600., Delta_t = 0.1):
    """
    Generates a set of points generated with a bidimensional gaussian distribution
    with sigma^2: Delta_t/beta
            
    Inputs:
        beta: Parametter for the simulation of the RW
        T : TOtal time to simulate
        Delta_t : Interval of generation for time updates
        t0 : Initial time
    Outputs:
        3-item array (x-y-t) of length Npoints starting at time t0 in intervals.
    
    """
    Npoints = np.floor(T*1./Delta_t)
    sigma = np.sqrt(Delta_t*1./beta)
    x = np.cumsum(np.random.normal(loc=0,scale=sigma,size=Npoints))
    y = np.cumsum(np.random.normal(loc=0,scale=sigma,size=Npoints))
    t = np.arange(Npoints)*Delta_t
    return np.array([x,y,t]).T
    
def filter_RW_points(points,tau=15.):
    """
    Samples RW trajectory to intervals of tau seconds:
    Inputs:
        Points: 3-item array (x-y-t)
        Tau: Time interval to sample
    Outputs:
        3-item array with sampled points
    """
    # Get points to 0 #
    tfake = points.T[2]
    tfake = tfake - tfake[0]
    # Apply sampling
    tfake = np.floor(tfake*1./tau)
    inds = np.unique(tfake,return_index=True)[-1]
    return points[inds]
    
    
### Statistical function tests with analytical null model ###

def null_model_stops(Rmin=1,Rmax = 20, Np = 20, beta=1., T=3600., simulate=True, **kwargs):
    """
        Counts frac of stops for given params by kwargs under RW model
       Input:
            Rmin: Rs minimum to start computing
            Rmax: Rs max to end computing
            Np : Number of points in graph
            simulate : If True, returns simulation, else, returns theoretical values
            kwargs:
                beta, Delta_t, tau, samples
        Output:
            4-item tuple Rs,Frac_stops,avg_stops,std_stops
    """
    tau     = kwargs.get('tau',15.)
    Rf      = kwargs.get('Rf',8.)
    Delta_t = kwargs.get('Delta_t',0.1)
    frac    = np.zeros(Np)
    s_avg     = np.zeros(Np)
    #s_std     = np.zeros(Np)
    RRs     = np.linspace(Rmin,Rmax,Np)
    
    if simulate:
        for i,RR in enumerate(RRs):
            #dum_p = np.zeros(samples)
            #dum_s = np.zeros(samples)
            #dum_st = []
            #for ii,s in enumerate(range(samples)):
                # generate points
            xyt = generate_RW_points(beta,T,Delta_t)
            xytt = filter_RW_points(xyt,tau)
            p = models.stopormove(xytt,RR)
            Tr = models.rectangularmodel(p,Rf,RR)
            info = Tr.N_points(True)
            allt = [e for e in Tr.iter_stops('Delta_t')] 
            #ll = [e[np.where(e.T[-1]==0)] for e in np.split(p,np.where(np.ediff1d(p.T[-1])==-1)[0])]
            #allt = np.array([e.T[-2][-1]-e.T[-2][0] for e in ll]) + 15
               # dum_p[ii] = info[2]+info[3] 
                #dum_s[ii] = info[3]
                #dum_st = np.r_[dum_st,allt]
            #frac[i] = np.sum(dum_s)*1./np.sum(dum_p)
            frac[i] = len(np.where(p.T[-1]==0)[0])*1./len(p)
            s_avg[i]  = np.mean(allt)
            #s_std[i]  = np.std(allt)
    else:
        Np2 = np.floor(T/tau)
        p1 = p_n_stops(RRs, beta, tau, Np2)
        p2 = av_std_stop_length(RRs, beta, tau)
        frac = p1.T[1]
        s_avg = p2.T[1]
        #s_std = p2.T[2]
    return np.array([RRs,frac,s_avg]).T
        
### Aux functions: Stops ###

def prob(Rs, beta, tau):
    """ Computes stop-or-move probability """
    R20 = tau/beta
    return np.exp(-Rs * Rs / (2*R20))
### Stops ####
def p_n_stops(Rs, beta, tau, Ntot):
    """ Computes expected fraction of stops over total and std """
    p = prob(Rs, beta ,tau)
    avg =  (Ntot-1.)*1./Ntot*(1.-p)+1./Ntot
    std = np.sqrt((Ntot-1.)*1./Ntot*p*(1.-p))
    return np.array([Rs,avg,std]).T

def av_std_stop_length(Rs, beta , tau):
    """ Computes expected stop length and std """
    p = prob(Rs,beta,tau)
    avg = tau*(1.-p)/p
    std = np.sqrt((1-p)/(p*p))*tau
    return np.array([Rs,avg,std]).T

### Aux functions: Flights ###
from scipy import special as sp

def CDF_flight_all(Rs,Rf,beta,tau=15,hmin=1,hmax=20,Np=1000,pmin=1e-10):
    """ Computes unnormalized CDF of flight lengths """
    if hmin<Rs:
        hmin = Rs
    hh = np.linspace(hmin,hmax,Np)
    CDF = np.zeros(Np)
    for i,h in enumerate(hh):
        CDF[i] = CDF_flight_L(h,Rs,Rf,beta,tau,pmin)
    return np.array([hh,CDF]).T
        
    

def CDF_flight_L(L,Rs,Rf,beta,tau=15,pmin=1e-10):
    """
        Computes flight length CDF theoretical, unnormalized
    """
    if L<Rs: return 0
    p = prob(Rs,beta,tau)
    R0 = np.sqrt(tau/beta)
    suma = 0
    Nf = 1
    In=100
    while In>pmin:
        In = I_box(L,Rf,R0,Nf)*(p/2.)**(Nf-1)
        #print I_box(L,Rf,R0,Nf)
        #print In
        suma += In
        Nf+=1
    return suma

def I_box(hh,Rf,R0,Nf):
    N = np.sqrt(np.arange(Nf))
    I = sp.erf(hh/(2*R0*N))*sp.erf(Rf/(2*R0*N))
    return np.prod(I)
    
    
    