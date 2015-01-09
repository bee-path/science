# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:05:56 2013

@author: Oleguer Sagarra <osagarra@ub.edu>

This modules helps in quick mapping for results, requires "simplekml module" to work
"""
## Internal modules
#import stats as st
## External modules
import pylab as pl
#import numpy as np

def gmaps(Trace,path,details=True,stops=False,det_points=False):
    """ Creates google map with given trace indicating updates w real points, flights with arrows
    if details is set to False, then only arrows
    input:
        - Trace (class Trace)
        - path (string): path to save kml file
        - Stops (boolean): If True, plot also stops
    Note:
        - Needs simplekml to work
    """
    try:
        import simplekml as spk
    except ImportError:
        print "Please install simplekml module"
        return
    hubmap=spk.Kml() # We create the map
    for num,fli in enumerate(Trace.flights): #for each flight
        if details:
            lin=hubmap.newlinestring(name="Flight %i"%(num), description="Real flight %i" % (num),coords=fli.lat_lon('UTM').T[::-1].T.tolist())
            lin.style.linestyle.color = 'ffff0000'  # Blue
            lin.style.linestyle.width= 3  # 10 pixels
            if det_points:
                for i,point in enumerate(fli.lat_lon('UTM')):  # for each point
                    hubmap.newpoint(name="Flight %i, update %i" % (num,i), coords=[point[::-1]])
        lin=hubmap.newlinestring(name="Flight %i"%(num), description="Filtered flight %i" % (num),coords=fli.lat_lon('StartEnd').T[::-1].T.tolist())
        lin.style.linestyle.color = 'ff00000f'  # Red
        lin.style.linestyle.width= 10  # 10 pixels
    if stops:
        for num,fli in enumerate(Trace.stops): #for each flight
            if details:
                lin=hubmap.newlinestring(name="Stop %i"%(num), description="Real Stop %i" % (num),coords=fli.lat_lon('UTM').T[::-1].T.tolist())
                lin.style.linestyle.color = 'fffff001'  # Blue
                lin.style.linestyle.width= 3  # 10 pixels
                if det_points:
                    for i,point in enumerate(fli.lat_lon):  # for each point
                        hubmap.newpoint(name="Stop %i, update %i" % (num,i), coords=[point[::-1]])
            hubmap.newpoint(name="Mean Stop %i" % (num), coords=[fli.lat_lon('mean_point')[::-1]])
    hubmap.save(path)
    return



def quickview(points,tracs,path=None):
    """ generates quick plot of points and flights 
    Input:
        Points: Raw points in x,y
        Tracs: LIst of 3 model tracks [Rect, Angular, Pause]
        Path: Path to save image
    
    """
    styles=('-','--',':')
    names=('rect','angular','pause')
    colors=('y','r','g')
    pl.clf()
    pl.plot(points.T[0],points.T[1],'o',label='raw_points')
    i=0
    for tr in tracs:
        for fl in tr.iter_flights():
            pl.plot(fl.StartEnd.T[0],fl.StartEnd.T[1],styles[i],label=names[i],color=colors[i],lw=3.)
        if names[i]=='pause':
            for fl in tr.iter_stops():
                pl.plot(fl.lat_lon('mean_point').T[0],fl.lat_lon('mean_point').T[1],'x',label=names[i],color=colors[i],markersize=12)
        i+=1
    pl.xlabel('x UTM')
    pl.ylabel('y UTM')
    #pl.legend(loc='best')
    if path: pl.savefig(path)
    else: pl.show()



#def map_stats(Trace,i): # not finished
#    """ Maps trace stats """
#    pp=Trace.allpoints(stops=True,time=True)
#    pp=pp-pp[0] # center time
#    r=st.r_sigma(pp,1)[:-2] # r_sigma
#    rg=st.gyrad(pp,1) # radius of gyration
#    
#    #pl.clf()
#    fig=pl.figure(i)
#    ax=pl.subplot(2,2,1)
#    ax.plot(pp.T[0],pp.T[1],'.-',label='pos')
#    ax.plot(rg.T[1],rg.T[2],'.-',label='CM')
#    ax.set_title('Position')
#    ax.set_xlabel('x [UTM]')
#    ax.set_ylabel('y [UTM]')
#    ax.set_title('Position and CM position')    
#    ax.legend(loc='best')
#
#
#    ax=pl.subplot(2,2,2)    
#    #ax.plot(r.T[0],r.T[6],'.-',label='sigma R increment')
#    ax.plot(r.T[0],r.T[5],'.-',label='sigma x')
#    ax.plot(r.T[0],r.T[4],'.-',label='sigma y')
#    ax.plot(r.T[0],np.mean(r.T[4][:10]/r.T[0][:10])*r.T[0],'--',label=None)
#    ax.set_xscale('log')
#    ax.set_yscale('log')
#    ax.set_xlabel('t')
#    ax.set_ylabel(r'\sigma^2_r')
#    ax.set_title('Diffusion Coeff')    
#    ax.legend(loc='best')    
##    ax.legend(bbox_to_anchor=(0.5,0))    
#
#    
#    ax=pl.subplot(2,2,4)
#    ax.plot(r.T[0],r.T[3],'.-',label='R increment')
#    ax.plot(r.T[0],r.T[1],'.-',label='x increment')
#    ax.plot(r.T[0],r.T[2],'.-',label='y increment')
#    ax.set_xlabel('t')
#    ax.set_ylabel(r'<r>')
#    ax.set_title('Average Increment Position')    
#    ax.set_xlim(xmin=0,xmax=1000)
#    ax.legend(loc='best')
#    
#    
#    ax=pl.subplot(2,2,3)
#    ax.plot(rg.T[0],rg.T[-2],'.-',label='rgyr')
#    ax.plot(rg.T[0],np.sqrt(rg.T[1]*rg.T[1]+rg.T[2]*rg.T[2]),'.-',label='cm')
#    ax.set_xlabel('t')
#    ax.set_ylabel('r_gyr')
#    ax.set_title('Radius of gyration')
#    ax.legend(loc='best')
#
#    return fig
