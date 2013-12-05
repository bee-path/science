# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:08:55 2013


@author: Oleguer Sagarra <osagarra@ub.edu>

Model module: Implements the models to process traces (see the docs for more info)
"""

# External modules #
import numpy as np
#from scipy import stats as sts


# Internal modules #
import classes as cls
import geometrics as geom
#import constants as cnt

######



### Stop or move ###

def stopormove(Points,R): # works fine!
    """Calculates which update corresponds to a stop or to a movement
    input:
        - x-y-t GPS UTM converted points
        - tolerance parameter (in meters)
    output:
        - x-y-t-s Points, where s is status: 0 for stoped, 1 for walking
    """
    if len(Points)>1:
        mov=np.ones((len(Points)),dtype=np.int)
        for i in iter(range(len(Points)-1)):
            r=Points[i][0:2]-Points[i+1][0:2]
            d=np.linalg.norm(r)
            #if (d-R)<eps:
            if d<R:
                mov[i]=0
        mov[-1]=0 # last point is a stop by default
        mov[0]=0 # first point is a stop by default
    else:
        mov=[0]
    return np.hstack((Points,np.reshape(mov,(len(Points),1))))

def stopormove_v(Points,r_min,t_min,t_max):
    """Calculates which update corresponds to a stop or to a movement
    input:
        - x-y-t GPS UTM converted points
        - tolerance parameter (in meters)
    output:
        - x-y-t-s Points, where s is status: 0 for stoped, 1 for walking
    """
    v_min=r_min*1./t_min
    if len(Points)>1:
        mov=np.ones((len(Points)),dtype=np.int)
        for i in iter(range(len(Points)-1)):
            r=Points[i][0:2]-Points[i+1][0:2]
            t=Points[i+1][2]-Points[i][2]
            d=np.linalg.norm(r)
            if t<t_min: #not enough time or too much time
                print "# Too short or too long, t=%f" % (t)
                if i!=0:
                    mov[i]=mov[i-1] # copy last
                else:
                    mov[i]=0
            elif t>t_max: mov[i]=1; # if too long, moving as default (malfoctioning GPS)
            else:
                if d*1./t<v_min: # not enough velocity, but enough time
                    if d*1./t>0.3: print "# Too slow, v=%f" % (d*1./t)                
                    mov[i]=0    
        mov[-1]=0 # last point is a stop by default
        mov[0]=0 # first point is a stop by default
    else:
        mov=[0]
    return np.hstack((Points,np.reshape(mov,(len(Points),1))))


### Split flights ###


def separe_flights(oTrace):
    """from a trace gets a list of lists, where each list is a set of consecutive flights between stops"""
    alls=np.sort(list(oTrace.iter_all()))
    separe=[]
    i=0
    for e in alls:
        if isinstance(e,cls.Stop):
            separe.append(i)
        i+=1
    fl_group=[]
    # group flights
    for g in np.split(alls,separe): 
        try:
            fl_group.append(g[1:])
        except IndexError: pass
    return [e for e in fl_group if any(e)]


def flightworks(Points,i,j,R): #works fine
    """ From a group of points, select minimal group that corresponds a flight
    input: 
        - Group of points (x,y,s) [list]
        - Starting index of flight
        - Ending index of flight
        - Tolerance parameter
    output:
        - Index of minimal group of updates that correspond to a flight -> points[i:j+1]
    """
    # Possibly here add box condition !
    flag=False
    while j>i+1 and not flag:
        flag=True
        for k in iter(range(i+1,j)):
            d=geom.distancepointtoline(Points[i][0:2],Points[j][0:2],Points[k][0:2])
            if d>R:
                flag=False
                j-=1
                break
    return j





### Models ####
    
def rectangularmodel(Points,R):
    """Split the trace in flights according to rectangular model (see aditional docs)
    input:
        - UTM points x-y-t-s as obtained from stopormove()
        - R parameter tolerance
    output:
        - Track Class with Flights and Stops
    """
    Track=cls.Trace() # create Trace
    dummy=[] # create dummy list
    #if Points[0][3]==1: print "first a flight?" # first is a flight (this is impossible) --> checked already
    for i in iter(range(len(Points)-1)): # For all points except the last one (since we check forward)
        point=Points[i][0:3] # store point
        dummy.append(point) # add to dummy list
        if Points[i][3]!=Points[i+1][3]: #if change (last point has to be a stop by default (see else at the bottom)
            if Points[i][3]==0: # stop -> add stop from dummy array
                s=cls.Stop(np.array(dummy))
                Track.add_stop(s)
                dummy=[point] # keep last point () (to add to next flight)
                #raw_input('added stop %s \n index: %i' % (s,i))
            else: # flight
                dummy.append(Points[i+1][0:3]) # append next stop [end point, flight, next-stop point]
                #if len(dummy)<3: print "Warning!!" #(this is impossible, already checked!)
                if len(dummy)==3: # means just one update with flight [end point, flight] or simply [flight]
                    f=cls.Flight(np.array(dummy))
                    Track.add_flight(f)
                    #raw_input('added single flight %s \n index: %i' % (f,i))
                else: # more than 1 update
                    aux_fl=[] #list of flights to add                   
                    while len(dummy)>2: # means not only one or two update: detect flights
                        j=len(dummy)-1 # last point
                        k=flightworks(dummy,0,j,R) # get minimum flight from dummy                        
                        aux_fl.append(cls.Flight(np.array(dummy[0:k+1]))) # append flight
                        dummy=dummy[k:] # remove points from dummy list except last point
                    if len(dummy)==2: # two updates surviving: last flight of the set
                        aux_fl.append(cls.Flight(np.array(dummy)))
                    #if len(dummy)==1: pass  #one update survives by construction
                    for fl in aux_fl: Track.add_flight(fl) #add flights
                dummy=[] # empty list
        else: # nothing to do except if last point
            if i==len(Points)-2: # two last points are stops 
                dummy.append(Points[i+1][0:3]) # add last point to dummy
                s=cls.Stop(np.array(dummy)) # add stop
                Track.add_stop(s)
    return Track



def pausemodel(oTrace):
    """ merges successive flights from the rectangular model if their relative angle does not exceed theta_max
    input: theta_max, Trace containing stops and flights
    output: new Trace containing merged flights and same stops
    """
    new_Trace=cls.Trace()
    # add stops
    for s in oTrace.iter_stops():
        new_Trace.add_stop(s) 
    fl_group=separe_flights(oTrace)
    for fl in fl_group: new_Trace.add_flight(np.sum(fl))
    return new_Trace


        
def angularmodel(oTrace,theta_max):
    """ merges successive flights from the rectangular model if their relative angle does not exceed theta_max
    input: theta_max, Trace containing stops and flights
    output: new Trace containing merged flights and same stops
    """
    new_Trace=cls.Trace()
    # add stops
    for s in oTrace.iter_stops():
        new_Trace.add_stop(s) 
    fl_group=separe_flights(oTrace)
    # and now merge if needed
    for gr in fl_group:
        if len(gr)>1: # if more than 1 flight
            r=np.array([f.UTM[-1][:-1]-f.UTM[0][:-1] for f in gr])
            dif_angle=geom.angle_diff(r)
            #dif_angle=np.abs(np.diff(np.arctan2(r.T[1],r.T[0]))) # not really working
            #print np.rad2deg(dif_angle)
            if any(dif_angle<0): print dif_angle
            if len(np.where(dif_angle>np.deg2rad(theta_max))[0])!=0: # if something to merge
                new_fl=np.split(gr,np.where(dif_angle>np.deg2rad(theta_max))[0]+1)
                for fl in new_fl: 
                    if len(fl)!=0: new_Trace.add_flight(np.sum(fl))
            else: new_Trace.add_flight(np.sum(gr)) #if not, merge all!
        else: 
            if len(gr)!=0: new_Trace.add_flight(gr[0])
    return new_Trace





