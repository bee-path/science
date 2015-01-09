# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:08:55 2013


@author: Oleguer Sagarra <osagarra@ub.edu> and Mario Gutiérrez

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

def stopormove(Points,R):
    """Calculates which update corresponds to a stop or to a movement
    input:
        - x-y-t GPS UTM converted points
        - tolerance parameter (in meters)
    output:
        - x-y-t-s Points, where s is status: 0 for stoped, 1 for walking
    Warning: Performs a stopormove conditoin for point i based on distance to point i+1
    """
    Points=np.array(Points)
    Points = sorted(Point,key = lambda x:x[-1]) # make sure points are in decreasing order in time!
    if len(Points) > 1:
        mov = np.ones (len (Points), dtype = int)
#        for i in iter (range (len (Points) - 2)): # all points e
        for i,p in enumerate(Points[:-1]): # all points (except last)
            r = Points[i+1][0:2] - p[0:2] # delta_x,delta_y
            d = np.linalg.norm (r) # norm
            if d < np.float(R):
                mov[i+1] = 0 # update
        if mov[-1] == 1: # if last point moving
            pnew= Points[-1]
            mov=np.append(mov,[0])
            Points2=np.append(Points,pnew)
            Points=np.reshape(Points2,(len(Points)+1,3))
        mov[0] = 0 # first point is a stop by default
    else:
        mov = [0]
    return np.hstack ((Points, np.reshape(mov, (len (Points), 1))))

def stopormove_mean(Points,R): # Experimental stuff
    # Needs to be rewritten (forward stop or move)
    """Calculates which update corresponds to a stop or to a movement
    input:
        - x-y-t GPS UTM converted points
        - tolerance parameter (in meters)
    output:
        - x-y-t-s Points, where s is status: 0 for stoped, 1 for walking
    """
    Points=np.array(Points)
    Points = sorted(Point,key = lambda x:x[2]) # make sure points are in decreasing order in time! 
    if len(Points) > 1:
        mov = np.ones((len (Points)), dtype = np.int)
        st = np.array([Points[0][0:2]])
        for i in iter (range (1, len (Points))):
            c = np.average (st, axis=0)
            r = Points[i][0:2] - c
            d = np.linalg.norm (r)
            if d < np.float(R):
                mov[i] = 0
                st = np.append (st, [Points[i][0:2]], axis=0)
            else:
                st = np.array([Points[i][0:2]])
        mov[-1] = 0 # last point is a stop by default
        mov[0]  = 0 # first point is a stop by default
    else:
        mov=[0]
    return np.hstack((Points,np.reshape(mov,(len(Points),1))))

def stopormove_log(Points,R): # Experimental stuff
    # Needs to be rewritten (forward stop or move)
    """Calculates which update corresponds to a stop or to a movement
    input:
        - x-y-t GPS UTM converted points
        - tolerance parameter (in meters)
    output:
        - x-y-t-s Points, where s is status: 0 for stoped, 1 for walking
    """
    Points=np.array(Points)
    Points = sorted(Point,key = lambda x:x[2]) # make sure points are in decreasing order in time! 
    if len(Points)>1:
        mov = np.ones((len(Points)),dtype=np.int)
        st = np.array([Points[0][0:2]])
        for i in iter (range (1, len (Points))):
            c = np.average (st, axis=0)
            r = Points[i][0:2] - c
            d = np.linalg.norm (r)
            if d < (R * (1. + np.log (len(st)))):
                mov[i] = 0
                st = np.append (st, [Points[i][0:2]], axis=0)
            else:
                st = np.array([Points[i][0:2]])
        mov[-1]=0 # last point is a stop by default
        mov[0]=0 # first point is a stop by default
    else:
        mov=[0]
    return np.hstack((Points,np.reshape(mov,(len(Points),1))))







def stopormove_v(Points,r_min,t_min,t_max, verbose=False):
    # Needs to be rewritten (forward stop or move)
    """Calculates which update corresponds to a stop or to a movement
    input:
        - x-y-t GPS UTM converted points
        - tolerance parameter (in meters)
        If verbose: Prints verbose
    output:
        - x-y-t-s Points, where s is status: 0 for stoped, 1 for walking
    """
    Points=np.array(Points)
    Points = sorted(Point,key = lambda x:x[2]) # make sure points are in decreasing order in time! 
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
                    if d*1./t>0.3 and verbose: print "# Too slow, v=%f" % (d*1./t)                
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


#def flightworks_old(Points,i,j,R): #wrong!
#    """ From a group of points, select minimal group that corresponds a flight
#    input: 
#        - Group of points (x,y,s) [list]
#        - Starting index of flight
#        - Ending index of flight
#        - Tolerance parameter
#    output:
#        - Index of minimal group of updates that correspond to a flight -> points[i:j+1]
#    """
#    # This is wrong
#    flag=False
#    while j>i+1 and not flag:
#        flag=True
#        for k in iter(range(i+1,j)):
#            dy=geom.distancepointtoline(Points[i][0:2],Points[j][0:2],Points[k][0:2])
#            rik=np.linalg.norm(Points[i][0:2]-Points[k][0:2])*np.linalg.norm(Points[i][0:2]-Points[k][0:2])
#            rij2=np.linalg.norm(Points[i][0:2]-Points[j][0:2])*np.linalg.norm(Points[i][0:2]-Points[j][0:2])
#            dx2=rik-dy*dy
#            if dy>R and dx2<rij2 :
#                flag=False
#                j-=1
#                break
#    return j

def flightworks(Points,R,R_stop): #works indeed!
    """ From a group of points, select minimal group that corresponds a flight (see paper)
    input: 
        - Group of points (x,y,s) [list]
        - Tolerance parameter
        - Stop parameter
    output:
        - Index of minimal group of updates that correspond to a flight -> points[i:j+1]
    """
    flag=False
    Points=np.array(Points).T[:-1].T
    j=len(Points)
    while j>1 and not flag:
        if np.linalg.norm(Points[j-1]-Points[0]) >= R_stop: #check sufficient distance
            flag=geom.pointsinbox(Points[:j],R) # check in box
            #print flag,j
        else:
            flag=False
        j-=1            
    return j




### Models ####
def rectangularmodel(Points,R,R_stop):
    """Split the trace in flights according to rectangular model (see aditional docs)
    input:
        - UTM points x-y-t-s as obtained from stopormove()
        - R parameter tolerance
        - R_stop parameter (from stopormove)
    output:
        - Track Class with Flights and Stops
    """
    Points=np.array(Points)
    Points = sorted(Point,key = lambda x:x[2]) # make sure points are in decreasing order in time! 
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
                aux_fl=[] #list of flights to add
                while len(dummy)>2: # means not only one or two update: detect flights
                    k=flightworks(dummy,R,R_stop) # get minimum flight from dummy                        
                    aux_fl.append(cls.Flight(np.array(dummy[0:k+1]))) # append flight
                    dummy=dummy[k:] # remove points from dummy list except last point
                if len(dummy)==2: # two updates surviving: last flight of the set 
                    q = stopormove(dummy,R_stop)
                    if q[1][-1] == 1: #(or not, check condition)
                        #qq = raw_input("aaa")
                        aux_fl.append(cls.Flight(np.array(dummy)))
                for fl in aux_fl: 
                    ######  Check function ###### (some flights might be smaller than R_stop by bad luck)
                    #if fl.Delta_r<R_stop:
                        #print fl
                        #q= raw_input("Proceed?")
                        #if q=='N':
                            #return fl
                    Track.add_flight(fl) #add flights
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





