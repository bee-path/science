#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       aux_module.py
#
#       Copyright 2012 Oleguer Sagarra <oleguer.sagarra@gmail.com>
#       and Mario Gutiérrez-Roig <mariogutierrezroig@gmai.com>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#
#       Contains:
#             General functions
#             Class definitions:
#                Class Trace: GPS Trace by each user, composed of Flights and Stops
#                Class Flight: Groups of updates consituting a flight
#                Class Stop: Inherited class from Flight, special type of flight where user does not move
#      requires several modules to work: numpy, geopy, simplekml, time, calendar...
#

### General variables definition ###

zone=31 # zone utm 31T
origin=(431769,4582211) # map origin #Coordenades geodèsiqueS: N41º23.312 E02º11.034
eps=0.0000001 # very small value (to avoid problems)

### Modules ###

import numpy as np
from scipy import stats
try:
    import pyproj as pr
except ImportError:
    print "Please install pyproj module"
utm = pr.Proj(ellips="WGS84", proj="utm", zone=zone, preserve_units=False) # Define ellipsoid and UMT zone

### Functions ###

##### Flight selections functions #####

def filter_v_updates(points,v_max=50/3.6,zone=zone,UTM=True):
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
        utm=pr.Proj(ellips="WGS84",proj="utm",zone=zone,preserve_units=False) # Define ellipsoid and UMT zone
        xy=np.array(utm(points.T[0],points.T[1])).T # transformacio geo (lon-lat)
        xy=xy-xy[0] # refer to arbitrary origin (to ease calculations)
    else:
        xy=points.T[:-1].T
    deltar=np.array(map(np.linalg.norm,np.diff(xy.T).T),dtype=float)
    deltat=np.array(np.diff(points.T[-1]),dtype=int) # time diff
    v=deltar/deltat
    group_points=[e for e in np.split(points,np.where(v>v_max)[0]+1) if len(e)>1] # only sets of points bigger than one (if not, means isolated bad point)
    return group_points


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
            if (d-R)<eps:
                mov[i]=0
        mov[-1]=0 # last point is a stop by default
        mov[0]=0 # first point is a stop by default
    else:
        mov=[0]
    return np.hstack((Points,np.reshape(mov,(len(Points),1))))

def separe_flights(oTrace):
    """from a trace gets a list of lists, where each list is a set of consecutive flights between stops"""
    alls=np.sort(list(oTrace.iter_all()))
    separe=[]
    i=0
    for e in alls:
        if isinstance(e,Stop):
            separe.append(i)
        i+=1
    fl_group=[]
    # group flights
    for g in np.split(alls,separe): 
        try:
            fl_group.append(g[1:])
        except IndexError: pass
    return [e for e in fl_group if any(e)]

def pausemodel(oTrace):
    """ merges successive flights from the rectangular model if their relative angle does not exceed theta_max
    input: theta_max, Trace containing stops and flights
    output: new Trace containing merged flights and same stops
    """
    new_Trace=Trace()
    # add stops
    for s in oTrace.iter_stops():
        new_Trace.add_stop(s) 
    fl_group=separe_flights(oTrace)
    for fl in fl_group: new_Trace.add_flight(np.sum(fl))
    return new_Trace


def angle_diff(points):
    """ given a set of vectors defined by delta_x,delta_y, get their angles """
    angles=np.arctan2(points.T[1],points.T[0])
    norms=np.hypot(points.T[1],points.T[0])
    ang_dif=np.empty(len(points)-1)
    for i in range(len(points)-1):
        dot=np.sum(points[i]*points[i+1])/(norms[i]*norms[i+1])
        ang_dif[i]=np.arccos(dot)
    return ang_dif
        
def angularmodel(oTrace,theta_max):
    """ merges successive flights from the rectangular model if their relative angle does not exceed theta_max
    input: theta_max, Trace containing stops and flights
    output: new Trace containing merged flights and same stops
    """
    new_Trace=Trace()
    # add stops
    for s in oTrace.iter_stops():
        new_Trace.add_stop(s) 
    fl_group=separe_flights(oTrace)
    # and now merge if needed
    for gr in fl_group:
        if len(gr)>1: # if more than 1 flight
            r=np.array([f.atts[2][-1][:-1]-f.atts[2][0][:-1] for f in gr])
            dif_angle=angle_diff(r)
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
def rectangularmodel(Points,R):
    """Split the trace in flights according to rectangular model (see aditional docs)
    input:
        - UTM points x-y-t-s as obtained from stopormove()
        - R parameter tolerance
    output:
        - Track Class with Flights and Stops
    """
    Track=Trace() # create Trace
    dummy=[] # create dummy list
    #if Points[0][3]==1: print "first a flight?" # first is a flight (this is impossible) --> checked already
    for i in iter(range(len(Points)-1)): # For all points except the last one (since we check forward)
        point=Points[i][0:3] # store point
        dummy.append(point) # add to dummy list
        if Points[i][3]!=Points[i+1][3]: #if change (last point has to be a stop by default (see else at the bottom)
            if Points[i][3]==0: # stop -> add stop from dummy array
                s=Stop(np.array(dummy))
                Track.add_stop(s)
                dummy=[point] # keep last point () (to add to next flight)
                #raw_input('added stop %s \n index: %i' % (s,i))
            else: # flight
                dummy.append(Points[i+1][0:3]) # append next stop [end point, flight, next-stop point]
                #if len(dummy)<3: print "Warning!!" #(this is impossible, already checked!)
                if len(dummy)==3: # means just one update with flight [end point, flight] or simply [flight]
                    f=Flight(np.array(dummy))
                    Track.add_flight(f)
                    #raw_input('added single flight %s \n index: %i' % (f,i))
                else: # more than 1 update
                    aux_fl=[] #list of flights to add                   
                    while len(dummy)>2: # means not only one or two update: detect flights
                        j=len(dummy)-1 # last point
                        k=flightworks(dummy,0,j,R) # get minimum flight from dummy                        
                        aux_fl.append(Flight(np.array(dummy[0:k+1]))) # append flight
                        dummy=dummy[k:] # remove points from dummy list except last point
                    if len(dummy)==2: # two updates surviving: last flight of the set
                        aux_fl.append(Flight(np.array(dummy)))
                    #if len(dummy)==1: pass  #one update survives by construction
                    for fl in aux_fl: Track.add_flight(fl) #add flights
                dummy=[] # empty list
        else: # nothing to do except if last point
            if i==len(Points)-2: # two last points are stops 
                dummy.append(Points[i+1][0:3]) # add last point to dummy
                s=Stop(np.array(dummy)) # add stop
                Track.add_stop(s)
    return Track

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
    flag=False
    while j>i+1 and not flag:
        flag=True
        for k in iter(range(i+1,j)):
            d=distancepointtoline(Points[i][0:2],Points[j][0:2],Points[k][0:2])
            if d>R:
                flag=False
                j-=1
                break
    return j

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

##### Format handling functions ####

def format_input(raw_data):
    """ Transform raw_input: update id lat lon str_day str_time accuracy (896 68 0.0000000 0.0000000 2012-06-08 18:42:38 0.0)
    into appropiate lat lon Unix_time
    input:
        - matrix with one value per column
    ouput:
        - matrix with id lat lon unix_time at each line)
    """
    np_data=np.array(raw_data)
    np_data.T[-1]=string_to_UNIX([r[-3]+" "+r[-2] for r in raw_data])
    return np.array([np.array([r[1],r[2],r[3],r[-1]]) for r in np_data],dtype=float)

def string_to_UNIX(times):
    """ Transform input times points into unix time
    input:
        - list with times as item strings: a time string of the form yyyy-mm-dd hh:mm:ss (for instnace: 2012-06-08 18:42:38)
    ouput:
        - list with unix times as float point numbers
    Note: Needs calendar and time modules to work
    """
    try:
        import time
        import calendar
    except ImportError:
        " Needs time and calendar modules to work "
        return
    return np.array([calendar.timegm(time.strptime(a,"%Y-%m-%d %H:%M:%S")) for a in times])

def unique(a):
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1) 
    return a[ui]
##### Geodesic functions ######

def UTM_2_latlon(points,origin,zone):
    """ Converts (x-y) input to lat-lon files on given UTM zone with given origin
    input:
        - origin must be a 2 index tupple (x,y)
        - np.array with at least 2 columns x-y
        - Scalar UTM zone
    output:
        - array with lat,lon
    """
    points=np.array(points)
    coords=utm(points.T[0]+origin[0],points.T[1]+origin[1],inverse=True)
    return np.array(coords)[::-1].T # transformacio geo (lat-lon)

def latlon_2_UTM(points,origin,zone,time=False):
    """ Converts lat-lon-time input to x-y-files on given UTM zone with given origin
    input:
        - origin must be a 2 index tupple (x,y)
        - np.array with at least 2 columns lat-lon
        - Scalar UTM zone
    output:
        - array with lat,lon
    """
    points=np.array(points)
    xy=np.array(utm(points.T[0],points.T[1])).T-origin # transformacio geo (lon-lat)
    if time:
        return np.array([xy.T[0],xy.T[1],points.T[-1]]).T # transformacio geo (lon-lat)
    else:
        return np.array(xy) # transformacio geo (lon-lat)

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

def gmaps(trace,path,details=True,stops=False,det_points=False):
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
    for num,fli in trace.flights.iteritems(): #for each flight
        if details:
            lin=hubmap.newlinestring(name="Flight %i"%(num), description="Real flight %i" % (num),coords=fli.atts[fli.att_dic['latlon']].T[::-1].T.tolist())
            lin.style.linestyle.color = 'ffff0000'  # Blue
            lin.style.linestyle.width= 3  # 10 pixels
            i=0 # Update number
            if det_points:
                for point in fli.atts[fli.att_dic['latlon']]:  # for each point
                    i+=1
                    hubmap.newpoint(name="Flight %i, update %i" % (num,i), coords=[point[::-1]])
        lin=hubmap.newlinestring(name="Flight %i"%(num), description="Filtered flight %i" % (num),coords=fli.atts[fli.att_dic['StartEnd_latlon']].T[::-1].T.tolist())
        lin.style.linestyle.color = 'ff00000f'  # Red
        lin.style.linestyle.width= 10  # 10 pixels
    if stops:
        for num,fli in trace.stops.iteritems(): #for each flight
            if details:
                lin=hubmap.newlinestring(name="Stop %i"%(num), description="Real Stop %i" % (num),coords=fli.atts[fli.att_dic['latlon']].T[::-1].T.tolist())
                lin.style.linestyle.color = 'fffff001'  # Blue
                lin.style.linestyle.width= 3  # 10 pixels
                i=0 # Update number
                if det_points:
                    for point in fli.atts[fli.att_dic['latlon']]:  # for each point
                        i+=1
                        hubmap.newpoint(name="Stop %i, update %i" % (num,i), coords=[point[::-1]])
            hubmap.newpoint(name="Mean Stop %i" % (num), coords=[fli.atts[fli.att_dic['mean_point_latlon']][::-1]])
    hubmap.save(path)
    return

### Classes ###

class Trace(object):
    """ This is a GPS trace object container of Flights and Stops
    type Trace.func_list for a list of available functions
    """
    att_dic={'N_flights':1,'N_stops':2}
    num_atts=len(att_dic)

    ## Basic funcs of the class ##

    def __init__(self,flights=None):
        self.atts={}
        self.atts[self.att_dic['N_flights']]=0
        self.atts[self.att_dic['N_stops']]=0
        self.flights={}
        self.stops={}
        if not flights:
            pass
            #print "Your Trace upon initialization is empty, be aware of it!"
        else:
            for fli in flights:
                if isinstance(fli,Flight):
                    if isinstance(fli,Stop):
                        self.add_stop(fli)
                    else:
                        self.add_flight(fli)
                else:
                    print "Need to pass an object flight or stop as argument to initialize"

    def __str__(self):
        return " This is a Trace object consistent of %i flights and %i stops " % (self.atts[self.att_dic['N_flights']],self.atts[self.att_dic['N_stops']])
    def __add__(self,other):
        """ merges traces into single trace """
        all_1f=list(self.iter_flights())
        all_2f=list(other.iter_flights())
        all_1s=list(self.iter_stops())
        all_2s=list(other.iter_stops())
        allf=np.hstack((all_1f,all_2f))
        alls=np.hstack((all_1s,all_2s))
        new_Trace=Trace()
        for f in allf: new_Trace.add_flight(f)
        for f in alls: new_Trace.add_stop(f)
        return new_Trace
    def add_flight(self,flight):
        """ Adds flight to trace, can pass flight object or simply points x-y-t"""
        if not isinstance(flight,Flight):
            flight=Flight(flight)
        self.atts[self.att_dic['N_flights']]+=1
        self.flights[self.atts[self.att_dic['N_flights']]]=flight

    def add_stop(self,stop):
        """ Adds stop to trace, can pass stop object or simply points x-y-t"""
        if not isinstance(stop,Stop):
            stop=Stop(stop)
        self.atts[self.att_dic['N_stops']]+=1
        self.stops[self.atts[self.att_dic['N_stops']]]=stop

    def add_atts(self,att,value):
        """ Add new atts to the Trace and to the att_dict """
        self.num_atts+=1
        self.att_dic[att]=self.num_atts
        self.atts[self.num_atts]=value

    def del_flight(self,flight_id):
        """ Deletes flight from trace according to flight_id"""
        self.flights.pop(flight_id)
        self.atts[self.att_dic['N_flights']]-=1

    def del_stop(self,stop_id):
        """ Deletes stop from trace according to stop_id """
        self.stops.pop(stop_id)
        self.atts[self.att_dic['N_stops']]-=1

    def del_atts(self,att):
        """ deletes atts from Trace """
        self.num_atts-=1
        self.atts.pop(self.att_dic[att])
        self.att_dic.pop(att)

    def length(self,att='UTM',stops=False):
        """ Counts lenght of all items of iterator
        if stops: iterates over stops """
        leng=0
        if stops:
            for ele in self.iter_stops(att=att): leng+=len(ele)
        else:
            for ele in self.iter_flights(att=att): leng+=len(ele)
        return leng

    ## Iterators ##
    def iter_all(self,att=None,count=False):
        """ Generates an iterator over all flights and stops sorted, same input as iter_flights"""
        fl=[f for f in self.iter_flights(att=att,count=count)]
        st=[f for f in self.iter_stops(att=att,count=count)]
        all_fs=np.hstack([fl,st])
        try:
            for e in np.sort(all_fs): yield e
        except ValueError:
            yield None
    def iter_flights(self,att=None,count=False,stops=False):
        """ Generates an iterator over flights
        OPt attributeS: att (gets the value of 'att' if exists in flights'
                        count: gets a list of tuples labelled (flight num, att)
        """
        if stops:
            itera=self.stops
        else:
            itera=self.flights
        if count:
            if att:
                for num,fli in itera.iteritems():
                    try:
                        yield (num,fli.atts[fli.att_dic[att]])
                    except KeyError:
                        print 'This att does not exist, None returned'
                        yield None
            else:
                for num,fli in itera.iteritems():
                    yield (num,fli)
        else:
            if att:
                for fli in itera.itervalues():
                    try:
                        yield fli.atts[fli.att_dic[att]]
                    except KeyError:
                        print 'This att does not exist, None returned'
                        yield None
            else:
                for fli in itera.itervalues():
                    yield fli

    def iter_stops(self,att=None,count=False,stops=False):
        """ Generates an iterator over stops
        OPt attributeS: att (gets the value of 'att' if exists in flights'
                        count: gets a list of tuples labelled (flight num, att)
        """
        return self.iter_flights(att,count,stops=True)

    def allpoints(self,stops=False):
        """ Generates an array of all points in the trace with UTM
        if stops: counts also points considered as stopped
        """
        if stops:
            l=(self.length(att='UTM',stops=stops)+self.length(att='UTM',stops=False),2)
        else:
            l=(self.length(att='UTM',stops=stops),2)
        all_points = np.empty(l)
        i=0
        for fli in self.iter_flights('UTM'):
            all_points[i:i+len(fli)]=fli.T[:-1].T
            i+=len(fli)
        if stops:
            for fli in self.iter_stops('UTM'):
                all_points[i:i+len(fli)]=fli.T[:-1].T
                i+=len(fli)
        return all_points

    #### General funcs #####

    func_list=['N_points','r_gyr','total_time','total_time_flights','total_time_stops','total_length','start_end','v','update_freq','r_sigma','tortuosity','type_mov']
    def N_points(self,extended=False,add=True):
        """ Returns total number of points in flights,stops in the trace. Also returns number of flights, stops in the trace
        If add: True, also adds to trace stats
        If extended:
            Output: Tuple of number of points per flight, Tuple of number of points per stop, total number of points in all flights, total number of points in all stops, total number of flights, total number of stops
        else:
            Output: total number of points (int) (flights and stops), total number of flights(stops) (int)
        """
        a=[fli for fli in self.iter_flights(att='N_points')]
        b=[fli for fli in self.iter_stops(att='N_points')]
        if add:
            self.add_atts('N_points_flights',np.sum(a))
            self.add_atts('N_points_stops',np.sum(b))
        if extended:
            return tuple(a),tuple(b),np.sum(a),np.sum(b),self.atts[self.att_dic['N_flights']],self.atts[self.att_dic['N_stops']]
        else:
            return np.sum(a),np.sum(b),self.atts[self.att_dic['N_flights']],self.atts[self.att_dic['N_stops']]

    def rgyr(self,add=True,stops=True):
        """ Computes CM and rgyr of the Trace using all points in UTM
        If add: True, adds also as flight att
        if stops: True counts as well stopped points
        """
        a=(gyrad(self.allpoints(stops=stops))[-2],gyrad(self.allpoints(stops=stops))[-1])
        b=tuple(gyrad(self.allpoints(stops=stops))[0])
        if add:
            self.add_atts('r_gyr',a)
            self.add_atts('CM',b)
        return a,b

    def start_end(self,unit='UTM'):
        """ Obtain the list of starting-ending positions of flights in the trace
        in units: 'latlon' or 'UTM'
        output: list x-y-t or list lat-lon
        """
        mode={'UTM':'StartEnd_UTM','latlon':'StartEnd_latlon'}
        if unit=='UTM':
            st_end=np.empty((self.N_points()[-2]*2,3))
        else:
            st_end=np.empty((self.N_points()[-2]*2,2))
        att=mode[unit]
        for num,fli in self.iter_flights(att,count=True):
            st_end[2*num-2]=fli[0]
            st_end[2*num-1]=fli[1]
        return st_end

    def start(self,unit='UTM',stops=True):
        """ Obtain the list of starting positions of flights and stops (if stops=True) in ordered manner in the trace
        in units: 'latlon' or 'UTM'
        output: list x-y-t or list lat-lon
        """
        mode={'UTM':'StartEnd_UTM','latlon':'StartEnd_latlon'}
        if stops:
            len_arr=self.N_points()[-2]+self.N_points()[-1]
        else:
            len_arr=self.N_points()[-2]
        if unit=='UTM':
            st=np.empty((len_arr,3))
        else:
            st=np.empty((len_arr,2))
        att=mode[unit]
        for num,fli in self.iter_flights(att,count=True):
            st[num-1]=fli[0]
        num_0=self.N_points()[-2]-1
        if stops:
            for num,fli in self.iter_stops(att,count=True):
                st[num+num_0]=fli[0]
        return st[st[:,-1].argsort()]

    def total_time(self,add=True):
        """ Compute total time of trace
        If add: True, then add also to flight atts
        """
        a=self.total_time_flights()+self.total_time_stops()
        if add: self.add_atts('Delta_t_total',a)
        return a

    def total_time_flights(self,add=True):
        """ Compute total time of flights in the trace
        If add: True, then add also to flight atts
        """
        a=np.sum(self.iter_flights('Delta_t'))
        if add: self.add_atts('Delta_t_flights',a)
        return a

    def total_time_stops(self,add=True):
        """ Compute total time of stops in the trace
        If add: True, then add also to flight atts
        """
        a=np.sum(self.iter_stops('Delta_t'))
        if add: self.add_atts('Delta_t_stops',a)
        return a

    def total_length(self,add=True):
        """ Compute total length of flights in the trace
        If add: True, then add also to flight atts
        """
        a=np.sum(self.iter_flights('Delta_r'))
        if add: self.add_atts('Delta_r',a)
        return a

    def v(self,add=True,kmh=False):
        """ Compute mean and std of mean velocities of flights in the trace
        If add: True, then add also to flight atts
        if kmh: True, result in km/h
        """
        if kmh:
            norm=3.6
        else:
            norm=1.
        a=np.mean([ele.atts[ele.att_dic['Delta_r']] / float(ele.atts[ele.att_dic['Delta_t']]) for ele in self.iter_flights()]),\
                np.std([ele.atts[ele.att_dic['Delta_r']] / float(ele.atts[ele.att_dic['Delta_t']]) for ele in self.iter_flights()])
        if not np.isfinite(a[0]): a=np.array([0,0])
        if add:
            self.add_atts('v',(a[0]*norm,a[1]*norm))
        return (a[0]*norm,a[1]*norm)

    def update_freq(self,add=True):
        """ Computes the mean time and std frequency of sucessive updates
            for flights and stops
        If add: True, then add also to flight atts
        output:
            - flights mean
            - flights std
            - stops mean
            - stops std
        """
        a=np.mean([fli.update_freq(add=False) for fli in self.iter_flights()],axis=0)
        b=np.mean([fli.update_freq(add=False) for fli in self.iter_stops()],axis=0)
        try:
            if not any(np.isfinite(a)): a=np.array([0,0])
        except TypeError:
            if not np.isfinite(a): a=np.array([0,0])
        try:
            if not any(np.isfinite(b)): b=np.array([0,0])
        except TypeError:
            if not np.isfinite(b): b=np.array([0,0])
        if add:
            self.add_atts('t_freq',(a,b))
        return a,b

    def r_sigma(self,add=True,div=1.):#not implemented
        """ Calculates exponent gamma of sigma(r)~t^{gamma} of the given points in UTM format (x,y,t)
        input:
            - array with (x,y,t), which alternatively are the start and the end of flights & stops sorted in time
        output:
            - gamma value,goodness of fit value,real points,fitted func
        """
        if not self.N_points()[-2]: # not flights
            gamma,p_value,rt,x,y=0,0,np.array([0,0]),0,0
        else:
            try:
                from collections import defaultdict
            except ImportError:
                print "Needs collections module "
                return
            rt=defaultdict(float)
            rcount=defaultdict(float)
            points=self.start(stops=False)
            N=len(points)
            for d in iter(range(1,N)):
                for i in iter(range(0,N-d)):
                    r0=points[i][0:2]
                    r1=points[i+d][0:2]
                    dr2=np.linalg.norm(r1-r0)*np.linalg.norm(r1-r0)
                    dt=points[i+d][2]-points[i][2]
                    rt[dt]+=dr2
                    rcount[dt]+=1
            for key,val in rt.items():
                if val==0:
                    del rt[key]
                else:
                    rt[key]=rt[key]/float(rcount[key])
            rt=np.array([ (k,rt[k]) for k in sorted(rt.keys())])
            rt_log=np.log10(rt[:int(len(rt)/div)])
            if N>1:
                gamma, intercept, r_value, p_value, std_err = stats.linregress(rt_log.T[0],rt_log.T[1])
            else: # It did not work
                print "Not enough data"
                return 0,0,0,0,0
            #for a,b in rt:
            #    print a,b
            #print gamma, intercept, r_value, p_value, std_err
            if add:
                self.add_atts('gamma',(gamma,p_value))
            x=np.linspace(min(rt.T[0]),max(rt.T[0]))
            y=x**gamma*10.**intercept
        return gamma,p_value,rt,np.array([x,y])

    def tortuosity(self,r=10.,unit='UTM'): #not implemented
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
        return 0

    def type_move(self,params={'v':(4,5),'r_gyr':(75)},mov_types={(1,1):'shark',(1,-1):'fox',(0,1):'dollar',(0,-1):'bee',(-1,1):'pollen',(-1,-1):'cell'},rand=False):
        """ Gets a type of animanl movement
        input:
            - params: set threshold for calssif (velocity and r_gyr) in (km/h and m)
            - mov_types: types of movement
            - rand: If true, make up the type (in case something fails)
        output:
            - type_mov: a string selected according to pre-set parametters dic
            - status: If rand True, else false
        """
        if rand:
            return mov_types.values().pop(np.random.randint(0,high=len(mov_types))),False
        else:
            v=self.v(kmh=True)[0]
            rgyr=self.rgyr()[0][0]
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

    def stop_flight_sequence(self):
        """Returns the sequence of flights and stops ordered according to the time"""
        fl = np.array([f for f in self.iter_flights()])
        st = np.array([s for s in self.iter_stops()])
        all = np.hstack([fl, st])
        all_sorted = np.sort(all)

        return all_sorted

    #### Statistical Functions  -> han de cridar funcions totals ####

    def stats_pdf_int(self,att): #not finished
        """ Returns histogram of given att over flights, converting elements to integers
        input:
            - att of flights
        output:
            - 2-D array (x_values,counts)
        """
        try:
            data=np.array(list(self.iter_flights(att)),dtype=int)
        except ValueError:
            print "This is not a scalar quantity for each flight, cannot make histrogram"
            return None
        # not finished
        return

    def stats_pdf_float(self,att):
        """ Return histogram of given att over flights, with float element using binning
        input:
            - att of flights
            - minbin: minimum number of data per bin
        output:
            - 2-D array (x_values (center of bin), counts)
        """
        #not implemented
        pass

    def stats_CDF_int(self,att,compl=False):
        """ Returns CDF of given att over flights, probability of finding an att < than a given element
        if compl is True, then returns CCDF """
        #not implemented
        pass

    def stats_stops_duration_PDF(self, hbin=10, hmin=0, hmax=100):
        """ Returns an histogram of stops duration
        input:
            - Number of bins
            - Minimum of the range
            - Maximum of the range
        output:
            - Array with counts
            - Array with bins
        """
        x = np.array([s for s in self.iter_stops(att='Delta_t')])
        hist, bins = np.histogram(x, hbin, [hmin,hmax])
        pdf = hist/float(max(np.cumsum(hist)))
        return pdf, bins

    def stats_flights_length_PDF(self, hbin=10, hmin=0, hmax=100):
        """ Returns histogram of flight length
         input:
            - Number of bins
            - Minimum of the range
            - Maximum of the range
        output:
            - Array with counts
            - Array with bins
        """
        x = np.array([f for f in self.iter_flights(att='Delta_r')])
        hist, bins = np.histogram(x, hbin, [hmin,hmax])
        pdf = hist/float(max(np.cumsum(hist)))
        return pdf, bins

    def stats_flights_time_PDF(self, hbin=10, hmin=0, hmax=100):
        """ Returns an histogram of flight times
        input:
            - Number of bins
            - Minimum of the range
            - Maximum of the range
        output:
            - Array with counts
            - Array with bins
        """
        x = np.array([f for f in self.iter_flights(att='Delta_t')])
        hist, bins = np.histogram(x, hbin, [hmin,hmax])
        pdf = hist/float(max(np.cumsum(hist)))
        return pdf, bins

    def stats_mean_std_stops_duration(self):
        """ Returns the mean and the standard deviation of stops duration
        input:
        output:
            - Mean
            - Standard deviation
        """
        x = np.array([s for s in self.iter_stops(att='Delta_t')])
        mean = np.mean(x)
        std = np.std(x)
        return mean, std

    def stats_mean_std_flights_duration(self):
        """ Returns the mean and the standard deviation of flights duration
        input:
        output:
            - Mean
            - Standard deviation
        """
        x = np.array([s for s in self.iter_flights(att='Delta_t')])
        mean = np.mean(x)
        std = np.std(x)
        return mean, std

    def stats_mean_std_flights_length(self):
        """ Returns the mean and the standard deviation of flights length
        input:
        output:
            - Mean
            - Standard deviation
        """
        x = np.array([s for s in self.iter_flights(att='Delta_r')])
        mean = np.mean(x)
        std = np.std(x)
        return mean, std

class Flight(object):
    """ This is a flight, needs numpy array as input with 3 cols: x,y,time """

    # Attributes

    att_dic={'latlon':1,'UTM':2,'N_points':3,'Delta_r':4,'Delta_t':5,'StartEnd_latlon':6,'StartEnd_UTM':7}
    num_atts=len(att_dic)

    # General methods

    def __init__(self,points=None):
        if len(points)==0:
            print "Need an array of x-y-time points to work"
        else:
            self.atts={}
            self.atts[1]=UTM_2_latlon(points,origin,zone)
            self.atts[2]=points
            self.atts[3]=len(self.atts[2])
            self.atts[4]=np.linalg.norm(self.atts[2][-1]-self.atts[2][0])
            self.atts[5]=self.atts[2][-1][-1]-self.atts[2][0][-1]
            self.atts[7]=np.array([self.atts[2][0],self.atts[2][-1]]) # in x y t
            self.atts[6]=UTM_2_latlon(self.atts[7],origin,zone) #in lat lon

    def __str__(self):
        return "## Flight ##: \n N_points=%i \t Length=%f \t Time=%f" % (self.atts[3],self.atts[4],self.atts[5])
    def __eq__(self,other):
        """ determines if two flights/stops are equal based on time of start and type """
        if  not isinstance(self,Stop) == isinstance(other,Stop): return False
        else:
            return self.atts[self.att_dic['UTM']][0][-1]==other.atts[self.att_dic['UTM']][0][-1]
    def __cmp__(self, other):
        """ Compares flights in temporal order always putting stops in front of flights in case of same time record"""
        if self.atts[self.att_dic['UTM']][0][-1] == other.atts[self.att_dic['UTM']][0][-1]:
            if  not isinstance(self,Stop) == isinstance(other,Stop):
                if isinstance(self,Stop): return -1
                else: return 1
            else: print "Warning! Two equal flights or stops at same time?"
        else: return cmp(self.atts[self.att_dic['UTM']][0][-1],other.atts[self.att_dic['UTM']][0][-1])
    def __add__(self,other):
        """ sums flights: Assigns all positions of the flight to a single flight only if they are consecutive"""
        points=np.vstack((self.atts[2],other.atts[2]))
        return Flight(unique(points))
    def add_atts(self,att,value):
        """ Add new atts to the flight and to the att_dict """
        self.num_atts+=1
        self.att_dic[att]=self.num_atts
        self.atts[self.num_atts]=value
        return

    def del_atts(self,att):
        """ Deletes atts  of the flight """
        self.atts.pop(self.att_dic[att])
        self.att_dic.pop(att)
        self.num_atts-=1
        return

    # Functions

    def rgyr(self,add=True):
        """ Computes CM and rgyr of the flight (just for checking
        If add=True, then adds also to flight attributes
        """
        a=(gyrad(self.atts[2])[-2],gyrad(self.atts[2])[-1])
        b=tuple(gyrad(self.atts[2])[0])
        if add:
            self.add_atts('r_gyr',a)
            self.add_atts('CM',b)
        return a,b

    def v(self,kmh=False,add=True):
        """ Compute mean velocity of flight
        If kmh is True, then output in  km/h, otherwise in metters/s
        If add=True, then adds also to flight attributes
        """
        if kmh:
            norm=3.6
        else:
            norm=1.
        a=self.atts[self.att_dic['Delta_r']] * norm / float(self.atts[self.att_dic['Delta_t']])
        if add:
            self.add_atts('v',a)
        return a

    def update_freq(self,add=True):
        """ Computes update freq in time of flight (mean)
        output:
            - 2-item tuple (freq_t,freq_t_sigma) with mean and std
        If add=True, then adds also to flight attributes
        """
        if len(self.atts[self.att_dic['UTM']].T[-1])!=1:
            a=(np.mean(np.diff(self.atts[self.att_dic['UTM']].T[-1])),np.std(np.diff(self.atts[self.att_dic['UTM']].T[-1])))
        else:
            a=(0,0)
        if add:
            self.add_atts('t_freq',a)
        return a

class Stop(Flight):
    """ Special type of flight: Stopped flight """
    att_dic={'latlon':1,'UTM':2,'N_points':3,'Delta_t':4,'mean_point_latlon':5,'mean_point_UTM':6,'StartEnd_latlon':7,'StartEnd_UTM':8}

    def __init__(self,points=None,Stop=False):
        if len(points)==0:
            print "Need an array of lat-lon-time points to work"
        else:
            self.atts={}
            self.atts[2]=points
            self.atts[1]=UTM_2_latlon(self.atts[2],origin,zone)
            self.atts[3]=len(points)
            self.atts[4]=points[-1][-1]-points[0][-1]
            self.atts[6]=np.mean(points,axis=0) # Mean point
            self.atts[5]=UTM_2_latlon(self.atts[6],origin,zone)
            self.atts[8]=np.array([np.array([self.atts[6][0],self.atts[6][1],self.atts[2][0][-1]]),np.array([self.atts[6][0],self.atts[6][1],self.atts[2][-1][-1]])])
            self.atts[7]=UTM_2_latlon(self.atts[8],origin,zone) #in lat lon

    def __str__(self):
        return "## Stoping point  ##: \n N_points=%i \t Position=(%f,%f) \n Time=%f " % (self.atts[3],self.atts[6][0],self.atts[6][1],self.atts[6][-1])

    def rgyr(self):
        """ Not defined for stops! """
        print "Not defined for stopped points"
        return None

    def v(self):
        """ Not defined for stops! """
        print "Not defined for stopped points"
        return None

def main(GPSpoints,R1=1.,R2=10.,rand=False):
    """ gets the group of GPSpoints (lat,lon,t) of a given user and returns a Trace object with flights and stops. It also returns a tuple of values
    input:
        - GPS raw points (lat,lon,unixtime)
        - R1: tolerance index stop-move
        - R2: tolerance index flights
        - rand: Randomization of movent stats...
    output:
        - Trace
        - Tuple: Total_time,Total_time_moving, total_time_stopped,total_lenght,mean velocity (when moving), radius of gyration, gamma exponent (of MSD),list of starting-ending points of flights (to represent),movement type
    """
    Trac=rectangularmodel(stopormove(latlon_2_UTM(GPSpoints,origin,zone,time=True),R1),R2)
    return Trac

if __name__ == '__main__':
    input_data=np.genfromtxt('../samples/sample.latlon') # here input data
    main(input_data)

#### Log #####
# TO DO:
# Escriure r_sigma
# Falta stats, pdf, CCDF i fitter.
# Falta server issues : written format_handling functions

# To do:
# Stopormove -> average positions ?
# angularmodel
# pause model

