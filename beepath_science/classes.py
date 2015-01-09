#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       classes.py
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
#             Class definitions:
#                Class Trace: GPS Trace by each user, composed of Flights and Stops
#                Class Flight: Groups of updates consituting a flight
#                Class Stop: Inherited class from Flight, special type of flight where user does not move
#      requires several modules to work: numpy, geopy, simplekml, time, calendar...
#



### Modules ###

# Internal #
import constants as cnt
import geodesics as geod
import stats as st
import misc
import geometrics as geom
# External #
import numpy as np
#from scipy import stats as sts
#import pyproj as pr


### General variables definition ###

zone=cnt.zone # zone utm 31T
origin=cnt.origin # map origin #Coordenades geodèsiqueS: N41º23.312 E02º11.034
eps=cnt.eps # very small value (to avoid problems)
utm = cnt.utm() # Define ellipsoid and UMT zone



### Classes ###

class Trace(object):
    """ This is a GPS trace object container of Flights and Stops
    type Trace.func_list for a list of available functions
    """
        
    
    ## Basic funcs of the class ##

    def __init__(self,flights=None):
        self.__setattr__('N_flights',0) #atts[self.att_dic['N_flights']]=0
        self.__setattr__('N_stops',0) #atts[self.att_dic['N_flights']]=0
        self.flights=[]
        self.stops=[]
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
        return " This is a Trace object consistent of %i flights and %i stops " % (self.N_flights,self.N_stops)
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
        self.N_flights+=1
        self.flights.append(flight)

    def add_stop(self,stop):
        """ Adds stop to trace, can pass stop object or simply points x-y-t"""
        if not isinstance(stop,Stop):
            stop=Stop(stop)
        self.N_stops+=1
        self.stops.append(stop)

    def add_atts(self,att,value,subst=True):
        """ Add new atts to the Trace and to the att_dict """
        try:        
            z=self.__getattribute__(att)
            if subst:            
                #print "Warning, substituting attribute!"
                self.__setattr__(att,value)
                return True
            else:
                #print "Att Already exists, nothing done"
                return False
        except AttributeError:
            self.__setattr__(att,value)
            return True

    def del_flight(self,flight_id):
        """ Deletes flight from trace according to flight_id"""
        self.flights.pop(flight_id)
        self.N_flights-=1

    def del_stop(self,stop_id):
        """ Deletes stop from trace according to stop_id """
        self.stops.pop(stop_id)
        self.N_stops-=1

    def del_atts(self,att):
        """ deletes atts from Trace """
        self.__delattr__(att)

    def iter_length(self,att='UTM',stops=False):
        """ Counts lenght of all items of iterator
        if stops: iterates over stops """
        if stops:
            obj=self.__getattribute__('iter_stops')
        else:
            obj=self.__getattribute__('iter_flights')            
        return sum(len(ele) for ele in obj(att))

    ## Iterators ##
    def iter_all(self,att=None,count=False):
        """ Generates an iterator over all flights and stops sorted, same input as iter_flights"""
        fl=[f for f in self.iter_flights(att=att,count=count)]
        sto=[f for f in self.iter_stops(att=att,count=count)]
        all_fs=np.hstack([fl,sto])
        try:
            for e in sorted(all_fs): yield e
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
                for num,fli in enumerate(itera):
                    try:
                        yield (num,fli.__getattribute__(att))
                    except AttributeError:
                        print 'This att does not exist, None returned'
                        yield None
            else:
                for num,fli in enumerate(itera):
                    yield (num,fli)
        else:
            if att:
                for fli in itera:
                    try:
                        yield fli.__getattribute__(att)
                    except AttributeError:
                        print 'This att does not exist, None returned'
                        yield None
            else:
                for fli in itera:
                    yield fli

    def iter_stops(self,att=None,count=False,stops=False):
        """ Generates an iterator over stops
        OPt attributeS: att (gets the value of 'att' if exists in flights'
                        count: gets a list of tuples labelled (flight num, att)
        """
        return self.iter_flights(att,count,stops=True)

    def allpoints(self,stops=False,time=False):
        """ Generates an array of all points in the trace with UTM
        if stops: counts also points considered as stopped
        """
        if stops:
            l=(self.iter_length(att='UTM',stops=stops)+self.iter_length(att='UTM',stops=False),3)
        else:
            l=(self.iter_length(att='UTM',stops=stops),3)
        all_points = np.empty(l)
        i=0
        for fli in self.iter_flights('UTM'):
            all_points[i:i+len(fli)]=fli
            i+=len(fli)
        if stops:
            for fli in self.iter_stops('UTM'):
                all_points[i:i+len(fli)]=fli
                i+=len(fli)
        order=np.argsort(all_points.T[-1])
        if time:
            return all_points[order]
        else:
            return all_points[order].T[:-1].T

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
            return tuple(a),tuple(b),np.sum(a),np.sum(b),self.N_flights,self.N_stops
        else:
            return np.sum(a),np.sum(b),self.N_flights,self.N_stops

    def rgyr(self,add=True,stops=True):
        """ Computes CM and rgyr of the Trace using all points in UTM
        If add: True, adds also as flight att
        if stops: True counts as well stopped points
        """
        a=st.gyrad(self.allpoints(stops=stops))[-1]
        b=tuple(st.gyrad(self.allpoints(stops=stops))[0])
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

    def v(self,kmh=False, add=True):
        """ Compute mean and std of mean velocities of flights in the trace
        If add: True, then add also to flight atts
        if kmh: True, result in km/h
        """
        if kmh:
            norm=3.6
        else:
            norm=1.
        a=np.mean([ele.Delta_r / float(ele.Delta_t) for ele in self.iter_flights()]),\
                np.std([ele.Delta_r / float(ele.Delta_t) for ele in self.iter_flights()])
        if not np.isfinite(a[0]): a=np.array([0,0])
        if add: self.add_atts('vel',a)
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

    def stop_flight_sequence(self):
        """Returns the sequence of flights and stops ordered according to the time"""
        fl = np.array([f for f in self.iter_flights()])
        st = np.array([s for s in self.iter_stops()])
        all = np.hstack([fl, st])
        all_sorted = np.sort(all)
        return all_sorted
    
    ############# Statistical functions ################
    def stats_PDF(self, att, opt='stops', hbin=10, hmin=0, hmax=100):
        """ Returns an histogram of att
        input:
            - Number of bins
            - Minimum of the range
            - Maximum of the range
            - option: flights (flights) or stops (stops) or all(all)
        output:
            - Array with counts
            - Array with bins 
        Warning: Does not check if att is numerical 
        """    
        obj=self.__getattribute__('iter_'+opt)
        if np.isfinite(obj(att).next()):
            return st.data_PDF(np.array([s for s in obj(att)]), hbin, hmin, hmax)    
        else:
            print "Not scalar quantity, nothing done"
            return 0

    def stats_CCDF(self, att, opt='stops', hbin=10, hmin=0, hmax=100):
        """ Returns CDF of att
        input:
            - Number of bins
            - Minimum of the range
            - Maximum of the range
            - option: flights (flights) or stops (stops) or all(all)
        output:
            - Array with counts
            - Array with bins
            Warning: Does not check if att is numerical 
        """
        obj=self.__getattribute__('iter_'+opt)
        if np.isfinite(obj(att).next()):
            return st.data_CCDF(np.array([s for s in obj(att)]), hbin, hmin, hmax)    
        else:
            print "Not scalar quantity, nothing done"
            return 0

    def r_sigma(self):
        """ Calculates exponent gamma of sigma(r)~t^{gamma} of the given points in UTM format (x,y,t)
        input:
        output:
            - gamma value,goodness of fit value,real points,fitted func
        """
        all_pts=self.allpoints(stops=True,time=True) # with stops
        return st.r_sigma(all_pts,1)

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
        return st.tortuosity()
        
    def turning_angle(self,consecutive=False):
        """ Returns collection of turning angles in the trace (in radians)
            If consecutive=True, only counts angles between consecutive flights
        """
        if self.N_flights < 2:
            return [False]
        else:
            if not consecutive:
                dp = np.array([g for g in self.iter_flights('StartEnd')]) #(pairs of points)
                angs = np.ediff1d(tuple(geom.angle_vect(g[0][:-1],g[1][:-1]) for g in dp))
            else:
                angs = []
                angs_d = []
                for f in self.iter_all():
                    if not isinstance(f,Stop):
                        angs_d.append(geom.angle_vect(f.StartEnd[0][:-1],f.StartEnd[-1][:-1]))
                    else:
                        if len(angs_d)>1:
	                        angs.extend(np.ediff1d(angs_d))
                        angs_d = []
                if len(angs_d)>1:
                    angs.extend(np.ediff1d(angs_d))
            for i,a in enumerate(angs):
                if a>np.pi:
                    angs[i] = - (2*np.pi-a)
                if a<-np.pi:
                    angs[i] = 2*np.pi+a
            return np.array(angs)

    def type_move(self,params={'v':(4,5),'r_gyr':(75)},mov_types={(1,1):'shark',(1,-1):'fox',(0,1):'dollar',(0,-1):'bee',(-1,1):'pollen',(-1,-1):'cell'},rand=False):
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
        raise DeprecationWarning(" Just in use for the 'festa' ")
        #return misc.type_move(Trace,params,rand):
        return None




class Flight(object):
    """ This is a flight, needs numpy array as input with 3 cols: x,y,time """

    # Attributes

    #att_dic={'latlon':1,'UTM':2,'N_points':3,'Delta_r':4,'Delta_t':5,'StartEnd_latlon':6,'StartEnd_UTM':7} # deprecated

    # General methods

    def __init__(self,points=None):
        if len(points)==0:
            print "Need an array of x-y-time points to work"
        else:
            points = sorted(points,key=lambda x: x[2]) # sort by increasing time
            self.__setattr__('UTM',points)
            ### Should implement the rest as methods
            self.__setattr__('N_points',len(self.UTM)) # as method
            self.__setattr__('Delta_r',np.linalg.norm(self.UTM[-1][:-1]-self.UTM[0][:-1])) # exclude time! (again method)
            self.__setattr__('Delta_t',self.UTM[-1][-1]-self.UTM[0][-1]) # method
            self.__setattr__('StartEnd',np.array([self.UTM[0],self.UTM[-1]])) #  method
    def __str__(self):
        return "## Flight ##: \n N_points=%i \t Length=%f \t Time=%f" % (self.N_points,self.Delta_r,self.Delta_t)
    def __eq__(self,other):
        """ determines if two flights/stops are equal based on time of start and type """
        if  not isinstance(self,Stop) == isinstance(other,Stop): return False
        else:
            return self.UTM[0][-1]==other.UTM[0][-1]
    def __cmp__(self, other):
        """ Compares flights in temporal order always putting stops in front of flights in case of same time record"""
        if self.UTM[0][-1] == other.UTM[0][-1]:
            if  not isinstance(self,Stop) == isinstance(other,Stop):
                if isinstance(self,Stop): return -1
                else: return 1
            else: 
                print "Warning! Two equal flights or stops at same time?"
                return 0
        else: return cmp(self.UTM[0][-1],other.UTM[0][-1])
    def __add__(self,other):
        """ sums flights: Assigns all positions of the flight to a single flight only if they are consecutive"""
        points=np.vstack((self.UTM,other.UTM))
        return Flight(misc.unique(points))

    
    def add_atts(self,att,value,subst=True):
        """ Add new atts to the flight """
        try:        
            z=self.__getattribute__(att)
            if subst:            
                #print "Warning, substituting attribute!"
                self.__setattr__(att,value)                
                return True
            else:
                #print "Att Already exists, nothing done"
                return False
        except AttributeError:
            return True

    def del_atts(self,att):
        """ deletes atts from flight """
        self.__delattr__(att)


    def lat_lon(self,att='UTM'):
        """ Coverts the given attribute to latlon"""
        if att!= 'UTM' and att!='StartEnd':
            print "Nothing done, this att is not a UTM position"
            return False
        else:
            return geod.UTM_2_latlon(self.__getattribute__(att),origin,zone)

    def rgyr(self,add=True):
        """ Computes CM and rgyr of the flight (just for checking
        If add=True, then adds also to flight attributes
        """
        a=st.gyrad(self.UTM)[-1]
        b=st.gyrad(self.UTM)[0]
        b=tuple(st.gyrad(self.UTM)[0])
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
        a=self.Delta_r * norm / float(self.Delta_t)
        if add:
            self.add_atts('vel',a)
        return a

    def update_freq(self,add=True):
        """ Computes update freq in time of flight (mean)
        output:
            - 2-item tuple (freq_t,freq_t_sigma) with mean and std
        If add=True, then adds also to flight attributes
        """
        if len(self.UTM.T[-1])!=1:
            a=(np.mean(np.diff(self.UTM.T[-1])),np.std(np.diff(self.UTM.T[-1])))
        else:
            a=(0,0)
        if add:
            self.add_atts('t_freq',a)
        return a

class Stop(Flight):
    """ Special type of flight: Stopped flight """
    
    def __init__(self,points=None,Stop=False):
        if len(points)==0:
            print "Need an array of lat-lon-time points to work"
        else:
            points = sorted(points,key=lambda x: x[2]) # sort by increasing time
            self.__setattr__('UTM',points)
            ### Should implement the rest as methods
            self.__setattr__('N_points',len(self.UTM)) # as method
            self.__setattr__('Delta_t',points[-1][-1]-points[0][-1]) # method
            self.__setattr__('mean_point',np.mean(points,axis=0) )
            self.__setattr__('StartEnd',
                            np.array([np.array([self.mean_point[0],self.mean_point[1],self.UTM[0][-1]]),
                                      np.array([self.mean_point[0],self.mean_point[1],self.UTM[-1][-1]])]))
    def __str__(self):
        return "## Stoping point  ##: \n N_points=%i \t Position=(%f,%f) \n Time=%f " % (self.N_points,self.mean_point[0],self.mean_point[1],self.mean_point[-1])

    def lat_lon(self,att):
        """ Coverts the given attribute to latlon"""
        if att!= 'UTM' and att!='StartEnd' and att!='mean_point':
            print "Nothing done, this att is not a UTM position"
            return False
        else:
            return geod.UTM_2_latlon(self.__getattribute__(att),origin,zone)
            
    def rgyr(self):
        """ Not defined for stops! """
        print "Not defined for stopped points"
        return None

    def v(self):
        """ Not defined for stops! """
        print "Not defined for stopped points"
        return None

#### Log #####
# TO DO:
# Escriure r_sigma, tortuosity, angle_diff


