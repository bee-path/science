# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 09:54:41 2013

@author: Oleguer Sagarra <osagarra@ub.edu>

This module handles the format transforms of the packages
"""

# External modules
import numpy as np

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
