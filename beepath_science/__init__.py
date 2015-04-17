import constants
import formats
import geodesics
import geometrics
import misc
#import mapping
import models
import filters
import stats
import classes
import tests

try:
    import numpy as np
    import scipy as sp
    import pyproj as pr
except ImportError:
    print """Install apropiate requirements:
        Requires:
        Pymongo
        Numpy/Scipy
        Pyproj
    Optionally:
        Geojson        
    """

__all__ = ["constants" ,"formats" ,"geodesics" ,"geometrics" ,"misc", "models", "filters", "stats", "classes", "tests"] #"mapping"




""" 
    Requires:
        Pymongo
        Numpy/Scipy
        Pyproj
    Optionally:
        Geojson        
"""

