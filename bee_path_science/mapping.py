# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:05:56 2013

@author: Oleguer Sagarra <osagarra@ub.edu>

This modules helps in quick mapping for results, requires "simplekml module" to work
"""

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
