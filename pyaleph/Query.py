import numpy as np
from . import DataBase
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from astropy.time import Time
import sqlite3, time, warnings
import numpy as np
import multiprocessing as mp
from astropy.coordinates.angle_utilities import angular_separation
from .aleph_utils import *
from astropy.table import QTable, Table, hstack

def getint(s):
    try: return int(s)
    except: return np.nan

class Query:
    def __init__(self, open_db=True, dbfile='aleph_states.db'):
        """
        Initializing asteroids' data
        """
        alephs = DataBase.AlerceEphs('Lowell')
        self.asts = alephs.asts
        nums = []; nms = []; hs = []; gs = []
        for ast in self.asts:
            nums.append(getint(ast['num'])); nms.append(ast['name'].strip())
            hs.append(ast['H']); gs.append(ast['G'])
        self.numbers = np.array(nums)
        self.names = np.array(nms)
        self.H = np.array(hs)
        self.G = np.array(gs)
        del alephs, nums, nms, hs, gs
        
    
    def query(self, field_center, radius, **kwargs):
        """
        Query for asteroides visible in a field.
        
        Parameters:
                    field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
                    radius: astropy.units.Quantity with radius of field.
                    
        Optional Parameters:
                    epoch: astropy.time.Time with epoch for query. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    njobs: Number of parallel processes to calculate asteroids' apparent positions. Default: 1.
                    method: '2B' or 'NB'. Default 'NB'.
                           '2B': Calculates all bodies positions using 2-Bodies equations and retrieves bodies in
                                 requested field.
                           'NB': Calculates all bodies positions using N-Bodies integration (Sun, planets, Pluto and 
                                 massless asteroids) and retrieves bodies in requested field.
                    ast_idxs: Asteroids' indexes (starting from 0) from Query.asts. Ephemerides will be
                              calculated only for those asteroids. Only if 'method' is '2B' or 'NB'
                    
        Returns:
                  astropy.table with results
        """
        # Selecting method
        if not 'method' in kwargs: return self.query_nb(field_center, radius, **kwargs)
        else:
            if kwargs['method']=='2B': return self.query_2b(field_center, radius, **kwargs)
            elif kwargs['method']=='NB': return self.query_nb(field_center, radius, **kwargs)
            else: raise ValueError("'method' must be one of 'db', '2B' or 'NB'")
    
    def query_2b(self, field_center, radius, **kwargs):
        """
        Query for asteroides visible in a field solving the 2-body equation from the asteroids'
        orbital parameters.
        
        Parameters:
                    field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
                    radius: astropy.units.Quantity with radius of field.
                    
        Optional Parameters:
                    epoch: astropy.time.Time with epoch for query. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    ast_idxs: Asteroids' indexes (starting from 0) from Query.asts. Ephemerides will be
                              calculated only for those asteroids.
                    
        Returns:
                  astropy.table with computed ephemerides.
        """
        # Getting all necessary parameters
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if 'epoch' in kwargs: epoch = kwargs['epoch']
        else: epoch = Time.now()
        if 'ast_idxs' in kwargs: ast_idxs = kwargs['ast_idxs']
        else: ast_idxs = np.arange(len(self.asts))
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd
        tr = Time(self.asts[0]['epoch'], format='iso', scale='tt').tdb.jd
        
        coords = []
        for idx in ast_idxs:
            ast = self.asts[idx]
            coords.append(apparent_coords_equatorial_heliocentric(tq, obsx, obsy, obsz, tr, ast))
            
        # Getting final coordinates
        coords = np.array(coords)
        [ra,dec,x,y,z,vx,vy,vz,geodist] = coords.T # deg, au
        dists = angular_separation(field_center.ra, field_center.dec, ra*u.deg, dec*u.deg)
        infield = dists.to('rad').value<radius.to('rad').value
        coords = coords[infield]; nums = self.numbers[tuple(ast_idxs),][infield]; nms = self.names[tuple(ast_idxs),][infield]
        
        #return coords, nums, nms
        
        t0 = Table((nums, nms), names=('number','name'))
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        us = (u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        cols = [column*us[i] for i, column in enumerate(coords.T)]
        t = QTable(cols, names=cnames)
        return hstack([t0, t])
    
    def query_nb(self, field_center, radius, **kwargs):
        """
        Query for asteroides visible in a field. Each asteroids is initialized using its orbital parameters
        and then it is integrated until requested epoch using an n-body simulation including the Sun,
        the planets, Pluto and the asteroid as a massless particle.
        
        Parameters:
                    field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
                    radius: astropy.units.Quantity with radius of field.
                    
        Optional Parameters:
                    epoch: astropy.time.Time with epoch for query. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    ast_idxs: Asteroids' indexes (starting from 0) from Query.asts. Ephemerides will be
                              calculated only for those asteroids.
                    
        Returns:
                  astropy.table with results.
        """
        # Getting all necessary parameters
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if 'epoch' in kwargs: epoch = kwargs['epoch']
        else: epoch = Time.now()
        if 'ast_idxs' in kwargs: ast_idxs = kwargs['ast_idxs']
        else: ast_idxs = np.arange(len(self.asts))
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd
        t0_py = Time(self.asts[0]['epoch'], format='iso', scale='tt')
        t0 = t0_py.tdb.jd
        ss_states = get_SS_states(t0_py, frame='helio')
        dt = tq-t0
        
        coords = []
        for idx in ast_idxs:
            ast = self.asts[idx]
            state_ast=state_equatorial_heliocentric(t0, t0, ast)
            coords.append(get_apparent_coordinates(state_ast, ss_states, obsx, obsy, obsz, dt))
            
        # Getting final coordinates
        coords = np.array(coords)
        [ra_geo, dec_geo, x ,y, z, vx,vy, vz, geodist] = coords.T # rad, au, au/d
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra_geo, dec_geo)
        infield = dists<radius.to('rad').value
        coords = coords[infield]; coords[:,0] = coords[:,0]*180/np.pi; coords[:,1] = coords[:,1]*180/np.pi
        nums = self.numbers[tuple(ast_idxs),][infield]; nms = self.names[tuple(ast_idxs),][infield]
        #return coords, nums, nms
        
        t0 = Table((nums, nms), names=('number','name'))
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        us = (u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        cols = [column*us[i] for i, column in enumerate(coords.T)]
        t = QTable(cols, names=cnames)
        return hstack([t0, t])                
