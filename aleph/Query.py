import numpy as np
from . import DataBase
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from astropy.time import Time
import sqlite3, time, warnings, os
import numpy as np
import multiprocessing as mp
from astropy.coordinates.angle_utilities import angular_separation
from .aleph_utils import *
from astropy.table import QTable, Table, hstack
import pandas as pd

def getint(s):
    try: return int(s)
    except: return np.nan

PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
IS_SQL_DB = os.path.exists(PACKAGE_PATH+'/aleph_states.db')

class Query:
    def __init__(self, service='MPC', open_db=IS_SQL_DB, dbfile=PACKAGE_PATH+'/aleph_states.db', **kwargs):
        """
        Initializing 'Query' enviroment. Orbital parameters are loaded
        and the sql database with asteroid's states is connected.

        Parameters:
                    service: Sets origin of orbital parameters. Can be 'Lowell' or 'MPC'.
                             Default to 'Lowell'.
        Other Parameters:
                    open_db
                    dbfile
                    filename: Str. Specifies file (and path) with orbital parameters to
                              read. This file must be in Lowell's or MPC's format.
        Attributes:
                    service: 'Lowell' or 'MPC'
                    asts: List of dictionaries containing all asteroids' information.
                    ref_epochs
                    ref_jdtdb
                    numbers
                    names
                    H
                    G
                    connection
                    cursor
                    asts_db
                    numbers_db
                    names_db
                    H_db
                    G_db
        """
        if (service=='MPC' and open_db==True):
            w = "You are using the 'MPC' survey AND opening the 'dbfile'='%s'" %dbfile
            w+= ". Only 'query_2b_cat()' and 'query_nb_cat()' will use 'MPCORB.DAT' data."
            warnings.warn(w)
        self.service = service
        if 'filename' in kwargs: alephs = DataBase.OrbParams(service, filename=kwargs['filename'])
        else: alephs = DataBase.OrbParams(service)
        self.asts = alephs.asts
        self.ref_epochs = alephs.dicc_epochs_M
        self.ref_jdtdb = {t:self.ref_epochs[t].tdb.jd for t in self.ref_epochs}
        nums = []; nms = []; hs = []; gs = []
        for ast in self.asts:
            nums.append(getint(ast['num'])); nms.append(ast['name'].strip())
            hs.append(ast['H']); gs.append(ast['G'])
        self.numbers = np.array(nums)
        self.names = np.array(nms)
        self.H = np.array(hs)
        self.G = np.array(gs)
        del alephs, nums, nms, hs, gs
        
        # Opening database
        if open_db:
            self.connection = sqlite3.connect(dbfile)
            self.cursor = self.connection.cursor()
            
            self.asts_db = pd.read_sql('SELECT * FROM asts_info', self.connection, index_col='name')
            
            self.numbers_db = np.array([getint(num) for num in self.asts_db.number.values])
            self.names_db = self.asts_db.index.values; self.H_db = self.asts_db.H.values; self.G_db = self.asts_db.G.values
    
    #def query(self, field_center, radius, **kwargs):
    #    """
    #    Query for asteroides visible in a field.
    #    
    #    Parameters:
    #                field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
    #                radius: astropy.units.Quantity with radius of field.
    #                
    #    Optional Parameters:
    #                epoch: astropy.time.Time with epoch for query. Default to 'now'.
    #                observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
     #               njobs: Number of parallel processes to calculate asteroids' apparent positions. Default: 1.
     #               method: 'db', '2B' or 'NB'. Default 'db'.
     #                      'db': Queries database for bodies near field and integrates them.
     #                      '2B': Calculates all bodies positions using 2-Bodies equations and retrieves bodies in
     #                            requested field.
     #                      'NB': Calculates all bodies positions using N-Bodies integration (Sun, planets, Pluto and 
     #                            massless asteroids) and retrieves bodies in requested field.
     #               db_query_radius: Radius to query database. Only usefull if 'method'='db'.
     #                                It should be greater than 'radius'. Default: 5 degs.
     #               ast_names: Asteroids' name from Query.asts. Ephemerides will be calculated only for those
     #                          asteroids.
     #               
     #   Returns:
     #            astropy.table with computed ephemerides.
     #             column 'number': Asteroid's number (from MPC).
     #             column 'name': Asteroid's name (from MPC).
     #             column 'ra': Asteroid's apparent right ascension (as seen by the observer).
     #             column 'dec': Asteroid's apparent declination (as seen by the observer).
     #             column 'x': Asteroid's heliocentric cartesian 'x' coordinate at light-emission time.
     #             column 'y': Asteroid's heliocentric cartesian 'y' coordinate at light-emission time.
     #             column 'z': Asteroid's heliocentric cartesian 'z' coordinate at light-emission time.
     #             column 'vx': Asteroid's heliocentric cartesian 'vx' velocity at light-emission time.
     #             column 'vy': Asteroid's heliocentric cartesian 'vy' velocity at light-emission time.
     #             column 'vz': Asteroid's heliocentric cartesian 'vz' velocity at light-emission time.
     #             column 'lightdist': Distance travel by light from asteroid's coordinate until observer.
     #   """
     #   # Selecting method
     #   if not 'method' in kwargs: return self.query_db(field_center, radius, **kwargs)
     #   else:
     #       if kwargs['method']=='db': return self.query_db(field_center, radius, **kwargs)
     #       elif kwargs['method']=='2B': return self.query_2b(field_center, radius, **kwargs)
     #       elif kwargs['method']=='NB': return self.query_nb(field_center, radius, **kwargs)
     #       else: raise ValueError("'method' must be one of 'db', '2B' or 'NB'")
    
    
    def query_2b_cat(self, field_center, radius, **kwargs):
        """
        Query for asteroides visible in a field solving the 2-body equation from the asteroids'
        orbital parameters stored in an 'astorb' or an 'MPCORB' catalogue.
        
        Parameters:
                    field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
                    radius: astropy.units.Quantity with radius of field.
                    
        Optional Parameters:
                    epoch: astropy.time.Time with epoch for query. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    njobs: Number of parallel processes to calculate asteroids' apparent positions. Default: 1.
                    ast_idxs: Asteroids' indexes (starting from 0) from Query.asts. Ephemerides will be
                              calculated only for those asteroids.
                    full: bool. If True, ephemeris are calculated for all asteroids (or those specified in
                          'ast_idxs') before selectig those in the requested field. If False, instantaneous
                          coordinates for all asteroids (or those specified in 'ast_idxs') are calculated first
                          and then apparent coordinates are calculated only for those bodies with instantaneos
                          coordinates inside a wider field. Default: False.
                    confidence_radius: Only when 'full=False'. Specifies the radius around requested 'field_center'
                                       to select the asteroids to calculate their apparent ephemeris. Default:
                                       radius+(3 degrees).

                    
        Returns:
                 astropy.table with computed ephemerides.
                  column 'number': Asteroid's number (from MPC).
                  column 'name': Asteroid's name (from MPC).
                  column 'ra': Asteroid's apparent right ascension (as seen by the observer).
                  column 'dec': Asteroid's apparent declination (as seen by the observer).
                  column 'x': Asteroid's heliocentric cartesian 'x' coordinate at light-emission time.
                  column 'y': Asteroid's heliocentric cartesian 'y' coordinate at light-emission time.
                  column 'z': Asteroid's heliocentric cartesian 'z' coordinate at light-emission time.
                  column 'vx': Asteroid's heliocentric cartesian 'vx' velocity at light-emission time.
                  column 'vy': Asteroid's heliocentric cartesian 'vy' velocity at light-emission time.
                  column 'vz': Asteroid's heliocentric cartesian 'vz' velocity at light-emission time.
                  column 'lightdist': Distance travel by light from asteroid's coordinate until observer.
        """
        # Getting all necessary parameters
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if 'epoch' in kwargs: epoch = kwargs['epoch']
        else: epoch = Time.now()
        if 'njobs' in kwargs: njobs = kwargs['njobs']
        else: njobs = 1
        if 'ast_idxs' in kwargs: ast_idxs = np.array(kwargs['ast_idxs'])
        else: ast_idxs = np.arange(len(self.asts))
        if 'full' in kwargs: full = kwargs['full']
        else: full = False
        if 'confidence_radius' in kwargs: confidence_radius = kwargs['confidence_radius'].to('deg').value
        else: confidence_radius = radius.to('deg').value + 3
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd

        if not full:
            # Instantaneous states are calculated for all asteroids.
            states = []
            for idx in ast_idxs:
                ast = self.asts[idx]
                tr = self.ref_jdtdb[ast['epoch']]
                r = pycutils.state_equatorial_heliocentric(tq, tr, ast['a'],ast['e'],ast['incl'],ast['Node'],ast['Arg_Peri'],ast['M_epoch'])
                states.append(r)
            [xh ,yh, zh] = np.array(states)[:,:3].T
            # Instantaneous ephemerides are calculated for all asteroids.
            xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
            ra_geo  = np.arctan2( ye, xe ); dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) )
            dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra_geo, dec_geo)
            inconffield = dists<(confidence_radius*np.pi/180)
            # Ephemeris are calculated for those bodies inside confidence field
            ast_idxs = ast_idxs[inconffield]
        
        # Apparent coordinates are calculated
        if njobs==1:
            coords = []
            for idx in ast_idxs:
                ast = self.asts[idx]
                tr = self.ref_jdtdb[ast['epoch']]
                args =  tq, obsx, obsy, obsz, tr, ast['a'], ast['e'], ast['incl'], ast['Node'], ast['Arg_Peri'], ast['M_epoch']
                coords.append(pycutils.apparent_coords_equatorial_heliocentric(*args))
        else:
            args = []
            for idx in ast_idxs:
                ast = self.asts[idx]
                tr = self.ref_jdtdb[ast['epoch']]
                args.append((tq, obsx, obsy, obsz, tr, ast['a'], ast['e'], ast['incl'], ast['Node'], ast['Arg_Peri'], ast['M_epoch']))
            pool = mp.Pool(processes=njobs)
            coords = pool.starmap(pycutils.apparent_coords_equatorial_heliocentric, args)
            pool.close()

        # Getting final coordinates
        coords = np.array(coords); ast_idxs_mask = tuple(ast_idxs),
        #[ra,dec,x,y,z,vx,vy,vz,geodist] = coords.T # deg, au
        ra = coords[:,0]; dec = coords[:,1]
        dists = angular_separation(field_center.ra, field_center.dec, ra*u.deg, dec*u.deg)
        infield = dists.to('rad').value<radius.to('rad').value
        coords = coords[infield]; nums = self.numbers[ast_idxs_mask][infield]; nms = self.names[ast_idxs_mask][infield]
        hs = self.H[ast_idxs_mask][infield]; gs = self.G[ast_idxs_mask][infield]
        
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]   # Obs vector from ast
        ast_sun = [-x, -y, -z]      # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)  # Phase angle (in rads)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-gs)*phase_func1 + gs*phase_func2
        Vs = hs + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)
        #return coords, nums, nms
        
        geo_ast = [xe, ye, ze]            # Obs vector from ast
        geo_sun = [-obsx, -obsy, -obsz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)  # Elongation angle (in rads)
        to_deg = 180/np.pi
        
        t0 = Table((nums, nms), names=('number','name'))
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        us = (u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        cols = [column*us[i] for i, column in enumerate(coords.T)]
        t = QTable(cols, names=cnames)
        t2= QTable([alpha*to_deg*u.deg, elong*to_deg*u.deg, Vs], names=['phase', 'elongation', 'V'])
        return hstack([t0, t, t2])

    
    def query_nb_cat(self, field_center, radius, **kwargs):
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
                    njobs: Number of parallel processes to calculate asteroids' apparent positions. Default: 1.
                    ast_idxs: Asteroids' indexes (starting from 0) from Query.asts. Ephemerides will be
                              calculated only for those asteroids.
                    
        Returns:
                 astropy.table with results.
                  column 'number': Asteroid's number (from MPC).
                  column 'name': Asteroid's name (from MPC).
                  column 'ra': Asteroid's apparent right ascension (as seen by the observer).
                  column 'dec': Asteroid's apparent declination (as seen by the observer).
                  column 'x': Asteroid's heliocentric cartesian 'x' coordinate at light-emission time.
                  column 'y': Asteroid's heliocentric cartesian 'y' coordinate at light-emission time.
                  column 'z': Asteroid's heliocentric cartesian 'z' coordinate at light-emission time.
                  column 'vx': Asteroid's heliocentric cartesian 'vx' velocity at light-emission time.
                  column 'vy': Asteroid's heliocentric cartesian 'vy' velocity at light-emission time.
                  column 'vz': Asteroid's heliocentric cartesian 'vz' velocity at light-emission time.
                  column 'lightdist': Distance travel by light from asteroid's coordinate until observer.
        """
        # Getting all necessary parameters
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if 'epoch' in kwargs: epoch = kwargs['epoch']
        else: epoch = Time.now()
        if 'njobs' in kwargs: njobs = kwargs['njobs']
        else: njobs = 1
        if 'ast_idxs' in kwargs: ast_idxs = kwargs['ast_idxs']
        else: ast_idxs = np.arange(len(self.asts))
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd

        # Getting Solar System's states in reference epochs
        set_ref_epochs = set([self.asts[idx]['epoch'] for idx in ast_idxs])
        dicc_ss_states = {epoch:get_SS_states(self.ref_epochs[epoch], frame='helio') for epoch in set_ref_epochs}
        
        if njobs==1:
            coords = []
            for idx in ast_idxs:
                ast = self.asts[idx]
                tr = self.ref_jdtdb[ast['epoch']]; dt = tq-tr
                ss_states = dicc_ss_states[ast['epoch']]
                args =  tr, tr, ast['a'], ast['e'], ast['incl'], ast['Node'], ast['Arg_Peri'], ast['M_epoch']
                state_ast=pycutils.state_equatorial_heliocentric(*args)
                coords.append(get_apparent_coordinates(state_ast, ss_states, obsx, obsy, obsz, dt))
        else:
            arguments = []
            for idx in ast_idxs:
                ast = self.asts[idx]
                tr = self.ref_jdtdb[ast['epoch']]; dt = tq-tr
                ss_states = dicc_ss_states[ast['epoch']]
                args =  tr, tr, ast['a'], ast['e'], ast['incl'], ast['Node'], ast['Arg_Peri'], ast['M_epoch']
                state = pycutils.state_equatorial_heliocentric(*args)
                arguments.append((state, ss_states, obsx, obsy, obsz, dt))
            pool = mp.Pool(processes=njobs)
            coords = pool.starmap(get_apparent_coordinates, arguments)
            pool.close()
            
        # Getting final coordinates
        coords = np.array(coords); ast_idxs_mask = tuple(ast_idxs),
        #[ra_geo, dec_geo, x ,y, z, vx,vy, vz, geodist] = coords.T # rad, au, au/d
        ra_geo = coords[:,0]; dec_geo = coords[:,1]
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra_geo, dec_geo)
        infield = dists<radius.to('rad').value
        coords = coords[infield]; coords[:,0] = coords[:,0]*180/np.pi; coords[:,1] = coords[:,1]*180/np.pi
        nums = self.numbers[ast_idxs_mask][infield]; nms = self.names[ast_idxs_mask][infield]
        hs = self.H[ast_idxs_mask][infield]; gs = self.G[ast_idxs_mask][infield]
        
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]   # Obs vector from ast
        ast_sun = [-x, -y, -z]      # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-gs)*phase_func1 + gs*phase_func2
        Vs = hs + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)
        #return coords, nums, nms
        
        geo_ast = [xe, ye, ze]            # Obs vector from ast
        geo_sun = [-obsx, -obsy, -obsz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)  # Elongation angle (in rads)
        to_deg = 180/np.pi
        
        t0 = Table((nums, nms), names=('number','name'))
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        us = (u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        cols = [column*us[i] for i, column in enumerate(coords.T)]
        t = QTable(cols, names=cnames)
        t2= QTable([alpha*to_deg*u.deg, elong*to_deg*u.deg, Vs], names=['phase', 'elongation', 'V'])
        return hstack([t0, t, t2]) 

    def query_mixed_cat(self, field_center, radius, **kwargs):
        """
        Query for asteroides visible in a field. Asteroids ephemeris are pre-calculated using
        their orbital parameters. The ephemeris of tose asteroids that are within a 'confidence
        radius' (larger than the query's radius) are then integrated to increase their accuracy.
        Ephemeris integrations consider asteroids as massless particles under the influence of
        the Sun, the planets and Pluto.
    
        Parameters:
                    field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
                    radius: astropy.units.Quantity with radius of field.
                    
        Optional Parameters:
                    epoch: astropy.time.Time with epoch for query. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    njobs: Number of parallel processes to calculate asteroids' apparent positions. Default: 1.
                    confidence_radius: The 2-Body ephemeris of the asteroids within this radius (around
                                       'field_center') are selected to be integrated. It should be greater
                                       than 'radius'. Default: 5 degs.
                    ast_idxs: Asteroids' indexes (starting from 0) from Query.asts. Ephemerides will be
                              calculated only for those asteroids. It works as a mask of the asteroids data.
                    
        Returns:
                 astropy.table with computed ephemerides:
                  column 'number': Asteroid's number (from MPC).
                  column 'name': Asteroid's name (from MPC).
                  column 'ra': Asteroid's apparent right ascension (as seen by the observer).
                  column 'dec': Asteroid's apparent declination (as seen by the observer).
                  column 'x': Asteroid's heliocentric cartesian 'x' coordinate at light-emission time.
                  column 'y': Asteroid's heliocentric cartesian 'y' coordinate at light-emission time.
                  column 'z': Asteroid's heliocentric cartesian 'z' coordinate at light-emission time.
                  column 'vx': Asteroid's heliocentric cartesian 'vx' velocity at light-emission time.
                  column 'vy': Asteroid's heliocentric cartesian 'vy' velocity at light-emission time.
                  column 'vz': Asteroid's heliocentric cartesian 'vz' velocity at light-emission time.
                  column 'lightdist': Distance travel by light from asteroid's coordinate until observer.
        """
        # Getting all necessary parameters
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if 'epoch' in kwargs: epoch = kwargs['epoch']
        else: epoch = Time.now()
        if 'njobs' in kwargs: njobs = kwargs['njobs']
        else: njobs = 1
        if 'ast_idxs' in kwargs: ast_idxs = np.array(kwargs['ast_idxs'])
        else: ast_idxs = np.arange(len(self.asts))
        if 'confidence_radius' in kwargs: confidence_radius = kwargs['confidence_radius'].to('deg').value
        else: confidence_radius = radius.to('deg').value + 3
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd

        # Instantaneous states are calculated for all asteroids.
        states = []
        for idx in ast_idxs:
            ast = self.asts[idx]
            tr = self.ref_jdtdb[ast['epoch']]
            r = pycutils.state_equatorial_heliocentric(tq, tr, ast['a'],ast['e'],ast['incl'],ast['Node'],ast['Arg_Peri'],ast['M_epoch'])
            states.append(r)
        [xh ,yh, zh] = np.array(states)[:,:3].T
        # Instantaneous ephemerides are calculated for all asteroids.
        xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
        ra_geo  = np.arctan2( ye, xe ); dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) )
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra_geo, dec_geo)
        inconffield = dists<(confidence_radius*np.pi/180)
        # Ephemeris are calculated for those bodies inside confidence field
        ast_idxs = ast_idxs[inconffield]
        
        # Getting Solar System's states in pertinent reference epochs
        set_ref_epochs = set([self.asts[idx]['epoch'] for idx in ast_idxs])
        dicc_ss_states = {epoch:get_SS_states(self.ref_epochs[epoch], frame='helio') for epoch in set_ref_epochs}
        
        # Apparent coordinates are calculated integrating the N-body problem
        if njobs==1:
            coords = []
            for idx in ast_idxs:
                ast = self.asts[idx]
                tr = self.ref_jdtdb[ast['epoch']]; dt = tq-tr
                ss_states = dicc_ss_states[ast['epoch']]
                args =  tr, tr, ast['a'], ast['e'], ast['incl'], ast['Node'], ast['Arg_Peri'], ast['M_epoch']
                state_ast=pycutils.state_equatorial_heliocentric(*args)
                coords.append(get_apparent_coordinates(state_ast, ss_states, obsx, obsy, obsz, dt))
        else:
            arguments = []
            for idx in ast_idxs:
                ast = self.asts[idx]
                tr = self.ref_jdtdb[ast['epoch']]; dt = tq-tr
                ss_states = dicc_ss_states[ast['epoch']]
                args =  tr, tr, ast['a'], ast['e'], ast['incl'], ast['Node'], ast['Arg_Peri'], ast['M_epoch']
                state = pycutils.state_equatorial_heliocentric(*args)
                arguments.append((state, ss_states, obsx, obsy, obsz, dt))
            pool = mp.Pool(processes=njobs)
            coords = pool.starmap(get_apparent_coordinates, arguments)
            pool.close()

        # Getting final coordinates
        coords = np.array(coords); ast_idxs_mask = tuple(ast_idxs),
        #[ra,dec,x,y,z,vx,vy,vz,geodist] = coords.T # deg, au
        ra = coords[:,0]; dec = coords[:,1]
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra, dec)
        infield = dists<radius.to('rad').value
        coords = coords[infield]
        coords[:,0] = coords[:,0]*180/np.pi; coords[:,1] = coords[:,1]*180/np.pi
        nums = self.numbers[ast_idxs_mask][infield]
        nms = self.names[ast_idxs_mask][infield]
        hs = self.H[ast_idxs_mask][infield]
        gs = self.G[ast_idxs_mask][infield]
        
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]   # Obs vector from ast
        ast_sun = [-x, -y, -z]      # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-gs)*phase_func1 + gs*phase_func2
        Vs = hs + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)
        #return coords, nums, nms
        
        geo_ast = [xe, ye, ze]         # Obs vector from ast
        geo_sun = [-obsx, -obsy, -obsz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)  # Elongation angle (in rads)
        to_deg = 180/np.pi
        
        t0 = Table((nums, nms), names=('number','name'))
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        us = (u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        cols = [column*us[i] for i, column in enumerate(coords.T)]
        t = QTable(cols, names=cnames)
        t2= QTable([alpha*to_deg*u.deg, elong*to_deg*u.deg, Vs], names=['phase', 'elongation', 'V'])
    
        return hstack([t0, t, t2])  

    def query_2b_db(self, field_center, radius, **kwargs):
        """
        Query for asteroides visible in a field solving the 2-body equation from the asteroids'
        orbital parameters saved in the SQL database.
        
        Parameters:
                    field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
                    radius: astropy.units.Quantity with radius of field.
                    
        Optional Parameters:
                    epoch: astropy.time.Time with epoch for query. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    njobs: Number of parallel processes to calculate asteroids' apparent positions. Default: 1.
                    ast_names: Asteroids' names (from Query.asts_db). Ephemerides will be calculated only for
                               those asteroids.
                    
        Returns:
                 astropy.table with computed ephemerides.
                  column 'number': Asteroid's number (from MPC).
                  column 'name': Asteroid's name (from MPC).
                  column 'ra': Asteroid's apparent right ascension (as seen by the observer).
                  column 'dec': Asteroid's apparent declination (as seen by the observer).
                  column 'x': Asteroid's heliocentric cartesian 'x' coordinate at light-emission time.
                  column 'y': Asteroid's heliocentric cartesian 'y' coordinate at light-emission time.
                  column 'z': Asteroid's heliocentric cartesian 'z' coordinate at light-emission time.
                  column 'vx': Asteroid's heliocentric cartesian 'vx' velocity at light-emission time.
                  column 'vy': Asteroid's heliocentric cartesian 'vy' velocity at light-emission time.
                  column 'vz': Asteroid's heliocentric cartesian 'vz' velocity at light-emission time.
                  column 'lightdist': Distance travel by light from asteroid's coordinate until observer.
        """
        if self.service == 'MPC': raise ValueError("'service' must be 'Lowell'.")
        # Getting all necessary parameters
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if 'epoch' in kwargs: epoch = kwargs['epoch']
        else: epoch = Time.now()
        if 'njobs' in kwargs: njobs = kwargs['njobs']
        else: njobs = 1
        if 'ast_names' in kwargs:
            if isinstance(kwargs['ast_names'], (list,np.ndarray)): ast_names = list(kwargs['ast_names'])
            elif isinstance(kwargs['ast_names'], (str)): ast_names = [kwargs['ast_names']]
            else:
                print("ERROR in 'ast_names'"); return
            allasts = False
        else: allasts = True
        c = self.cursor; conn = self.connection
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd

        # Looking for parameters table at nearest time
        c.execute("SELECT table_title, time, ABS(time-?) FROM param_times ORDER BY ABS(time-?)", (tq,tq))
        (table, tr, terr) = c.fetchone()

        # Loading orbital parameters from database
        if allasts: qw = 'SELECT * FROM %s' %table
        else: qw = 'SELECT * FROM %s WHERE name IN (%s)' %(table,str(ast_names)[1:-1])
        tb = pd.read_sql(qw, conn, index_col='name')
        params = tb.values
        
        # Getting apparent coordinates
        args0 = tq, obsx, obsy, obsz, tr
        if njobs==1:
            coords = np.array([pycutils.apparent_coords_equatorial_heliocentric(*args0,*args1) for args1 in params])
        else:
            args = [(*args0,*args1) for args1 in params]
            pool = mp.Pool(processes=njobs)
            coords = np.array(pool.starmap(pycutils.apparent_coords_equatorial_heliocentric, args))
            pool.close()
            
        # Getting final coordinates
        ra = coords[:,0]; dec = coords[:,1]
        dists = angular_separation(field_center.ra, field_center.dec, ra*u.deg, dec*u.deg)
        infield = dists.to('rad').value<radius.to('rad').value
        coords = coords[infield]
        nms = tb.index.values[infield]; info = self.asts_db.loc[nms]
        nums = info.number.values; hs = info.H.values; gs = info.G.values
        
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]  # Obs vector from ast
        ast_sun = [-x, -y, -z]     # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)  # Phase angle (in rads)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-gs)*phase_func1 + gs*phase_func2
        Vs = hs + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)
        
        geo_ast = [xe, ye, ze]            # Obs vector from ast
        geo_sun = [-obsx, -obsy, -obsz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)  # Elongation angle (in rads)
        to_deg = 180/np.pi
        
        t0 = Table((nums, nms), names=('number','name'))
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        us = (u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        cols = [column*us[i] for i, column in enumerate(coords.T)]
        t = QTable(cols, names=cnames)
        t2= QTable([alpha*to_deg*u.deg, elong*to_deg*u.deg, Vs], names=['phase', 'elongation', 'V'])
        return hstack([t0, t, t2])


    def query_nb_db(self, field_center, radius, **kwargs):
        """
        Query for asteroides visible in a field. Asteoids ephemerides are calculated integrating them
        as massless particles under the influence of the Sun, the planets and Pluto. 
        Asteroids integrated are only those from nearby the requested field at the closest epoch 
        stored in the database.
    
        Parameters:
                    field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
                    radius: astropy.units.Quantity with radius of field.
                    
        Optional Parameters:
                    epoch: astropy.time.Time with epoch for query. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    njobs: Number of parallel processes to calculate asteroids' apparent positions. Default: 1.
                    confidence_radius: Radius to query database. It should be greater than 'radius'. 
                                       Default: 5 degs.
                    ast_names: Asteroids' names (from Query.asts_db). Ephemerides will be calculated only for
                               those asteroids.
                    
        Returns:
                 astropy.table with computed ephemerides:
                  column 'number': Asteroid's number (from MPC).
                  column 'name': Asteroid's name (from MPC).
                  column 'ra': Asteroid's apparent right ascension (as seen by the observer).
                  column 'dec': Asteroid's apparent declination (as seen by the observer).
                  column 'x': Asteroid's heliocentric cartesian 'x' coordinate at light-emission time.
                  column 'y': Asteroid's heliocentric cartesian 'y' coordinate at light-emission time.
                  column 'z': Asteroid's heliocentric cartesian 'z' coordinate at light-emission time.
                  column 'vx': Asteroid's heliocentric cartesian 'vx' velocity at light-emission time.
                  column 'vy': Asteroid's heliocentric cartesian 'vy' velocity at light-emission time.
                  column 'vz': Asteroid's heliocentric cartesian 'vz' velocity at light-emission time.
                  column 'lightdist': Distance travel by light from asteroid's coordinate until observer.
        """
        if self.service == 'MPC': raise ValueError("'service' must be 'Lowell'.")
        # Getting all necessary parameters
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if 'epoch' in kwargs: epoch = kwargs['epoch']
        else: epoch = Time.now()
        if 'confidence_radius' in kwargs: confidence_radius = kwargs['confidence_radius'].to('deg').value
        else: confidence_radius = radius.to('deg').value + 5
        if 'njobs' in kwargs: njobs = kwargs['njobs']
        else: njobs = 1
        if 'ast_names' in kwargs:
            if isinstance(kwargs['ast_names'], (list,np.ndarray)): ast_names = list(kwargs['ast_names'])
            elif isinstance(kwargs['ast_names'], (str)): ast_names = [kwargs['ast_names']]
            else:
                print("ERROR in 'ast_names'"); return
            allasts = False
        else: allasts = True
        c = self.cursor; conn = self.connection
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd
        
        # Looking for table at nearest time
        #c.execute("select table_title, time from times where abs(time-%f)<.25" %(tq))
        #[(table, t0)] = c.fetchall()
        c.execute("SELECT table_title, time, ABS(time-?) FROM state_times ORDER BY ABS(time-?)", (tq,tq))
        (table, t0, terr) = c.fetchone()
        if terr>0.25: 
            w = "Queried epoch is %.2f days away from closest epoch in the ephemeris database. " %terr
            w+= "User may want to increase 'db_query_radius' value."
            warnings.warn(w)
        
        # Looking for bodies near field
        ## Bodies' cartesian coordinates
        
        # Loading states from database
        if allasts: qw = 'SELECT * FROM %s' %table
        else: qw = 'SELECT * FROM %s WHERE name IN (%s)' %(table,str(ast_names)[1:-1])
        tb = pd.read_sql(qw, conn, index_col='name')
        states = tb.values; nms = tb.index.values
        
        x,y,z=states.T[:3]
        ## Bodies' equatorial coordinates from observer
        xe = x-obsx; ye = y-obsy; ze = z-obsz
        ra_geo  = np.arctan2( ye, xe ) * 180/np.pi; dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) ) * 180/np.pi
        coords_db = SkyCoord(ra_geo, dec_geo, unit='deg')
        dists = field_center.separation(coords_db)
        near = dists.deg < confidence_radius
        states = states[near]; nms = nms[near]
        
        # Calculating bodies' positions at requested epoch
        ## Solar system states
        t0_py = Time(t0, format='jd', scale='tdb')
        ss_states = get_SS_states(t0_py, frame='helio')
    
        dt = tq-t0
        args_ = ss_states, obsx, obsy, obsz, dt
        if njobs == 1:
            ephs = []
            for i, state in enumerate(states):
                ephs.append(get_apparent_coordinates(state, *args_))
            ephs = np.array(ephs)
        else:
            pool = mp.Pool(processes=njobs)
            arguments = []
            for state in states:
                arguments.append((state, *args_))
            ephs = np.array(pool.starmap(get_apparent_coordinates, arguments))
            pool.close()
        
        # Getting final coordinates
        ephs = np.array(ephs)
        ra = ephs[:,0]; dec = ephs[:,1]
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra, dec)
        infield = dists<radius.to('rad').value
        ephs = ephs[infield]; ephs[:,:2] = ephs[:,:2]*180/np.pi
        nms = nms[infield]; info = self.asts_db.loc[nms]
        nums = info.number.values; hs = info.H.values; gs = info.G.values
        
        [x ,y, z] = ephs[:,2:5].T; geodist = ephs[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]  # Obs vector from ast
        ast_sun = [-x, -y, -z]     # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)  # Phase angle (in rads)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-gs)*phase_func1 + gs*phase_func2
        Vs = hs + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)
        #return coords, nums, nms
        
        geo_ast = [xe, ye, ze]            # Obs vector from ast
        geo_sun = [-obsx, -obsy, -obsz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)  # Elongation angle (in rads)
        to_deg = 180/np.pi
        
        t0 = Table((nums, nms), names=('number','name'))
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        us = (u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        cols = [column*us[i] for i, column in enumerate(ephs.T)]
        t = QTable(cols, names=cnames)
        t2= QTable([alpha*to_deg*u.deg, elong*to_deg*u.deg, Vs], names=['phase', 'elongation', 'V'])
        return hstack([t0, t, t2])
    
    def query_mixed_db(self, field_center, radius, **kwargs):
        """
        Query for asteroides visible in a field. Asteroids ephemeris are pre-calculated using
        their orbital parameters. The ephemeris of tose asteroids that are within a 'confidence
        radius' (larger than the query's radius) are then integrated to increase their accuracy.
        Ephemeris integrations consider asteroids as massless particles under the influence of
        the Sun, the planets and Pluto.
    
        Parameters:
                    field_center: astropy.coordinates.SkyCoord element with center coordinate of field.
                    radius: astropy.units.Quantity with radius of field.
                    
        Optional Parameters:
                    epoch: astropy.time.Time with epoch for query. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    njobs: Number of parallel processes to calculate asteroids' apparent positions. Default: 1.
                    confidence_radius: The 2-Body ephemeris of the asteroids within this radius (around
                                       'field_center') are selected to be integrated. It should be greater
                                       than 'radius'. Default: 5 degs.
                    ast_names: Asteroids' names (from Query.asts_db). Ephemerides will be calculated only for
                               those asteroids.
                    
        Returns:
                 astropy.table with computed ephemerides:
                  column 'number': Asteroid's number (from MPC).
                  column 'name': Asteroid's name (from MPC).
                  column 'ra': Asteroid's apparent right ascension (as seen by the observer).
                  column 'dec': Asteroid's apparent declination (as seen by the observer).
                  column 'x': Asteroid's heliocentric cartesian 'x' coordinate at light-emission time.
                  column 'y': Asteroid's heliocentric cartesian 'y' coordinate at light-emission time.
                  column 'z': Asteroid's heliocentric cartesian 'z' coordinate at light-emission time.
                  column 'vx': Asteroid's heliocentric cartesian 'vx' velocity at light-emission time.
                  column 'vy': Asteroid's heliocentric cartesian 'vy' velocity at light-emission time.
                  column 'vz': Asteroid's heliocentric cartesian 'vz' velocity at light-emission time.
                  column 'lightdist': Distance travel by light from asteroid's coordinate until observer.
        """
        if self.service == 'MPC': raise ValueError("'service' must be 'Lowell'.")
        # Getting all necessary parameters
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if 'epoch' in kwargs: epoch = kwargs['epoch']
        else: epoch = Time.now()
        if 'njobs' in kwargs: njobs = kwargs['njobs']
        else: njobs = 1
        if 'confidence_radius' in kwargs: confidence_radius = kwargs['confidence_radius'].to('deg').value
        else: confidence_radius = radius.to('deg').value + 3
        if 'ast_names' in kwargs:
            if isinstance(kwargs['ast_names'], (list,np.ndarray)): ast_names = list(kwargs['ast_names'])
            elif isinstance(kwargs['ast_names'], (str)): ast_names = [kwargs['ast_names']]
            else:
                print("ERROR in 'ast_names'"); return
            allasts = False
        else: allasts = True
        c = self.cursor; conn = self.connection
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd
        
        # Looking for orbital parameters table and states table at nearest time
        c.execute("SELECT table_title, time, ABS(time-?) FROM param_times ORDER BY ABS(time-?)", (tq,tq))
        (table_params, tr, terr) = c.fetchone()
        table_states = 'state_'+table_params.split('_')[1]
        
        # Loading orbital parameters from database
        if allasts: qw = "SELECT * FROM %s" %table_params
        else: qw = "SELECT * FROM %s WHERE name IN (%s)" %(table_params,str(ast_names)[1:-1])
        tb_params = pd.read_sql(qw, conn, index_col='name')
        params = tb_params.values; nms = tb_params.index.values
        
        # Getting apparent coordinates using 2-body solution
        #args0 = tq, obsx, obsy, obsz, tr
        #coords = np.array([pycutils.apparent_coords_equatorial_heliocentric(*args0,*args1) for args1 in params])
        #ra = coords[:,0]; dec = coords[:,1]
        
        # Getting instantaneous states using 2-body solution (about 3s faster than getting apparent coordinates)
        args0 = tq, tr
        states_ = np.array([pycutils.state_equatorial_heliocentric(*args0,*args1) for args1 in params])
        [xh ,yh, zh] = np.array(states_)[:,:3].T
        # Instantaneous ephemerides are calculated for all asteroids.
        xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
        ra  = np.arctan2( ye, xe ); dec = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) )
        
        # Selecting bodies in confidence radius for n-body integration
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra, dec)
        inconffield = dists<(confidence_radius*np.pi/180)
        nms = nms[inconffield]
        
        # Loading states from database
        if allasts: qw = 'SELECT * FROM %s' %table_states
        else: qw = 'SELECT * FROM %s WHERE name IN (%s)' %(table_states,str(list(nms))[1:-1])
        tb_states = pd.read_sql(qw, conn, index_col='name').loc[nms]
        states = tb_states.values
        
        # Calculating bodies' positions at requested epoch
        ## Solar system states
        tr_py = Time(tr, format='jd', scale='tdb')
        ss_states = get_SS_states(tr_py, frame='helio')
        
        dt = tq-tr
        args0 = ss_states, obsx, obsy, obsz, dt
        if njobs == 1:
            ephs = np.array([get_apparent_coordinates(state, *args0) for state in states])
        else:
            pool = mp.Pool(processes=njobs)
            arguments = [(state, *args0) for state in states]
            ephs = np.array(pool.starmap(get_apparent_coordinates, arguments))
            pool.close()
        
        # Getting final coordinates
        ra = ephs[:,0]; dec = ephs[:,1]
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra, dec)
        infield = dists<radius.to('rad').value
        ephs = ephs[infield]; ephs[:,:2] = ephs[:,:2]*180/np.pi
        nms = nms[infield]; info = self.asts_db.loc[nms]
        nums = info.number.values; hs = info.H.values; gs = info.G.values
        
        [x, y, z] = ephs[:,2:5].T; geodist = ephs[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]  # Obs vector from ast
        ast_sun = [-x, -y, -z]     # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)  # Phase angle (in rads)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-gs)*phase_func1 + gs*phase_func2
        Vs = hs + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)
        #return coords, nums, nms
        
        geo_ast = [xe, ye, ze]            # Obs vector from ast
        geo_sun = [-obsx, -obsy, -obsz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)  # Elongation angle (in rads)
        to_deg = 180/np.pi
        
        t0 = Table((nums, nms), names=('number','name'))
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        us = (u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        cols = [column*us[i] for i, column in enumerate(ephs.T)]
        t = QTable(cols, names=cnames)
        t2= QTable([alpha*to_deg*u.deg, elong*to_deg*u.deg, Vs], names=['phase', 'elongation', 'V'])
        return hstack([t0, t, t2]) 
    
    def close(self): self.connection.close()
        
    def asteph(self, ast, **kwargs):
        """
        Returns ephemerides from specified asteroid at given epochs.
        Ephemerides are calculated performing an n-body simulation that considers
        the sun, all 8 planets and the asteroid as a massless particle, starting from
        a time where all states are known and integrating all bodies position until
        requested time (correcting by light aberration).
    
        Parameters:
                    ast: Name of the asteroid (str) or it index number (int). If a string is given,
                         a perfect match of the name is looked for at the Lowell's database. If a number
                         is given, it is considered as the index (starting from zero) of the asteroid's
                         in the Lowell' database (so, for numbered bodies, if you want '(1) Ceres',
                         the index is zero).
                    
        Optional Parameters:
                    epochs: astropy.time.Time with epoch(s) to calculate ephemerides. Default to 'now'.
                    observer: astropy.coordinates.EarthLocation with observer's position. Default to geocenter.
                    
        Returns:
                 astropy.table with computed ephemerides:
                  column 'epoch_jd': Epoch in julian days (in UTC).
                  column 'ra': Asteroid's apparent right ascension (as seen by the observer).
                  column 'dec': Asteroid's apparent declination (as seen by the observer).
                  column 'x': Asteroid's heliocentric cartesian 'x' coordinate at light-emission time.
                  column 'y': Asteroid's heliocentric cartesian 'y' coordinate at light-emission time.
                  column 'z': Asteroid's heliocentric cartesian 'z' coordinate at light-emission time.
                  column 'vx': Asteroid's heliocentric cartesian 'vx' velocity at light-emission time.
                  column 'vy': Asteroid's heliocentric cartesian 'vy' velocity at light-emission time.
                  column 'vz': Asteroid's heliocentric cartesian 'vz' velocity at light-emission time.
                  column 'lightdist': Distance travel by light from asteroid's coordinate until observer.
                  column 'phase': Angle between the Sun and the Observer from the asteroid.
                  column 'elongation': Angle between the Sun and the asteroid from the Obsever.
                  column 'V':  Apparent magnitude as seen by the observer in V band (calculated from the
                               asteroid's absolute magnitude in Lowell's database).
        """
        if self.service == 'MPC': raise ValueError("'service' must be 'Lowell'.")
        # Getting all necessary parameters
        if isinstance(ast, (int,np.int64)): ast_name = self.asts_db.iloc[ast].name
        if isinstance(ast, str): ast_name = ast
        
        if 'observer' in kwargs: observer = kwargs['observer']
        else: observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        
        if 'epochs' in kwargs:
            if not isinstance(kwargs['epochs'], Time):
                print('ERROR'); return
            if isinstance(kwargs['epochs'].jd, np.ndarray):
                epochs0 = kwargs['epochs']
            else: epochs0 = Time([kwargs['epochs']])
        else: epochs0 = Time([Time.now().jd], format='jd')
        
        c = self.cursor
        iepochs = np.argsort(epochs0.jd); epochs = epochs0[iepochs] # Epochs are sorted to speed up integration
        tqs  = epochs.tdb.jd
        tmin = min(tqs); tmax = min(tqs)
        
        # Looking for table at nearest time
        c.execute("SELECT table_title, time, ABS(time-?) FROM state_times ORDER BY ABS(time-?)", (tmin,tmin))
        (table_min, t0_min, terr_min) = c.fetchone()
        c.execute("SELECT table_title, time, ABS(time-?) FROM state_times ORDER BY ABS(time-?)", (tmax,tmax))
        (table_max, t0_max, terr_max) = c.fetchone()
        if terr_min < terr_max:
            t0 = t0_min; table = table_min; dir_integ = 1
        else:
            t0 = t0_max; table = table_max; dir_integ = -1; tqs = tqs[::-1]; epochs = epochs[::-1]; iepochs = iepochs[::-1]
        
        # Observatory coordinates
        [geopos,geovel] = observer.get_gcrs_posvel(epochs)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epochs)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        
        # Asteroid's state at nearest time
        c.execute('SELECT x,y,z,vx,vy,vz FROM %s WHERE name="%s"' %(table,ast_name))
        state_ast = c.fetchone()
        
        # Calculating bodies' positions at requested epoch
        ## Solar system states
        t0_py = Time(t0, format='jd', scale='tdb')
        ss_states = get_SS_states(t0_py, frame='helio')
            
        # ASTEROID'S INTEGRATION
        ## Creating simulation
        sim = get_SS_sim(ss_states); sim.N_active = sim.N # Solar system
        sim.add(x=state_ast[0], y=state_ast[1], z=state_ast[2], vx=state_ast[3], vy=state_ast[4], vz=state_ast[5]) # Asteroid
        sun = sim.particles[0]; ast = sim.particles[-1]; sim.dt = 0.1*dir_integ
        
        #return tqs, t0
        delta_tdbs = tqs - t0
        xhs, yhs, zhs, vxhs, vyhs, vzhs, lightdist = [], [], [], [], [], [], []
        for i, delta_tdb in enumerate(delta_tdbs):
            # Integration until observation time
            sim.integrate(delta_tdb)
            # Correcting by light's travel time
            gdist = np.sqrt((ast.x-sun.x-obsx[i])**2+(ast.y-sun.y-obsy[i])**2+(ast.z-sun.z-obsz[i])**2)
            lighttime = gdist*au_to_d
            t3 = delta_tdb - lighttime
            niter=0
            while niter<31:
                niter += 1; t2=t3
                sim.integrate(t2)
                gdist = np.sqrt((ast.x-sun.x-obsx[i])**2+(ast.y-sun.y-obsy[i])**2+(ast.z-sun.z-obsz[i])**2)
                lighttime = gdist*au_to_d
                t3 = delta_tdb - lighttime
                if abs (t3 - t2) < 1.0e-10: break
            xhs.append(ast.x-sun.x)   ; yhs.append(ast.y-sun.y)   ; zhs.append(ast.z-sun.z)
            vxhs.append(ast.vx-sun.vx); vyhs.append(ast.vy-sun.vy); vzhs.append(ast.vz-sun.vz)
            lightdist.append(lighttime/au_to_d)  # In au
            
        # Ordering data to original epochs' order
        N = len(epochs)
        obx,oby,obz,x,y,z = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
        vx,vy,vz,geodist = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
        obx[iepochs]=obsx; oby[iepochs]=obsy; obz[iepochs]=obsz; x[iepochs]=xhs; y[iepochs]=yhs; z[iepochs]=zhs
        vx[iepochs]=vxhs; vy[iepochs]=vyhs; vz[iepochs]=vzhs; geodist[iepochs]=lightdist
        
        # Getting final data
        c.execute('SELECT number,H,G FROM asts_info WHERE name="%s"' %(ast_name))
        num,h,g = c.fetchone()
        
        xe = x-obx; ye = y-oby; ze = z-obz
        ra  = np.arctan2( ye, xe ); dec = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) ) # In rads
        neg_ra = ra < 0; ra[neg_ra] += pycutils.TWOPI # In rads
        to_deg = 180/np.pi; ra *= to_deg; dec *= to_deg
        
        ast_geo = [-xe, -ye, -ze]  # Obs vector from ast
        ast_sun = [-x, -y, -z]     # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)  # Phase angle (in rads)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-g)*phase_func1 + g*phase_func2
        Vs = h + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)
        
        geo_ast = [xe, ye, ze]         # Obs vector from ast
        geo_sun = [-obx, -oby, -obz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)*to_deg  # Elongation angle (in degs)
        
        #t0 = Table((nums, nms), names=('number','name'))
        deg, au, au_d = u.deg, u.au, u.au/u.day
        cnames=('epoch_jd', 'ra','dec','x','y','z','vx','vy','vz','lightdist','phase','elongation', 'V')
        cols = [epochs0.jd, ra*deg, dec*deg, x*au, y*au, z*au, vx*au_d, vy*au_d, vz*au_d, geodist*au,
                alpha*to_deg*deg, elong*deg, Vs]
        t = QTable(cols, names=cnames)
        
        return t
