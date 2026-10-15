import sqlite3, warnings, os
from pathlib import Path
import multiprocessing as mp
import numpy as np
import pandas as pd
from collections.abc import Callable
from astropy.coordinates import SkyCoord, EarthLocation, angular_separation
from astropy.time import Time
from astropy.table import QTable, Table, hstack
from astropy import units as u
from .database import OrbParams
from .aleph_utils import *
from joblib import Parallel, delayed

def getint(s):
    try: return int(s)
    except: return np.nan

PACKAGE_PATH = Path(os.path.abspath(os.path.dirname(__file__)))
Parser = Callable[[str | Path], pd.DataFrame]

class Query:
    def __init__(
        self,
        filename: str|Path|None = None,
        parsing: str|Parser = 'MPC',
        open_sql_database: bool = False,
        sqlfile: str|Path|None = None,
        update: bool = False,
        n_jobs: int = 1,
        ):
        """
        Solar system minor bodies ephemeries calculator.
        It can use the MPC minor bodies orbital parameters catalogue (by default)
        or the Lowell's minor bodies orbital parameters catalogue.

        Parameters:
            filename : str, default=None
                Specifies file (and path) with orbital parameters to read.
            parsing : str, default='MPC'
                Sets parser to use. Can be 'MPC' or 'Lowell'.
            open_sql_database: bool, default=True if exists, False elsewhere
                If True, opens an SQL database with asteroids' states.
                If not specified, it looks for the existence of the default
                database file, 'aleph_states.db', in the package's directory.
            sqlfile: str, default='{package_path}/aleph_states.db'
                Path to the database file. If not provided, it looks for the
                default catalogue file name in the package's directory.
                
        Attributes:
            filename : str
                Path to used catalogue file.
            asts : list[dict]
                List with asteroids' information.
            ref_epochs: dict[str, astropy.time.Time]
                Dictionary with all reference epochs in the catalogue.
            connection: sqlite3.Connection
                Connection to the database (if database exists).
            cursor: sqlite3.Cursor
                Cursor to the database (if database exists).
            asts_db: pd.DataFrame
                DataFrame with asteroids' information (if database exists).
        """
        self.asts = pd.DataFrame()
        self.asts_db = pd.DataFrame()
        self.n_jobs = n_jobs
        self._parallel = Parallel(n_jobs=n_jobs) if n_jobs > 1 else None

        try:
            alephs = OrbParams(parsing=parsing, filename=filename, update=update)
            self.filename = alephs.filename
            self.asts = alephs.asts
            self._asts = self.asts.to_records()
            self.ref_epochs = alephs.dicc_epochs_M
            self.ref_jdtdb = {t:self.ref_epochs[t].tdb.jd for t in self.ref_epochs}
            self.names = self.asts['name'].to_numpy()
            del alephs
        except:
            if not open_sql_database:
                raise ValueError("No data available.")
            self.asts = pd.DataFrame()
            self.ref_epochs = {}
            self.ref_jdtdb = {}
            self.names = np.array([])
        
        # Opening database
        if open_sql_database:
            if sqlfile is None: 
                sqlfile = PACKAGE_PATH / "aleph_states.db"
            if not Path(sqlfile).exists():
                raise FileNotFoundError(f"SQL database file '{sqlfile}' not found.")
            self.connection = sqlite3.connect(sqlfile)
            self.cursor = self.connection.cursor()
            self.asts_db = pd.read_sql('SELECT * FROM asts_info', self.connection, index_col='name')
        else:
            self.connection = None
            self.cursor = None
            self.asts_db = pd.DataFrame()
        
        if len(self.asts)==0 and len(self.asts_db)==0:
            raise ValueError("No data available.")

    def skyfield(
        self, 
        field_center: SkyCoord,
        field_radius: u.Quantity,
        epoch: Time|None = None,
        observer: EarthLocation|None = None,
        source: str = 'catalogue',
        method: str = 'mixed',
        ast_names: list[str]|np.ndarray|None = None,
        confidence_radius: u.Quantity|None = None,
        n_jobs: int = 1
        ):
        """
        Query for asteroides visible in a field of the sky.

        The asteroides positions can be computed using orbital parameters
        from public catalogues (if "source='catalogue'") or from an SQL 
        database (if "source='database'").

        There are several methods to compute the asteroids positions, but 
        in all cases, the asteroids positions are computed in two steps:
        1. Computing instantaneous asteroids positions at the requested epoch.
        2. Computing apparent asteroids positions only for the asteroids with
        instantaneus positions within a confidence radius around the requested
        field center.

        The asteroids positions can be computed using three different methods:
        '2b' which uses the keplerian equations from the 2-body problem
        (fast but inacurate), 'nb' which uses the numerical integration of the 
        n-body problem including the Sun and the 8 planets (very accurate but 
        extremely slow and computational expensive), or 'mixed' which precomputes
        instantaneous positions using the keplerian equations, and then uses the 
        n-body problem to refine the positions of the bodies within the
        'field_radius'.
        
        Parameters:
            field_center : astropy.coordinates.SkyCoord 
                Center coordinates of queried field.
            field_radius: astropy.units.Quantity 
                Radius of queried field.
            epoch: astropy.time.Time, default=None
                Queried epoch. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.
            n_jobs: int, default=1
                Number of parallel processes to calculate asteroids' 
                apparent positions.
            source: str, default='catalogue'
                Source of data to obtain the asteroids positions. Options:
                    * 'catalogue' : Uses the catalogue file.
                    * 'database' : Uses the SQL database.
            method: str, default='mixed'
                Method to calculate asteroids' positions. Options:
                * '2b': Computes ephemeris using 2-Body equations.
                        Very fast but inacurate a few days beyond the orbit's
                        reference epoch.
                * 'nb': Computes ephemeris using N-Bodies integration 
                        (Sun, planets and massless asteroids). Very slow, but
                        very accurate.
                * 'mixed': First it calculates ephemeris using 2-Body equations.
                           Then recalculates ephemeris using the N-Bodies
                           equations for those asteroids in and around the 
                           queried field.
            ast_names: list[str], default=None
                Asteroids' names from Query.asts. If provided, ephemerides 
                will be calculated only for those asteroids.
            confidence_radius: astropy.units.Quantity, default=None
                Specifies the radius around requested 'field_center' to select 
                the asteroids to calculate their apparent ephemeris. If None, 
                it uses field_radius+(3 degrees).
                           
        Returns:
            astropy.table : Table with computed ephemerides containing:
                column 'number': Asteroid's number.
                column 'name': Asteroid's name.
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
        # Validating params
        if source not in ['catalogue', 'database']: 
            raise ValueError("'source' must be one of 'catalogue' or 'database'")
        if method not in ['2b', 'nb', 'mixed']: 
            raise ValueError("'method' must be one of '2b', 'nb' or 'mixed'")
        
        # Selecting method
        if source == 'catalogue':
            if method == '2b':
                return self.skyfield_2b_cat(
                    field_center, field_radius, 
                    epoch=epoch, observer=observer, n_jobs=n_jobs, 
                    ast_names=ast_names, confidence_radius=confidence_radius
                    )
            elif method == 'nb':
                return self.skyfield_nb_cat(
                    field_center, field_radius, 
                    epoch=epoch, observer=observer, n_jobs=n_jobs,
                    **kwargs
                    )
            elif method == 'mixed':
                return self.skyfield_mixed_cat(
                    field_center, field_radius, 
                    epoch=epoch, observer=observer, n_jobs=n_jobs, 
                    **kwargs
                    )

        elif source == 'database':
            if method == '2b':
                return self.skyfield_2b_db(
                    field_center, field_radius, 
                    epoch=epoch, observer=observer, n_jobs=n_jobs, 
                    ast_names=ast_names, confidence_radius=confidence_radius
                    )
            elif method == 'nb':
                return self.skyfield_nb_db(
                    field_center, field_radius, 
                    epoch=epoch, observer=observer, n_jobs=n_jobs, 
                    ast_names=ast_names, confidence_radius=confidence_radius
                    )
            elif method == 'mixed':
                return self.skyfield_mixed_db(
                    field_center, field_radius, 
                    epoch=epoch, observer=observer, n_jobs=n_jobs, 
                    ast_names=ast_names, confidence_radius=confidence_radius
                    )
    
    
    def skyfield_2b_cat(
        self,
        field_center: SkyCoord,
        field_radius: u.Quantity,
        epoch: Time|None = None,
        observer: EarthLocation|None = None,
        ast_names: list[str]|np.ndarray|None = None,
        confidence_radius: u.Quantity|None = None,
        n_jobs: int|None = None,
        ):
        """
        Query for asteroides visible in a field of the sky using the 2-body equation 
        on the orbital parameters stored in 'astorb' or 'MPCORB' catalogues.
        
        Parameters:
            field_center : astropy.coordinates.SkyCoord 
                Center coordinates of queried field.
            field_radius: astropy.units.Quantity 
                Radius of queried field.
            epoch: astropy.time.Time, default=None
                Queried epoch. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.
            ast_names: list[str], default=None
                Asteroids' names from Query.asts. If provided, ephemerides 
                will be calculated only for those asteroids.
            confidence_radius: astropy.units.Quantity, default=None
                Specifies the radius around 'field_center' to select 
                asteroids to calculate their ephemeris. If None, it uses 
                field_radius+(3 degrees).
            n_jobs: int, default=None
                Number of parallel processes to calculate asteroids' 
                apparent positions.

        Returns:
            astropy.table : Table with computed ephemerides containing:
                column 'number': Asteroid's number.
                column 'name': Asteroid's name.
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
                column 'V': Apparent magnitude as seen by the observer in V band (calculated from the
                            asteroid's absolute magnitude).
        """
        # Getting all necessary parameters
        if epoch is None:
            epoch = Time.now()
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if ast_names is None:
            ast_idxs = np.arange(len(self.asts))
        else:
            ast_idxs = indexes_of_a_in_b(ast_names, self.names)
        if confidence_radius is None:
            confidence_radius = field_radius.to('deg').value + 3
        if n_jobs is None:
            n_jobs = self.n_jobs
            parallel = self._parallel
        elif isinstance(n_jobs, int):
            parallel = Parallel(n_jobs=n_jobs) if n_jobs > 1 else None
        else:
            raise ValueError("'n_jobs' must be int or None")
        
        # Computing observer's position
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd

        # Instantaneous states are calculated for all asteroids.
        asts = self._asts[ast_idxs]
        params = ['semi_major_axis', 'eccentricity', 'inclination', 
                  'longitude_ascending_node', 'argument_perihelion', 'mean_anomaly']
        states = []
        for ast in asts:
            tr = self.ref_jdtdb[ast['epoch']]
            states.append(pycutils.state_equatorial_heliocentric(tq, tr, *ast[params]))
        [xh ,yh, zh] = np.array(states)[:,:3].T
        
        # Instantaneous ephemerides are calculated for all asteroids.
        xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
        ra_geo  = np.arctan2( ye, xe ); dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) )
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra_geo, dec_geo)
        inconffield = dists<(confidence_radius*np.pi/180)
        # Ephemeris are calculated for those bodies inside confidence field
        asts = asts[inconffield]

        # Apparent coordinates are calculated
        if parallel is None:
            coords = []
            for ast in asts:
                tr = self.ref_jdtdb[ast['epoch']]
                args =  (tq, obsx, obsy, obsz, tr, *ast[params])
                coords.append(pycutils.apparent_coords_equatorial_heliocentric(*args))
        else:
            args = []
            for ast in asts:
                tr = self.ref_jdtdb[ast['epoch']]
                args.append((tq, obsx, obsy, obsz, tr, *ast[params]))
            coords = parallel(
                delayed(pycutils.apparent_coords_equatorial_heliocentric)(*arg)
                for arg in args
            )

        # Getting final coordinates
        coords = np.array(coords) # [ra,dec,x,y,z,vx,vy,vz,geodist] = coords.T # deg, au
        ra = coords[:,0]; dec = coords[:,1]
        dists = angular_separation(field_center.ra, field_center.dec, ra*u.deg, dec*u.deg)
        infield = dists.to('rad').value<field_radius.to('rad').value
        coords = coords[infield]; asts = asts[infield]
        nums = asts['number']; nms = asts['name']; hs = asts['H']; gs = asts['G']
        
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]   # Obs vector from ast
        ast_sun = [-x, -y, -z]      # Sun vector from ast
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

    def skyfield_nb_cat(
        self,
        field_center: SkyCoord,
        field_radius: u.Quantity,
        epoch: Time|None = None,
        observer: EarthLocation|None = None,
        ast_names: list[str]|np.ndarray|None = None,
        confidence_radius: u.Quantity|None = None,
        n_jobs: int|None = None
        ):
        """
        Query for asteroides visible in a field. Each asteroids 
        is initialized using its orbital parameters stored in the
        'astorb' or 'MPCORB' catalogues and then it is integrated 
        until requested epoch using an n-body simulation including 
        the Sun, the planets and the asteroid (as a massless particle).
        
        Parameters:
            field_center : astropy.coordinates.SkyCoord 
                Center coordinates of queried field.
            field_radius: astropy.units.Quantity 
                Radius of queried field.
            epoch: astropy.time.Time, default=None
                Queried epoch. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.
            ast_names: list[str], default=None
                Asteroids' names from Query.asts. If provided, ephemerides 
                will be calculated only for those asteroids.
            confidence_radius: astropy.units.Quantity, default=None
                Specifies the radius around 'field_center' to select 
                asteroids to calculate their ephemeris. If None, it uses 
                field_radius+(3 degrees).
            n_jobs: int, default=1
                Number of parallel processes to calculate asteroids' 
                apparent positions.
                    
        Returns:
            astropy.table : Table with computed ephemerides containing:
                column 'number': Asteroid's number.
                column 'name': Asteroid's name.
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
                             asteroid's absolute magnitude).
        """
        # Getting all necessary parameters
        if epoch is None:
            epoch = Time.now()
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if ast_names is None:
            ast_idxs = np.arange(len(self.asts))
        else:
            ast_idxs = indexes_of_a_in_b(ast_names, self.names)
        if confidence_radius is None:
            confidence_radius = field_radius.to('deg').value + 3
        if n_jobs is None:
            n_jobs = self.n_jobs
            parallel = self._parallel
        elif isinstance(n_jobs, int):
            parallel = Parallel(n_jobs=n_jobs) if n_jobs > 1 else None
        else:
            raise ValueError("'n_jobs' must be int or None")

        # Computing observer's position
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd

        # Getting Solar System's states in reference epochs
        set_ref_epochs = set(self._asts[ast_idxs]['epoch'])
        dicc_ss_states = {epoch:get_SS_states(self.ref_epochs[epoch], frame='helio') for epoch in set_ref_epochs}
        
        # Instantaneous states are calculated for all asteroids.
        asts = self._asts[ast_idxs]
        params = ['semi_major_axis', 'eccentricity', 'inclination', 
                  'longitude_ascending_node', 'argument_perihelion', 'mean_anomaly']
        states = []
        if parallel is None:
            for ast in asts:
                tr = self.ref_jdtdb[ast['epoch']]
                ast_state_tr = np.array([pycutils.state_equatorial_heliocentric(tr, tr, *ast[params])])
                ss_states = dicc_ss_states[ast['epoch']]
                delta_tdbs = tq - tr; dt = 0.1 * np.sign(delta_tdbs)
                states.append(get_integration(ast_state_tr, ss_states, dt, [delta_tdbs])[delta_tdbs][0])
        else:
            arguments = []
            for ast in asts:
                tr = self.ref_jdtdb[ast['epoch']]
                ast_state_tr = np.array([pycutils.state_equatorial_heliocentric(tr, tr, *ast[params])])
                ss_states = dicc_ss_states[ast['epoch']]
                delta_tdbs = tq-tr; dt = 0.1 * np.sign(delta_tdbs)
                arguments.append((ast_state_tr, ss_states, dt, [delta_tdbs]))
            states = parallel( delayed(get_integration)(*arg) for arg in arguments)
            states = [list(stt.values())[0][0] for stt in states]
        states = np.array(states)
        [xh ,yh, zh] = states[:,:3].T

        # Instantaneous ephemerides are calculated for all asteroids.
        xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
        ra_geo  = np.arctan2( ye, xe ); dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) )
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra_geo, dec_geo)
        inconffield = dists<(confidence_radius*np.pi/180)
        # Ephemeris are calculated for those bodies inside confidence field
        asts = asts[inconffield]
        states = states[inconffield]
        
        # Apparent states are calculated
        ss_states_q = get_SS_states(epoch, frame='helio')
        if parallel is None:
            coords = [get_apparent_coordinates(state, ss_states_q, obsx, obsy, obsz, 0) for state in states]
        else:
            arguments = [(state, ss_states_q, obsx, obsy, obsz, 0) for state in states]
            coords = parallel( delayed(get_apparent_coordinates)(*arg) for arg in arguments)
            
        # Getting final coordinates
        coords = np.array(coords); # [ra_geo, dec_geo, x ,y, z, vx,vy, vz, geodist] = coords.T # rad, au, au/d
        ra_geo = coords[:,0]; dec_geo = coords[:,1]
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra_geo, dec_geo)
        infield = dists<field_radius.to('rad').value
        coords = coords[infield]; coords[:,0] = coords[:,0]*180/np.pi; coords[:,1] = coords[:,1]*180/np.pi
        asts = asts[infield]; nums = asts['number']; nms = asts['name']; hs = asts['H']; gs = asts['G']
        
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]   # Obs vector from ast
        ast_sun = [-x, -y, -z]      # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)
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

    def skyfield_mixed_cat(
        self,
        field_center: SkyCoord,
        field_radius: u.Quantity,
        epoch: Time|None = None,
        observer: EarthLocation|None = None,
        ast_names: list[str]|np.ndarray|None = None,
        confidence_radius: u.Quantity|None = None,
        n_jobs: int|None = None
        ):
        """
        Query for asteroides visible in a field of the sky using a mixture of 
        the 2-body equation and n-body simulations, based on the orbital 
        parameters stored in 'astorb' or 'MPCORB' catalogues. Asteroid
        positions are initially calculated using the 2-body equation and then, 
        the apparent positions seen by the observer for those asteroids 
        within a 'confidence radius' (larger than the query's radius) are 
        refined using n-body simulations. In this way, the accuracy of the 
        ephemeris is increased while maintaining the computational efficiency.
    
        Parameters:
            field_center : astropy.coordinates.SkyCoord 
                Center coordinates of queried field.
            field_radius: astropy.units.Quantity 
                Radius of queried field.
            epoch: astropy.time.Time, default=None
                Queried epoch. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.
            ast_names: list[str], default=None
                Asteroids' names from Query.asts. If provided, ephemerides 
                will be calculated only for those asteroids.
            confidence_radius: astropy.units.Quantity, default=None
                Specifies the radius around 'field_center' to select 
                asteroids to calculate their ephemeris. If None, it uses 
                field_radius+(3 degrees).
            n_jobs: int, default=1
                Number of parallel processes to calculate asteroids' 
                apparent positions.
                    
        Returns:
            astropy.table : Table with computed ephemerides containing:
                column 'number': Asteroid's number.
                column 'name': Asteroid's name.
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
        # Getting all necessary parameters
        if epoch is None:
            epoch = Time.now()
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if ast_names is None:
            ast_idxs = np.arange(len(self.asts))
        else:
            ast_idxs = indexes_of_a_in_b(ast_names, self.names)
        if confidence_radius is None:
            confidence_radius = field_radius.to('deg').value + 3
        if n_jobs is None:
            n_jobs = self.n_jobs
            parallel = self._parallel
        elif isinstance(n_jobs, int):
            parallel = Parallel(n_jobs=n_jobs) if n_jobs > 1 else None
        else:
            raise ValueError("'n_jobs' must be int or None")
        
        # Computing observer's position
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd

        # Instantaneous states are calculated for all asteroids.
        asts = self._asts[ast_idxs]
        params = ['semi_major_axis', 'eccentricity', 'inclination', 
                  'longitude_ascending_node', 'argument_perihelion', 'mean_anomaly']
        states = []
        for ast in asts:
            tr = self.ref_jdtdb[ast['epoch']]
            states.append(pycutils.state_equatorial_heliocentric(tq, tr, *ast[params]))
        [xh ,yh, zh] = np.array(states)[:,:3].T

        # Instantaneous ephemerides are calculated for all asteroids.
        xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
        ra_geo  = np.arctan2( ye, xe ); dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) )
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra_geo, dec_geo)
        inconffield = dists<(confidence_radius*np.pi/180)
        # Ephemeris are calculated for those bodies inside confidence field
        asts = asts[inconffield]
        
        # Getting Solar System's states in pertinent reference epochs
        set_ref_epochs = set(self._asts[ast_idxs]['epoch'])
        dicc_ss_states = {epoch:get_SS_states(self.ref_epochs[epoch], frame='helio') for epoch in set_ref_epochs}
        
        # Apparent coordinates are calculated integrating the N-body problem
        if parallel is None:
            coords = []
            for ast in asts:
                tr = self.ref_jdtdb[ast['epoch']]; dt = tq-tr
                ss_states = dicc_ss_states[ast['epoch']]
                state_ast = pycutils.state_equatorial_heliocentric(tr, tr, *ast[params])
                coords.append(get_apparent_coordinates(state_ast, ss_states, obsx, obsy, obsz, dt))
        else:
            arguments = []
            for ast in asts:
                tr = self.ref_jdtdb[ast['epoch']]; dt = tq-tr
                ss_states = dicc_ss_states[ast['epoch']]
                state = pycutils.state_equatorial_heliocentric(tr, tr, *ast[params])
                arguments.append((state, ss_states, obsx, obsy, obsz, dt))
            coords = parallel( delayed(get_apparent_coordinates)(*arg) for arg in arguments)

        # Getting final coordinates
        coords = np.array(coords) #[ra,dec,x,y,z,vx,vy,vz,geodist] = coords.T # deg, au
        ra = coords[:,0]; dec = coords[:,1]
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra, dec)
        infield = dists<field_radius.to('rad').value
        coords = coords[infield]
        coords[:,0] = coords[:,0]*180/np.pi; coords[:,1] = coords[:,1]*180/np.pi
        asts = asts[infield]; nums = asts['number']; nms = asts['name']; hs = asts['H']; gs = asts['G']
        
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]   # Obs vector from ast
        ast_sun = [-x, -y, -z]      # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-gs)*phase_func1 + gs*phase_func2
        Vs = hs + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)
        
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

    def skyfield_2b_db(
        self,
        field_center: SkyCoord,
        field_radius: u.Quantity,
        epoch: Time|None = None,
        observer: EarthLocation|None = None,
        ast_names: list[str]|np.ndarray|None = None,
        confidence_radius: u.Quantity|None = None,
        n_jobs: int|None = None,
    ):
        """
        Query for asteroides visible in a field solving the 2-body equation 
        from the asteroids' orbital parameters saved in the SQL database.
        
        Parameters:
            field_center : astropy.coordinates.SkyCoord 
                Center coordinates of queried field.
            field_radius: astropy.units.Quantity 
                Radius of queried field.
            epoch: astropy.time.Time, default=None
                Queried epoch. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.
            ast_names: list[str], default=None
                Asteroids' names from Query.asts. If provided, ephemerides 
                will be calculated only for those asteroids.
            confidence_radius: astropy.units.Quantity, default=None
                Specifies the radius around 'field_center' to select 
                asteroids to calculate their ephemeris. If None, it uses 
                field_radius+(3 degrees).
            n_jobs: int, default=None
                Number of parallel processes to calculate asteroids' 
                apparent positions.

        Returns:
            astropy.table : Table with computed ephemerides containing:
                column 'number': Asteroid's number.
                column 'name': Asteroid's name.
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
        # Getting all necessary parameters
        if epoch is None:
            epoch = Time.now()
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if ast_names is None:
            allasts = True
        else:
            allasts = False
        if confidence_radius is None:
            confidence_radius = field_radius.to('deg').value + 3
        if n_jobs is None:
            n_jobs = self.n_jobs
            parallel = self._parallel
        elif isinstance(n_jobs, int):
            parallel = Parallel(n_jobs=n_jobs) if n_jobs > 1 else None
        else:
            raise ValueError("'n_jobs' must be int or None")
        
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
        
        # Getting instantaneous states using 2-body solution
        args0 = tq, tr
        states_ = np.array([pycutils.state_equatorial_heliocentric(*args0,*args1) for args1 in params])
        [xh ,yh, zh] = np.array(states_)[:,:3].T
        # Instantaneous ephemerides are calculated for all asteroids.
        xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
        ra  = np.arctan2( ye, xe ); dec = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) )
        
        # Ephemeris are calculated for those bodies inside confidence field
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra, dec)
        inconffield = dists<(confidence_radius*np.pi/180)
        tb=tb[inconffield]; params=params[inconffield]

        # Getting apparent coordinates
        args0 = tq, obsx, obsy, obsz, tr
        if parallel is None:
            coords = np.array([pycutils.apparent_coords_equatorial_heliocentric(*args0,*args1) for args1 in params])
        else:
            args = [(*args0,*args1) for args1 in params]
            pool = mp.Pool(processes=n_jobs)
            coords = np.array(pool.starmap(pycutils.apparent_coords_equatorial_heliocentric, args))
            pool.close()
            
        # Getting final coordinates
        ra = coords[:,0]; dec = coords[:,1]
        dists = angular_separation(field_center.ra, field_center.dec, ra*u.deg, dec*u.deg)
        infield = dists.to('rad').value<field_radius.to('rad').value
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


    def skyfield_nb_db(
        self,
        field_center: SkyCoord,
        field_radius: u.Quantity,
        epoch: Time|None = None,
        observer: EarthLocation|None = None,
        ast_names: list[str]|np.ndarray|None = None,
        confidence_radius: u.Quantity|None = None,
        n_jobs: int|None = None,
    ):
        """
        Query for asteroides visible in a field. Asteroids are initialized
        using their orbital states stored in the SQL database and then they 
        are integrated until the requested epoch using an n-body simulation 
        including the Sun, the planets and the asteroid (as a massless particle).
    
        Parameters:
            field_center : astropy.coordinates.SkyCoord 
                Center coordinates of queried field.
            field_radius: astropy.units.Quantity 
                Radius of queried field.
            epoch: astropy.time.Time, default=None
                Queried epoch. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.
            ast_names: list[str], default=None
                Asteroids' names from Query.asts. If provided, ephemerides 
                will be calculated only for those asteroids.
            confidence_radius: astropy.units.Quantity, default=None
                Specifies the radius around 'field_center' to select 
                asteroids to calculate their ephemeris. If None, it uses 
                field_radius+(3 degrees).
            n_jobs: int, default=None
                Number of parallel processes to calculate asteroids' 
                apparent positions.

        Returns:
            astropy.table : Table with computed ephemerides containing:
                column 'number': Asteroid's number.
                column 'name': Asteroid's name.
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
        # Getting all necessary parameters
        if epoch is None:
            epoch = Time.now()
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if ast_names is None:
            allasts = True
        else:
            allasts = False
        if confidence_radius is None:
            confidence_radius = field_radius.to('deg').value + 3
        if n_jobs is None:
            n_jobs = self.n_jobs
            parallel = self._parallel
        elif isinstance(n_jobs, int):
            parallel = Parallel(n_jobs=n_jobs) if n_jobs > 1 else None
        else:
            raise ValueError("'n_jobs' must be int or None")
        
        # Computing observer's position
        c = self.cursor; conn = self.connection
        [geopos,geovel] = observer.get_gcrs_posvel(epoch)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epoch)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
        tq = epoch.tdb.jd
        
        # Looking for table at nearest time
        c.execute("SELECT table_title, time, ABS(time-?) FROM state_times ORDER BY ABS(time-?)", (tq,tq))
        (table, t0, terr) = c.fetchone()
        if terr>0.25: 
            w = "Queried epoch is %.2f days away from closest epoch in the ephemeris database. " %terr
            w+= "User may want to increase 'db_query_radius' value."
            warnings.warn(w)
        
        # Looking for bodies near field
        ## Loading states from database
        if allasts: qw = 'SELECT * FROM %s' %table
        else: qw = 'SELECT * FROM %s WHERE name IN (%s)' %(table,str(ast_names)[1:-1])
        tb = pd.read_sql(qw, conn, index_col='name')
        states = tb.values; nms = tb.index.values
        x,y,z=states.T[:3]

        ## Bodies' equatorial coordinates from observer
        xe = x-obsx; ye = y-obsy; ze = z-obsz
        ra_geo  = np.arctan2( ye, xe ) * 180/np.pi
        dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) ) * 180/np.pi
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
        if parallel is None:
            ephs = []
            for i, state in enumerate(states):
                ephs.append(get_apparent_coordinates(state, *args_))
            ephs = np.array(ephs)
        else:
            pool = mp.Pool(processes=n_jobs)
            arguments = []
            for state in states:
                arguments.append((state, *args_))
            ephs = np.array(pool.starmap(get_apparent_coordinates, arguments))
            pool.close()
        
        # Getting final coordinates
        ephs = np.array(ephs)
        ra = ephs[:,0]; dec = ephs[:,1]
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra, dec)
        infield = dists<field_radius.to('rad').value
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
    
    def skyfield_mixed_db(
        self,
        field_center: SkyCoord,
        field_radius: u.Quantity,
        epoch: Time|None = None,
        observer: EarthLocation|None = None,
        ast_names: list[str]|np.ndarray|None = None,
        confidence_radius: u.Quantity|None = None,
        n_jobs: int|None = None,
):
        """
        Query for asteroides visible in a field of the sky using a mixture of 
        2-body equation and n-body simulations, based on the orbital parameters 
        and states saved in the SQL database. Asteroid positions are initially 
        calculated using the 2-body equation and then, the apparent positions 
        seen by the observer for those asteroids within a 'confidence radius' 
        (larger than the query's radius) are refined using n-body simulations. 
        In this way, the accuracy of the ephemeris is increased while 
        maintaining the computational efficiency.
    
        Parameters:
            field_center : astropy.coordinates.SkyCoord 
                Center coordinates of queried field.
            field_radius: astropy.units.Quantity 
                Radius of queried field.
            epoch: astropy.time.Time, default=None
                Queried epoch. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.
            ast_names: list[str], default=None
                Asteroids' names from Query.asts. If provided, ephemerides 
                will be calculated only for those asteroids.
            confidence_radius: astropy.units.Quantity, default=None
                Specifies the radius around 'field_center' to select 
                asteroids to calculate their ephemeris. If None, it uses 
                field_radius+(3 degrees).
            n_jobs: int, default=None
                Number of parallel processes to calculate asteroids' 
                apparent positions.

        Returns:
            astropy.table : Table with computed ephemerides containing:
                column 'number': Asteroid's number.
                column 'name': Asteroid's name.
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
       # Getting all necessary parameters
        if epoch is None:
            epoch = Time.now()
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        if ast_names is None:
            allasts = True
        else:
            allasts = False
        if confidence_radius is None:
            confidence_radius = field_radius.to('deg').value + 3
        if n_jobs is None:
            n_jobs = self.n_jobs
            parallel = self._parallel
        elif isinstance(n_jobs, int):
            parallel = Parallel(n_jobs=n_jobs) if n_jobs > 1 else None
        else:
            raise ValueError("'n_jobs' must be int or None")
        
        # Computing observer's position
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
        if parallel is None:
            ephs = np.array([get_apparent_coordinates(state, *args0) for state in states])
        else:
            pool = mp.Pool(processes=n_jobs)
            arguments = [(state, *args0) for state in states]
            ephs = np.array(pool.starmap(get_apparent_coordinates, arguments))
            pool.close()
        
        # Getting final coordinates
        ra = ephs[:,0]; dec = ephs[:,1]
        dists = angular_separation(field_center.ra.rad, field_center.dec.rad, ra, dec)
        infield = dists<field_radius.to('rad').value
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

    def astephs(
        self, 
        ast_name:str, 
        epochs: Time|None = None,
        observer: EarthLocation|None = None,
        source: str = 'catalogue',
        method: str = 'nb',
        ):
        """
        Returns ephemerides from specified asteroid at given epochs.

        The asteroides ephemerides can be computed using orbital parameters
        from public catalogues (if "source='catalogue'") or from an SQL 
        database (if "source='database'").

        The asteroids positions can be computed using two different methods:
        '2b' which uses the keplerian equations from the 2-body problem
        (fast but inacurate), and 'nb' which uses the numerical integration of the 
        n-body problem including the Sun and the 8 planets (very accurate much 
        slower and computational expensive).
        
        Parameters:
            ast_name: str 
                Asteroid's name.
            epochs: astropy.time.Time, default=None
                Queried epochs. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.
            source: str, default='catalogue'
                Source of data to obtain the asteroids positions. Options:
                    * 'catalogue' : Uses the catalogue file.
                    * 'database' : Uses the SQL database.
            method: str, default='mixed'
                Method to calculate asteroids' positions. Options:
                * '2b': Computes ephemeris using 2-Body equations.
                        Very fast but inacurate a few days beyond the orbit's
                        reference epoch.
                * 'nb': Computes ephemeris using N-Bodies integration 
                        (Sun, planets and massless asteroids). Very slow, but
                        very accurate.
                           
        Returns:
            astropy.table : Table with computed ephemerides containing:
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
        """
        # Validating params
        if source not in ['catalogue', 'database']: 
            raise ValueError("'source' must be one of 'catalogue' or 'database'")
        if method not in ['2b', 'nb']: 
            raise ValueError("'method' must be one of '2b' or 'nb'")
        
        # Selecting method
        if source == 'catalogue':
            if method == '2b':
                return self.astephs_2b_cat(
                    ast_name=ast_name, epochs=epochs, observer=observer,
                    )
            elif method == 'nb':
                return self.astephs_nb_cat(
                    ast_name=ast_name, epochs=epochs, observer=observer,
                    )

        elif source == 'database':
            if method == '2b':
                return self.astephs_2b_db(
                    ast_name=ast_name, epochs=epochs, observer=observer,
                    )
            elif method == 'nb':
                return self.astephs_nb_db(
                    ast_name=ast_name, epochs=epochs, observer=observer,
                    )

    def astephs_2b_cat(
        self, 
        ast_name:str, 
        epochs: Time|None = None,
        observer: EarthLocation|None = None,
    ):
        """
        Returns ephemerides from specified asteroid at given epochs.
        Ephemerides are calculated performing an n-body simulation that considers
        the sun, all 8 planets and the asteroid as a massless particle, starting from
        a time where all states are known and integrating all bodies position until
        requested time (correcting by light aberration).
    
        Parameters:
            ast: str
                Name of the asteroid.
            epochs: astropy.time.Time, default=None
                Epoch(s) to calculate ephemerides. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.

        Returns:
            astropy.table : Table with computed ephemerides containing:
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
                column 'V': Apparent magnitude as seen by the observer in V band (calculated from the
                            asteroid's absolute magnitude in Lowell's database).
        """
        # Getting all necessary parameters
        if epochs is None:
            epochs0 = Time([Time.now().jd], format='jd')
        else:
            if not isinstance(epochs, Time):
                raise TypeError('Expected astropy.time.Time, received %s' %(type(epochs)))
            if isinstance(epochs.jd, np.ndarray):
                epochs0 = epochs
            else: epochs0 = Time([epochs])
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        
        iepochs = np.argsort(epochs0.jd); epochs = epochs0[iepochs] # Epochs are sorted to speed up integration
        tqs  = epochs.tdb.jd
        tmin = min(tqs); tmax = min(tqs)
        
        # Looking for asteroid's parameters
        ast_params = self.asts.set_index('name').loc[ast_name]
        ast_ref_epoch = Time(ast_params['epoch'], format='iso', scale='tt')
        dt_tmin = abs(ast_ref_epoch.tdb.jd - tmin); dt_tmax = abs(ast_ref_epoch.tdb.jd - tmax)
        
        # If asteroids ref. epoch is closer to the latest queried epoch,
        # positions are calculated from there, backwards
        if dt_tmin > dt_tmax:
            tqs = tqs[::-1]; epochs = epochs[::-1]; iepochs = iepochs[::-1]
        
        # Observatory coordinates
        [geopos,geovel] = observer.get_gcrs_posvel(epochs)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epochs)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value
            
        # ASTEROID'S INTEGRATION
        params = ['semi_major_axis', 'eccentricity', 'inclination', 
                  'longitude_ascending_node', 'argument_perihelion', 'mean_anomaly']
        params_prev = ast_params[params].to_numpy().astype(float)
        params_prev_epoch = ast_ref_epoch.tdb.jd
        coords = []
        for i, tq in enumerate(tqs):
            # Asteroid's state at queried time tq
            state_ast = pycutils.state_equatorial_heliocentric(tq, params_prev_epoch, *params_prev)

            # Asteroid's orbital parameters at queried time tq
            # a, e, i, longnode, argperi, meananom
            ast_params_tq = pycutils.orbparams_from_state(*state_ast)

            # Asteroids ephemeris from observer at queried time tq
            # ra, dec, *state, lightdist
            ast_ephs_tq = pycutils.apparent_coords_equatorial_heliocentric(
                tq, obsx[i], obsy[i], obsz[i], tq, *ast_params_tq)
            coords.append(ast_ephs_tq)

            # Saving asteroid's params for next queried epoch
            params_prev = ast_params_tq
            params_prev_epoch = tq
        coords = np.array(coords)

        # Getting phase, elongation & V
        h = ast_params['H']; g = ast_params['G']
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]   # Obs vector from ast
        ast_sun = [-x, -y, -z]      # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)  # Phase angle (in rads)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-g)*phase_func1 + g*phase_func2
        Vs = h + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)

        geo_ast = [xe, ye, ze]            # Obs vector from ast
        geo_sun = [-obsx, -obsy, -obsz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)  # Elongation angle (in rads)
        to_deg = 180/np.pi

        # Joining results into final table
        t0 = QTable(epochs.jd.reshape(-1,1)*u.day, names=['epoch_jd'])
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        cunits=(u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        t1 = QTable(coords, names=cnames, units=cunits)
        t2 = QTable([alpha*to_deg*u.deg, elong*u.deg, Vs], names=('phase','elongation', 'V'))
        t = hstack([t0, t1, t2])[iepochs]

        return t

    def astephs_nb_cat(
        self, 
        ast_name:str, 
        epochs: Time|None = None,
        observer: EarthLocation|None = None,
    ):
        """
        Returns ephemerides from specified asteroid at given epochs.
        Ephemerides are calculated performing an n-body simulation that considers
        the sun, all 8 planets and the asteroid as a massless particle, starting from
        a time where all states are known and integrating all bodies position until
        requested time (correcting by light aberration).
    
        Parameters:
            ast: str
                Name of the asteroid.
            epochs: astropy.time.Time, default=None
                Epoch(s) to calculate ephemerides. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.

        Returns:
            astropy.table : Table with computed ephemerides containing:
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
                column 'V': Apparent magnitude as seen by the observer in V band (calculated from the
                            asteroid's absolute magnitude in Lowell's database).
        """
        # Getting all necessary parameters
        if epochs is None:
            epochs0 = Time([Time.now().jd], format='jd')
        else:
            if not isinstance(epochs, Time):
                raise TypeError('Expected astropy.time.Time, received %s' %(type(epochs)))
            if isinstance(epochs.jd, np.ndarray):
                epochs0 = epochs
            else: epochs0 = Time([epochs])
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        
        iepochs = np.argsort(epochs0.jd); epochs = epochs0[iepochs] # Epochs are sorted to speed up integration
        tqs  = epochs.tdb.jd
        tmin = min(tqs); tmax = min(tqs)
        
        # Looking for asteroid's parameters
        ast_params = self.asts.set_index('name').loc[ast_name]
        ast_ref_epoch = Time(ast_params['epoch'], format='iso', scale='tt')
        dt_tmin = abs(ast_ref_epoch.tdb.jd - tmin); dt_tmax = abs(ast_ref_epoch.tdb.jd - tmax)
        
        # If asteroids ref. epoch is closer to the latest queried epoch,
        # positions are calculated from there, backwards
        if dt_tmin > dt_tmax:
            t0 = tmax; dir_integ = -1
            tqs = tqs[::-1]; epochs = epochs[::-1]; iepochs = iepochs[::-1]
        else:
            t0 = tmin; dir_integ = 1
        
        # Observatory coordinates
        [geopos,geovel] = observer.get_gcrs_posvel(epochs)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epochs)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value

        # SYSTEM AT ASTEROID'S REFERENCE EPOCH
        ## Asteroid's state at its reference epoch
        params = ['semi_major_axis', 'eccentricity', 'inclination', 
                  'longitude_ascending_node', 'argument_perihelion', 'mean_anomaly']
        ast_state = pycutils.state_equatorial_heliocentric(
            ast_ref_epoch.tdb.jd, ast_ref_epoch.tdb.jd, *ast_params[params])

        ## Solar system states at ast. params. reference epoch
        ss_states = get_SS_states(ast_ref_epoch, frame='helio')

        ## Creating simulation
        sim = get_SS_sim(ss_states); sim.N_active = len(ss_states) # Solar system
        sim.add(x=ast_state[0], y=ast_state[1], z=ast_state[2], vx=ast_state[3], vy=ast_state[4], vz=ast_state[5]) # Asteroid
        sun = sim.particles[0]; ast = sim.particles[-1]; 

        ## Moving to closest queried epoch
        #delta_tdbs_0 = t0 - ast_ref_epoch.tdb.jd
        #sim.dt = 0.1*np.sign(delta_tdbs_0) # direction to first queried epoch
        #sim.integrate(delta_tdbs_0)
        sim.dt = 0.1*dir_integ             # direction from first queried epoch onward

        # ASTEROID'S INTEGRATION
        delta_tdbs = tqs - ast_ref_epoch.tdb.jd #t0
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
        h = ast_params['H']; g = ast_params['G']
        
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
    
    def astephs_2b_db(
        self, 
        ast_name:str, 
        epochs: Time|None = None,
        observer: EarthLocation|None = None,
    ):
        """
        Returns ephemerides from specified asteroid at given epochs.
        Ephemerides are calculated performing an n-body simulation that considers
        the sun, all 8 planets and the asteroid as a massless particle, starting from
        a time where all states are known and integrating all bodies position until
        requested time (correcting by light aberration).
    
        Parameters:
            ast: str
                Name of the asteroid.
            epochs: astropy.time.Time, default=None
                Epoch(s) to calculate ephemerides. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.

        Returns:
            astropy.table : Table with computed ephemerides containing:
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
                column 'V': Apparent magnitude as seen by the observer in V band (calculated from the
                            asteroid's absolute magnitude in Lowell's database).
        """
        # Getting all necessary parameters
        if epochs is None:
            epochs0 = Time([Time.now().jd], format='jd')
        else:
            if not isinstance(epochs, Time):
                raise TypeError('Expected astropy.time.Time, received %s' %(type(epochs)))
            if isinstance(epochs.jd, np.ndarray):
                epochs0 = epochs
            else: epochs0 = Time([epochs])
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        
        iepochs = np.argsort(epochs0.jd); epochs = epochs0[iepochs] # Epochs are sorted to speed up integration
        tqs  = epochs.tdb.jd
        tmin = min(tqs); tmax = min(tqs)
        
        # Looking for table at nearest time
        c = self.cursor
        c.execute("SELECT table_title, time, ABS(time-?) FROM param_times ORDER BY ABS(time-?)", (tmin,tmin))
        (table_min, t0_min, terr_min) = c.fetchone()
        c.execute("SELECT table_title, time, ABS(time-?) FROM param_times ORDER BY ABS(time-?)", (tmax,tmax))
        (table_max, t0_max, terr_max) = c.fetchone()
        if terr_min < terr_max:
            ast_ref_epoch = t0_min; table = table_min; dir_integ = 1
        else:
            ast_ref_epoch = t0_max; table = table_max; dir_integ = -1
            tqs = tqs[::-1]; epochs = epochs[::-1]; iepochs = iepochs[::-1]
        
        # Observatory coordinates
        [geopos,geovel] = observer.get_gcrs_posvel(epochs)
        obscoords = SkyCoord(geopos, frame='gcrs', obstime=epochs)
        obsx = obscoords.hcrs.cartesian.x.to('au').value
        obsy = obscoords.hcrs.cartesian.y.to('au').value
        obsz = obscoords.hcrs.cartesian.z.to('au').value

         # Asteroid's orbital parameters at nearest time
        c.execute('SELECT * FROM %s WHERE name="%s"' %(table,ast_name))
        ast_params = c.fetchone()[1:]
        
        # ASTEROID'S INTEGRATION
        params_prev = np.array(ast_params)
        params_prev_epoch = ast_ref_epoch
        coords = []
        for i, tq in enumerate(tqs):
            # Asteroid's state at queried time tq
            state_ast = pycutils.state_equatorial_heliocentric(tq, params_prev_epoch, *params_prev)

            # Asteroid's orbital parameters at queried time tq
            # a, e, i, longnode, argperi, meananom
            ast_params_tq = pycutils.orbparams_from_state(*state_ast)

            # Asteroids ephemeris from observer at queried time tq
            # ra, dec, *state, lightdist
            ast_ephs_tq = pycutils.apparent_coords_equatorial_heliocentric(
                tq, obsx[i], obsy[i], obsz[i], tq, *ast_params_tq)
            coords.append(ast_ephs_tq)

            # Saving asteroid's params for next queried epoch
            params_prev = ast_params_tq
            params_prev_epoch = tq
        coords = np.array(coords)

        # Getting phase, elongation & V
        c.execute('SELECT H,G FROM asts_info WHERE name="%s"' %(ast_name))
        h,g = c.fetchone()
        [x ,y, z] = coords[:,2:5].T; geodist = coords[:,8] # au
        xe, ye, ze = x-obsx, y-obsy, z-obsz
        ast_geo = [-xe, -ye, -ze]   # Obs vector from ast
        ast_sun = [-x, -y, -z]      # Sun vector from ast
        alpha = angle_between_vectors(ast_geo, ast_sun)  # Phase angle (in rads)
        phase_func1 = np.exp(-3.33*np.tan(alpha/2)**0.36)
        phase_func2 = np.exp(-1.87*np.tan(alpha/2)**1.22)
        phase_func = (1-g)*phase_func1 + g*phase_func2
        Vs = h + 5*np.log10(geodist*np.sqrt(x*x+y*y+z*z)) - 2.5*np.log10(phase_func)

        geo_ast = [xe, ye, ze]            # Obs vector from ast
        geo_sun = [-obsx, -obsy, -obsz]   # Sun vector from ast
        elong = angle_between_vectors(geo_ast, geo_sun)  # Elongation angle (in rads)
        to_deg = 180/np.pi

        # Joining results into final table
        t0 = QTable(epochs.jd.reshape(-1,1)*u.day, names=['epoch_jd'])
        cnames=('ra','dec','x','y','z','vx','vy','vz','lightdist')
        cunits=(u.deg, u.deg, u.au, u.au, u.au, u.au/u.day, u.au/u.day, u.au/u.day, u.au)
        t1 = QTable(coords, names=cnames, units=cunits)
        t2 = QTable([alpha*to_deg*u.deg, elong*u.deg, Vs], names=('phase','elongation', 'V'))
        t = hstack([t0, t1, t2])[iepochs]

        return t

    def astephs_nb_db(
        self, 
        ast_name:str, 
        epochs: Time|None = None,
        observer: EarthLocation|None = None,
    ):
        """
        Returns ephemerides from specified asteroid at given epochs.
        Ephemerides are calculated performing an n-body simulation that considers
        the sun, all 8 planets and the asteroid as a massless particle, starting from
        a time where all states are known and integrating all bodies position until
        requested time (correcting by light aberration).
    
        Parameters:
            ast: str
                Name of the asteroid.
            epochs: astropy.time.Time, default=None
                Epoch(s) to calculate ephemerides. If None, it uses the current time.
            observer: astropy.coordinates.EarthLocation, default=None
                Observer's position. If None, it uses the geocenter.

        Returns:
            astropy.table : Table with computed ephemerides containing:
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
                column 'V': Apparent magnitude as seen by the observer in V band (calculated from the
                            asteroid's absolute magnitude in Lowell's database).
        """
        # Getting all necessary parameters
        if epochs is None:
            epochs0 = Time([Time.now().jd], format='jd')
        else:
            if not isinstance(epochs, Time):
                raise TypeError('Expected astropy.time.Time, received %s' %(type(epochs)))
            if isinstance(epochs.jd, np.ndarray):
                epochs0 = epochs
            else: epochs0 = Time([epochs])
        if observer is None:
            observer = EarthLocation.from_geocentric(0,0,0, unit='m')
        
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
            t0 = t0_max; table = table_max; dir_integ = -1
            tqs = tqs[::-1]; epochs = epochs[::-1]; iepochs = iepochs[::-1]
        
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
