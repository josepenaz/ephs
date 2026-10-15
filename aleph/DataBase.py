import urllib, gzip, shutil, os, time
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlretrieve
from collections.abc import Callable
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris
from joblib import Parallel, delayed
import sqlite3, warnings
from .aleph_utils import *
from . import pycutils

def unpack_mpc_epoch(packed_epoch):
    S = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    y = int(str(S.find(packed_epoch[0]))+packed_epoch[1:3])
    m = S.find(packed_epoch[3]); d = S.find(packed_epoch[4])
    return '%04i-%02i-%02i' %(y, m, d)

def unpack_lowell_epoch(packed_epoch):
    y = packed_epoch[:4]; m = packed_epoch[4:6]; d = packed_epoch[6:]
    return '%s-%s-%s' %(y, m, d)

def _astdata_from_mpc_line(s):
    readdes = s[166:194]
    div = readdes.split(')')
    name = div[-1].strip()
    if len(div)==2: num = div[0].strip()[1:]
    else: num = '-'
    return {
        'number':num, #unpack_number(s[:7])
        'H':_getfl(s[8:13]), 
        'G':_getfl(s[14:19]), 
        'epoch':unpack_mpc_epoch(s[20:25]), 
        'mean_anomaly':float(s[26:35]),
        'argument_perihelion':float(s[37:46]), 
        'longitude_ascending_node':float(s[48:57]), 
        'inclination':float(s[59:68]), 
        'eccentricity':float(s[70:79]), 
        'mean_motion':float(s[80:91]), 
        'semi_major_axis':float(s[92:103]),
        'U':s[105:106], 
        'ref':s[107:116], 
        'Nobs':_getfl(s[117:122]), 
        'Nops':float(s[123:126]), 
        'arc':s[127:136], 
        'rms':_getfl(s[137:141]), 
        'perts1':s[142:145],
        'perts2':s[146:149], 
        'compnm':s[150:160], 
        'flags':s[161:165], 
        'name':name, 
        'lastobs':s[194:202],
        }
           
def _astdata_from_astorb_line(s):
    return {
        'number':s[:6].strip() if s[:6].strip()!='' else '-', 
        'name':s[7:25].strip(), 
        'compnm':s[26:41], 
        'H':_getfl(s[42:47]), 
        'G':_getfl(s[48:53]), 
        'color_B_V':_getfl(s[54:58]),
        'iras_diam_km':_getfl(s[59:64]), 
        'iras_taxonomy':s[65:69], 
        'code6':s[70:94], 
        'arc_days':_getfl(s[95:100]),
        'Nobs':_getfl(s[101:105]), 
        'epoch':unpack_lowell_epoch(s[106:114]), 
        'mean_anomaly':float(s[115:125]), 
        'argument_perihelion':float(s[126:136]),
        'longitude_ascending_node':float(s[137:147]), 
        'inclination':float(s[147:157]), 
        'eccentricity':float(s[158:168]), 
        'semi_major_axis':float(s[168:181]), 
        'compdate':s[182:190], 
        'ceu_arcsec':s[191:198],
        'ceurate_arcsecday':s[199:207], 
        'ceudate':s[208:216], 
        'peu_arcsec':s[217:224], 
        'peudate':s[225:233], 
        'peu_gratest_ceunext10yr':s[234:241],
        'peudate_gratest_ceunext10yr':s[242:250], 
        'peu_gratest_peunext10yr':s[251:258],
        'peudate_gratest_peunext10yr':s[259:267],
        'mean_motion':360.0/(365.2568984 * float(s[168:181])**1.5)
        }

def parse_mpc(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"'{filename}' not found.")
    asts = []
    #is_numbered=True; is_unnumbered=False; is_1oposition=False
    in_header = True
    with open(filename,'r') as fi:
        for il, line in enumerate(fi):
            line = line[:-1]
            if set(line) == {'-'}:
                in_header = False
                continue
            if in_header: continue # Header and No content lines
            if len(line)==0: continue
            asts.append(_astdata_from_mpc_line(line))
    return pd.DataFrame(asts)

def parse_lowell(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"'{filename}' not found.")
    with open(filename,'r') as fi:
        asts = [_astdata_from_astorb_line(line) for line in fi]
    return pd.DataFrame(asts)

solar_system_ephemeris.set('de430') 

# https://stjarnhimlen.se/comp/ppcomp.html
# https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf

## FUNCTIONS TO READ ORBITAL PARAMETERS CATALOGUES
def unpack_number(packed_num):
    S = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    try: return int(packed_num)
    except:
        if '~' in packed_num:
            pw = np.array([62**3, 62**2, 62, 1])
            num = 620000
            for i, s in enumerate(packed_num[1:]): num += S.find(s)*pw[i]
        else:
            num = str(S.find(packed_num[0])) + packed_num[1:]
        return int(num)
def _getfl(s):
    try: return float(s)
    except: return np.nan

def _download_mpc(filename):
    """
    Download mpc orbital parameters catalogue from MPC.
    """
    url = 'https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz'
    filename = Path(filename)
    folder = filename.parent
    folder.mkdir(exist_ok=True)
    
    name_gz = Path(urlparse(url).path).name
    file_gz = folder / name_gz
    print(f"Downloading {file_gz}")
    response = urllib.request.urlopen(url)

    print(f"Uncompressing file into {filename}")
    with open(filename, 'wb') as fo:
        fo.write(gzip.decompress(response.read()))
    
def _download_lowell(filename):
    """
    Download lowell orbital parameters catalogue from Lowell.
    """
    url = "https://ftp.lowell.edu/pub/elgb/astorb.dat.gz"
    filename = Path(filename)
    folder = filename.parent
    folder.mkdir(exist_ok=True)

    name_gz = Path(urlparse(url).path).name
    file_gz = folder / name_gz

    # Downloads .gz file
    print(f"Downloading {file_gz}")
    urlretrieve(url, file_gz)

    # Uncompressing .gz file
    print(f"Uncompressing file into {filename}")
    with gzip.open(file_gz, "rb") as infile, open(filename, "wb") as outfile:
        shutil.copyfileobj(infile, outfile)

    print(f"Removing {file_gz}")
    file_gz.unlink()

PACKAGE_PATH = Path(os.path.abspath(os.path.dirname(__file__)))
Parser = Callable[[str | Path], pd.DataFrame]
PARSERS = {
    "Lowell": parse_lowell,
    "MPC": parse_mpc,
}
FILENAMES = {
    "Lowell": 'astorb.dat',
    "MPC": 'MPCORB.DAT',
}
DOWNLOADERS = {
    "Lowell": _download_lowell,
    "MPC": _download_mpc,
}
COMMON_COLUMNS = [
    "number",
    "name",
    "epoch",
    "mean_anomaly",
    "argument_perihelion",
    "longitude_ascending_node",
    "inclination",
    "eccentricity",
    #"mean_motion",  # Se calcula a partir de la masa del Sol
    "semi_major_axis",
    "H",
    "G",
]

class OrbParams:
    def __init__ (
        self,
        filename: str|Path|None = None,
        parsing: str|Parser = 'MPC',
        update: bool = False,
        n_jobs: int = 1,
        ):
        """
        Reads an orbital parameter file for asteroids and computes instantenous 
        osculating orbital elements.
        
        Parameters:
            filename : str, default=None
                Specifies file (and path) with orbital parameters to read.
            parsing : str, default='MPC'
                Sets parser to use. Can be 'MPC' or 'Lowell'.
            update : bool, default=False
                If True, it downloads the orbital parameter's file
                from the appropiate service and saves it in the package's
                directory.
        
        Attributes:
            asts: List of dictionaries containing all asteroids' information.
            dicc_epochs_M: Dictionary with reference epochs of orbital parameters.
            ref_epoch: The most common reference epoch to build SQL database.
        """
        available_parsers = ", ".join(PARSERS)
        if isinstance(parsing, str):
            if parsing not in PARSERS:
                raise ValueError(f"Unknown parser: {parsing}. Available: {available_parsers}")
            self.parser = PARSERS[parsing]
            self.parser_name = parsing
        elif callable(parsing):
            self.parser = parsing
            self.parser_name = parsing.__name__
        else:
            raise TypeError("'parsing' must be a string or a 'Callable[[str | Path], pd.DataFrame]' object.")

        downloader = None
        if isinstance(filename, (str, Path)):
            self.filename = Path(filename)
            if not self.filename.exists():
                raise FileNotFoundError(f"'{self.filename}' not found.")
        elif isinstance(parsing, str):
            self.filename = PACKAGE_PATH / FILENAMES[parsing]
            downloader = DOWNLOADERS[parsing]
        else:
            raise ValueError("Filename must be provided when 'parsing' is a 'Callable' object.")
        
        if not self.filename.exists() or update:
            if downloader is None:
                raise ValueError(f"If 'update' is True, 'parsing' must be a string in {available_parsers}")
            downloader(self.filename)
        
        self.asts = self.parser(self.filename)
        self._validate_common_schema()
        
        # Getting dicc of epochs
        epochs_M = self.asts['epoch'].values
        self.dicc_epochs_M = {t:Time(t, format='iso', scale='tt') for t in set(epochs_M)}
        # Getting orbital parameters reference epoch (the most frequent one)
        values, counts = np.unique(epochs_M, return_counts=True)
        self.ref_epoch = values[counts.argmax()]

        self.n_jobs = n_jobs
        self._parallel = Parallel(n_jobs=n_jobs) if n_jobs > 1 else None

    def _validate_common_schema(self):
        missing = set(COMMON_COLUMNS) - set(self.asts.columns)
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}")

    def __str__(self):
        msg = "OrbParams("
        msg+= f"parser : '{self.parser_name}', "
        msg+= f"number of bodies : {len(self.asts)}, "
        msg+= f"reference epoch : {self.ref_epoch}"
        msg+= ")"
        return msg

    def __repr__(self):
        return f"{self.__str__()}"
    
    def get_states_from_catalogue(
        self, 
        ast_names:list[float]|np.ndarray|None=None,
        ):
        """
        Returns states of asteroids at their reference epochs in 
        heliocentric equatorial coordinates.
        One state consists of a list containing positions and velocities: [x,y,z,vx,vy,z].
        
        Parameters:
            ast_names: list[str], default=None
                Asteroids' names from Query.asts. If provided, ephemerides 
                will be calculated only for those asteroids.
        """
        orb_param_cols = [
            'semi_major_axis', 'eccentricity', 'inclination', 'longitude_ascending_node', 
            'argument_perihelion', 'mean_anomaly'
            ]
        if ast_names is None:
            orb_params = self.asts[orb_param_cols].values
            orb_params_epochs = self.asts['epoch'].values
            orb_params_names = self.asts['name'].to_numpy()
        else:
            orb_params = self.asts.set_index('name').loc[ast_names]
            orb_params_epochs = orb_params['epoch'].values
            orb_params = orb_params[orb_param_cols].values
            orb_params_names = np.array(ast_names)
            
        if self._parallel is None:
            results = []
            for ast_idx, ast_params in enumerate(orb_params):
                t = self.dicc_epochs_M[orb_params_epochs[ast_idx]].tdb.jd
                r = pycutils.state_equatorial_heliocentric(t, t, *ast_params)
                results.append(r)
        else:
            arguments = []
            for ast_idx, ast_params in enumerate(orb_params):
                t = self.dicc_epochs_M[orb_params_epochs[ast_idx]].tdb.jd
                arguments.append((t, t, *ast_params))
            results = self._parallel(
                delayed(pycutils.state_equatorial_heliocentric)(*args)
                for args in arguments
            )
        results = np.hstack([orb_params_names.reshape(-1,1), results])
        
        return results

    def get_params_at_ref_epoch(
        self, 
        ast_names:list[str]|np.ndarray|None=None, 
        verbose:bool=False
        ):
        """
        Returns the equatorial heliocentric orbital parameters of all asteroids 
        from the specified catalogue in the common reference epoch 
        (attribute 'ref_epoch'). 
        
        Parameters:
            ast_names : list[str]|np.ndarray, default=None
                Names of asteroids to calculate the orbital parameters. 
                If None, the parameters of all asteroids in the catalogue 
                are calculated.
            verbose : bool, default=False
                If True, print messages indicating progress.
        """
        # Catalogue's orbital parameters
        orb_param_cols = [
            'semi_major_axis', 'eccentricity', 'inclination', 'longitude_ascending_node', 
            'argument_perihelion', 'mean_anomaly'
            ]
        if ast_names is None:
            orb_params = self.asts[orb_param_cols].to_numpy().copy()
            orb_params_epochs = self.asts['epoch'].to_numpy().copy()
            orb_params_names = self.asts['name'].to_numpy().copy()
        else:
            orb_params = self.asts.set_index('name').loc[ast_names].copy()
            orb_params_epochs = orb_params['epoch'].to_numpy().copy()
            orb_params_names = orb_params.index.to_numpy().copy()
            orb_params = orb_params[orb_param_cols].to_numpy().copy()
        ref_epoch = self.dicc_epochs_M[self.ref_epoch]
        ref_epoch = ref_epoch.tdb.jd

        if verbose: 
            task_str = '* Getting orbital parameters at reference epoch'
            print(task_str, end='\r')
        # Orbital parameters at catalogue reference epochs
        ti = time.time()

        # Moving orbital parameters to database reference epochs
        for epoch in set(orb_params_epochs):
            # Only bodies with other ref. epochs
            if epoch == self.ref_epoch: continue
            in_epoch = orb_params_epochs==epoch
            
            # Asteroids states at catalogue ref. epoch
            params0 = orb_params[in_epoch]
            t0 = Time(epoch, format='iso', scale='tt')
            states0 = np.array([pycutils.state_equatorial_heliocentric(
                t0.tdb.jd, t0.tdb.jd, *params) for params in params0])
            
            # SS states at catalogue ref. epoch
            ss_states = get_SS_states(t0, frame='helio')
            
            # Asteroids states and params at database ref. epoch
            delta_tdbs = ref_epoch - t0.tdb.jd; dt = 0.1 * np.sign(delta_tdbs)
            states_ = get_integration(states0, ss_states, dt, [delta_tdbs])[delta_tdbs]
            params_ = np.array([pycutils.orbparams_from_state(*state) for state in states_])
            
            # Savind states
            orb_params[in_epoch] = params_

        params = np.hstack([orb_params_names.reshape(-1,1), orb_params])
        if verbose:
            task_str = '* Getting orbital parameters at reference epoch'
            T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60 
            print(f'{task_str} (done in {h:02d}:{m:02d}:{s:05.2f})')
        return params
    
    def create_sql_database(
        self,
        sqlfile: str|Path|None = None,
        ast_names: list[str]|np.ndarray|None = None,
        verbose: bool = True
        ):
        """
        Creates an SQL database with the orbital parameters and states of 
        all asteroids from the specified catalogue in a common reference epoch.
        
        Since Lowell's catalogue has one common epoch for all asteroids, 
        this is the epoch used to initialize the SQL database. But MPC's
        catalogue has multiple epochs for different asteroids, so the most common
        epoch is used to initialize the SQL database (evolving the orbits to 
        this reference epoch using n-body simulations).

        This method creates the database with states and orbital parameters for
        the reference epoch only. Adding more epochs (states and orbital
        parameters) is possible using function 'adding_epochs_into_db'

        The SQL database contains the following tables:

        * 'asts_info'   : Contains asteroids' information that does not change over time. 
                          Columns: 'name', 'number', 'H', 'G'.
        * 'param_T0'    : Contains the orbital parameters of all asteroids 
                          as given in the source catalogue (Lowell's or MPC's).
                          Columns: 'name', 'a', 'e', 'i', 'longnode','argperi' 
                          and 'meananom' (heliocentric ecliptic coordinates).
        * 'param_TXXXX' : Same as 'param_T0', but with values in other epochs.
                          Only available after using function 'adding_epochs_into_db'.
        * 'param_times' : Summary table with all tables of orbital parameters
                          stored in the database.
                          Columns: 'table_title' (names of tables) and 'time' 
                          (the corresponding tdb time of those orbital 
                          parameters).
        * 'state_T0'    : Contains the states of all asteroids as given in the 
                          source catalogue (Lowell's or MPC's).
                          Columns: 'name', 'x', 'y', 'z', 'vx', 'vy' and 'vz' 
                          (heliocentric equatorial coordinates).
        * 'state_TXXXX' : Same as 'state_T0', but with values in other epochs.
                          Only available after using function 'adding_epochs_into_db'.
        * 'state_times' : Summary table with all tables of states stored in the 
                          database. 
                          Columns: 'table_title' (names of tables) and 'time' 
                          (the corresponding tdb time for those states).

        Parameters:
            sqlfile : str, default=None
                File name of SQL database. If None, it will create a new database
                with the name 'aleph_states.db' in the package directory.
            ast_names: list[str]|np.ndarray, default=None
                Names list of asteroids in AlerceEphs.asts to use.
                If not specified, states for all asteroids are calculated.
            verbose: bool, default=True
                If True, print messages indicating progress. 
        """
        # Setting arguments
        if sqlfile is None: 
            sqlfile = PACKAGE_PATH / 'aleph_states.db'
        else:
            sqlfile = Path(sqlfile)
        # Catalogue's orbital parameters
        orb_param_cols = [
            'semi_major_axis', 'eccentricity', 'inclination', 'longitude_ascending_node', 
            'argument_perihelion', 'mean_anomaly'
            ]
        if ast_names is None:
            orb_params = self.asts[orb_param_cols].to_numpy().copy()
            ast_names = self.asts['name'].to_numpy().copy()
            orb_params_info = self.asts[['number', 'H', 'G']].to_numpy().copy()
        else:
            orb_params = self.asts.set_index('name').loc[ast_names].copy()
            orb_params_info = orb_params[['number', 'H', 'G']].to_numpy().copy()
            orb_params = orb_params[orb_param_cols].to_numpy().copy()
        ref_epoch = self.dicc_epochs_M[self.ref_epoch]
        ref_epoch = ref_epoch.tdb.jd
        
        # Reference epoch
        epoch = self.dicc_epochs_M[self.ref_epoch]
        t0 = epoch.tdb.jd

        # Orbital parameters at reference epoch
        params = self.get_params_at_ref_epoch(ast_names=ast_names, verbose=verbose)

        # States at reference epoch
        if verbose: 
            print('* Calculating positions at reference epoch', end='\r')
        ti = time.time()
        states0 = np.array([pycutils.state_equatorial_heliocentric(
            t0, t0, *params_[1:]) for params_ in params])
        if verbose: 
            T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60
            verb_str = '* Calculating positions at reference epoch'
            print(f"{verb_str} (done in {h:02d} : {m:02d} : {s:05.2f})")
            #print('  Done in Time: %02i : %02i : %05.2f' %(h,m,s))
        
        # Creating SQL DB
        pkey = np.array(ast_names).reshape(-1,1) #[[self.asts[i]['name'].strip()] for i in ast_indexs]
        info = orb_params_info #[[self.asts[i]['number'], self.asts[i]['H'], self.asts[i]['G']] for i in ast_indexs]
        
        ## Starting database
        if verbose: 
            print('* Creating SQL Database')
            print('  * Creating informative tables', end='\r')
            ti = time.time()
        sqlfile.unlink(missing_ok=True) # Database is constructed again
        conn = sqlite3.connect(sqlfile)
        ## Create param_times table
        conn.execute("CREATE TABLE param_times(table_title, time)")
        ## Create state_times table
        conn.execute("CREATE TABLE state_times(table_title, time)")
        ## Create and fill asts_info table
        conn.execute("CREATE TABLE asts_info(name TEXT, number TEXT, H REAL, G REAL, PRIMARY KEY(name))")
        conn.executemany("INSERT INTO asts_info (name, number, H, G) VALUES (?,?,?,?)", np.hstack((pkey, info)) )
        if verbose:
            T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60
            verb_str = '  * Creating informative tables'
            print(f"{verb_str} (done in {h:02d} : {m:02d} : {s:05.2f})")

        ### Create state table
        if verbose: 
            print('  * Creating table with orbital parameters at reference epoch', end='\r')
            ti = time.time()
        table_name = "param_T0"
        conn.execute("CREATE TABLE param_T0(name TEXT,a REAL,e REAL,i REAL,longnode REAL,argperi REAL,meananom REAL, PRIMARY KEY(name))")
        ### Fill times table
        conn.execute("INSERT INTO param_times(table_title, time) VALUES (?, ?)", (table_name,t0))
        ### Fill state table
        conn.executemany("INSERT INTO param_T0(name,a,e,i,longnode,argperi,meananom) VALUES (?,?,?,?,?,?,?)", params)
        del params
        if verbose: 
            T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60
            verb_str = '  * Creating table with orbital parameters at reference epoch'
            print(f"{verb_str} (done in {h:02d} : {m:02d} : {s:05.2f})")
        
        ## Reference epoch
        if verbose: 
            print('  * Creating table with states at reference epoch', end='\r')
            ti = time.time()
        ### Create state table
        table_name = "state_T0"
        conn.execute("CREATE TABLE state_T0(name TEXT, x REAL,y REAL,z REAL,vx REAL,vy REAL,vz REAL, PRIMARY KEY(name))")
        ### Fill times table
        conn.execute("INSERT INTO state_times(table_title, time) VALUES (?, ?)", (table_name,t0))
        ### Fill state table
        conn.executemany("INSERT INTO state_T0(name,x,y,z,vx,vy,vz) VALUES (?,?,?,?,?,?,?)", np.hstack((pkey, states0)))
        if verbose: 
            T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60
            verb_str = '  * Creating table with states at reference epoch'
            print(f"{verb_str} (done in {h:02d} : {m:02d} : {s:05.2f})")

        # Closing db
        if verbose: 
            print('* Saving changes in SQL database', end='\r')
            ti = time.time()
        conn.commit()
        conn.close()
        if verbose: 
            T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60
            verb_str = '* Saving changes in SQL database'
            print(f"{verb_str} (done in {h:02d} : {m:02d} : {s:05.2f})")

    def dbconn(self, sqlfile=PACKAGE_PATH/'aleph_states.db'):
        """
        Returns an SQL connection (sqlite3.connect) to the database 
        (see 'sqlite3' documentation for all options). Remember to
        close the connection when finished or before using any 'aleph'
        function.
        """
        return sqlite3.connect(sqlfile)

    def adding_epochs_into_db (
        self, 
        epochs: float|list[float]|np.ndarray|Time, 
        sqlfile: str|Path|None = None,
        onebodysim: bool = True, 
        verbose: bool = True
        ):
        """
        Integrates states of all bodies in 'sqlfile' in the given 'tdb_jds' epochs
        starting from the states already stored in the database.

        Integration is made using `rebound.Simulation`, wich receives the Sun, all planets
        (and Pluto) and one or more asteroids (determined by 'onebodysim' parameter) starting
        from the states at a near time stored in the database ('sqlfile' parameter).
        
        Parameters:
            epochs : list[float]|np.ndarray
                Epoch or array of epochs to integrate orbits. If float or numpy array,
                the epochs must be in JD format and TDB scale. If astropy.Time object,
                the epochs are moved to TDB scale from the original scale.
            sqlfile: PathLike, default=None
                File name of SQL database. If not provided, it uses the default database
                created in the package directory.
            verbose: bool, default=True
                If True, print messages indicating progress.
            onebodysim: bool, default=True
                If True, each `rebound.Simulation` receives only one asteroid.
                If False, each `rebound.Simulation` receives several asteroids
                (determined by the attribute 'n_jobs').
        
        Setting 'onebodysim=True' occupies much more memory (required for all the planets in each
        simulation) but the simulations can run much faster. Setting 'onebodysim=False' saves much
        more memory, but if there are close encounters or other fenomena that could produce an 
        intensive integration (that affects not only to the asteroid that suffers it, but also all
        the other bodies in the simulation, that can be hundreds of thousands) the simulation can
        take an extremly long time to run.

        Take into account that you need to have around 300MB to load astorb.dat and around 115MB
        for one state of all asteroids (without taking into account the Sun and Planets' states)
        before you decide wether to set 'onebodysim' to True.
        """
        # Setting arguments
        if isinstance(epochs, Time):
            tdb_jds = epochs.tdb.jd
        else:
            tdb_jds = epochs
        if isinstance(tdb_jds, float):
            tdb_jds = [tdb_jds]
        if sqlfile is None: 
            sqlfile = PACKAGE_PATH / 'aleph_states.db'
        else:
            sqlfile = Path(sqlfile)
        # Opening SQL DB
        if verbose: print('* Opening SQL Database')
        conn = sqlite3.connect(sqlfile)
        c = conn.cursor()
        
        # Checking validity of times
        c.execute("SELECT time FROM state_times")
        ts_db = np.array(c.fetchall()).ravel()
        t_add = list(set(tdb_jds)-set(ts_db))
        if len(t_add)==0: 
            warnings.warn("All requested epochs are in the database. Ending task.")
            return
        
        # Getting nearest time in database
        if verbose: print('* Finding starting time in DB')
        min_t = min(t_add); max_t = max(t_add)
        c.execute("SELECT table_title, time, abs(time-?) FROM state_times ORDER BY ABS(time-?)", (min_t,min_t))
        (table_min, min_t_db, min_t_diff) = c.fetchone()
        c.execute("SELECT table_title, time, abs(time-?) FROM state_times ORDER BY ABS(time-?)", (max_t,max_t))
        (table_max, max_t_db, max_t_diff) = c.fetchone()
        if min_t_diff <= max_t_diff: 
            t_start=min_t_db; dt = .1
        else:
            t_start=max_t_db; dt = -.1
    
        # Getting states at nearest time
        if verbose: print('* Getting states at starting time in DB')
        c.execute("SELECT table_title FROM state_times WHERE time=?", (t_start,))
        (table,) = c.fetchone()
        df = pd.read_sql('SELECT * FROM %s' %table, conn, index_col='name')
        pkey = df.index.to_numpy().reshape(-1,1)
        states0 = df.to_numpy()
        ss_states = get_SS_states(Time(t_start, format='jd', scale='tdb'), frame='helio')
        
        # Integrating from starting states to requested times
        if dt>0: delta_tdbs = np.array(sorted(t_add)) - t_start
        else: delta_tdbs = np.array(sorted(t_add, reverse=True)) - t_start

        # Creating mapping dictionaries
        lt = len(ts_db)+1
        state_table_names = {it:"state_T%04i" %(j+lt) for j, it in enumerate(delta_tdbs)}
        param_table_names = {it:"param_T%04i" %(j+lt) for j, it in enumerate(delta_tdbs)}
        
        # Integrating from starting states to requested times
        if self._parallel is None:
            if verbose: print('* Integrating (can take several hours)')
            ti_ = time.time()
            if onebodysim:
                results = []
                for i in range(len(pkey)):
                    states0_i = np.array([states0[i]])
                    results.append(get_integration(states0_i, ss_states, dt, delta_tdbs))
                T = time.time()-ti_; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
                if verbose: print('    Done in Time: %02i : %02i : %05.2f' %(h_,m_,s_))
                if verbose: print('*   Joining results')
                dicc_integ = {}
                for dt in delta_tdbs:
                    sts = [dicc[dt] for dicc in results]
                    dicc_integ[dt] = np.vstack(sts)
                T = time.time()-ti_; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
                if verbose: print('        Done in Time: %02i : %02i : %05.2f' %(h_,m_,s_))
            else:
                dicc_integ = get_integration(states0, ss_states, dt, delta_tdbs)
                T = time.time()-ti_; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
                if verbose: print('    Done in Time: %02i : %02i : %05.2f' %(h_,m_,s_))
        else:
            if verbose: print('* Integrating in parallel (can take several hours)')
            ti_ = time.time()
            l = len(pkey); intervalo=l*1./self.n_jobs
            argums = []
            for i in range(self.n_jobs):
                i_0 = int(i*intervalo); i_f = int((i+1)*intervalo)
                if i==self.n_jobs-1: i_f = l
                argums.append((states0[i_0:i_f], ss_states, dt, delta_tdbs))
        
            if onebodysim:
                results = self._parallel(
                    delayed(get_multi_single_integrations)(*args)
                    for args in argums
                )
            else:
                results = self._parallel(
                    delayed(get_integration)(*args)
                    for args in argums
                )
            T = time.time()-ti_; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
            if verbose: print('    Done in Time: %02i : %02i : %05.2f' %(h_,m_,s_))
            if verbose: print('*   Joining parallel results')
            ti_ = time.time()
            dicc_integ = {}
            for dt in delta_tdbs:
                sts = [dicc[dt] for dicc in results]
                dicc_integ[dt] = np.vstack(sts)
            T = time.time()-ti_; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
            if verbose: print('        Done in Time: %02i : %02i : %05.2f' %(h_,m_,s_))

        # Saving Integration
        if verbose: print('* Saving states')
        ti = time.time()
        lt = len(ts_db)+1
        for j, it in enumerate(delta_tdbs):
            tdbs = t_start+it
            table_name = state_table_names[it]
            #if verbose: print('        ', j+lt, it, table_name, tdbs)
            data = np.hstack((pkey, dicc_integ[it]))
            ## Fill times table
            conn.execute("insert into state_times(table_title, time) values (?, ?)", (table_name,tdbs))
            ## Create state table    
            conn.execute("create table %s(name TEXT, x REAL,y REAL,z REAL,vx REAL,vy REAL,vz REAL, primary key(name))" %table_name)
            ## Fill state table
            conn.executemany("insert into "+table_name+"(name, x,y,z,vx,vy,vz) values (?,?,?,?,?,?,?)", data)
        
        # Calculating orbital parameters
        for it in dicc_integ:
            tdbs = t_start+it
            table_name = param_table_names[it]
            epoch_states = dicc_integ[it]
            epoch_params = [pycutils.orbparams_from_state(*state) for state in epoch_states]

            # Saving orbital parameters
            if verbose: print(f'    Saving orbital parameters for epoch = {tdbs}')
            ti_ = time.time()
            ## Fill times table
            conn.execute("INSERT INTO param_times(table_title, time) VALUES (?, ?)", (table_name,tdbs))
            ## Create state table
            data = np.hstack((pkey, epoch_params))
            conn.execute("create table %s(name TEXT,a REAL,e REAL,i REAL,longnode REAL,argperi REAL,meananom REAL, PRIMARY KEY(name))" %table_name)
            ## Fill state table
            conn.executemany("insert into "+table_name+"(name,a,e,i,longnode,argperi,meananom) values (?,?,?,?,?,?,?)", data)
            T = time.time()-ti_#; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
            if verbose: print('       Done in %05.2f s' %(T))
        
        T = time.time()-ti; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
        if verbose: print('    Done in Time: %02i : %02i : %05.2f' %(h_,m_,s_))
        del dicc_integ; del data

        # Closing db
        conn.commit()
        conn.close()

    def get_available_db_times(self, sqlfile=PACKAGE_PATH/'aleph_states.db'):
        """
        Returns an `astropy.time.Time` object with all available times stored in the database. 
        """
        conn = sqlite3.connect(sqlfile)
        c = conn.cursor()
        c.execute("SELECT time FROM state_times")
        ts = np.array(c.fetchall()).ravel()
        conn.close()
        return Time(ts, format='jd', scale='tdb')
