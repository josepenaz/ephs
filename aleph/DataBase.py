import urllib, gzip
import os, time
import numpy as np
from . import pycutils
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric, get_body_barycentric_posvel, solar_system_ephemeris, EarthLocation, GCRS
import multiprocessing as mp
import sqlite3, warnings
from .aleph_utils import *
import pandas as pd

solar_system_ephemeris.set('de432s') 

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
def getfl(s):
    try: return float(s)
    except: return np.nan
def get_astdata(s, source='MPCORB'):
    if source=='MPCORB':
        pack_num = s[:7]; H = s[8:13]; G = s[14:19]; epoch = s[20:25]; M = s[26:35]; peri = s[37:46]; node = s[48:57]
        incl = s[59:68]; e = s[70:79]; n = s[80:91]; a = s[92:103]; U = s[105:106]; ref = s[107:116]; nobs = s[117:122]
        nops = s[123:126]; arc = s[127:136]; rms = s[137:141]; perts1 = s[142:145]; perts2 = s[146:149]; compnm = s[150:160]
        flags = s[161:165]; readdes = s[166:194]; lastobs = s[194:202]
        div = readdes.split(')'); name = div[-1].strip()
        if len(div)==2: num = div[0].strip()[1:]
        else: num = '-'
        dat = {'num':num, 'H':getfl(H), 'G':getfl(G), 'epoch':unpack_epoch(epoch,source), 'M_epoch':float(M),
               'Arg_Peri':float(peri), 'Node':float(node), 'incl':float(incl), 'e':float(e), 'n':float(n), 'a':float(a),
               'U':U, 'ref':ref, 'Nobs':getfl(nobs), 'Nops':float(nops), 'arc':arc, 'rms':getfl(rms), 'perts1':perts1,
               'perts2':perts2, 'compnm':compnm, 'flags':flags, 'name':name, 'lastobs':lastobs}
    if source=='astorb':
        num = s[:6]; name = s[7:25]; compnm = s[26:41]; H = s[42:47]; G = s[48:53]; color_B_V = s[54:58]
        iras_diam_km = s[59:64]; iras_taxonomy = s[65:69]; code6 = s[70:94]; arc_days = s[95:100]; nobs = s[101:105]
        epoch = s[106:114]; M = s[115:125]; peri = s[126:136]; node = s[137:147]; incl = s[147:157]; e = s[158:168]
        a = float(s[168:181]); compdate = s[182:190]; ceu_arcsec = s[191:198]; ceurate_arcsecday = s[199:207]
        ceudate = s[208:216]; peu_arcsec = s[217:224]; peudate = s[225:233]; peu_gratest_ceunext10yr = s[234:241]
        peudate_gratest_ceunext10yr = s[242:250]; peu_gratest_peunext10yr = s[251:258]; peudate_gratest_peunext10yr = s[259:267]
        dat = {'num':num, 'name':name, 'compnm':compnm, 'H':getfl(H), 'G':getfl(G), 'color_B_V':getfl(color_B_V),
               'iras_diam_km':getfl(iras_diam_km), 'iras_taxonomy':iras_taxonomy, 'code6':code6, 'arc_days':getfl(arc_days),
               'Nobs':getfl(nobs), 'epoch':unpack_epoch(epoch,source), 'M_epoch':float(M), 'Arg_Peri':float(peri),
               'Node':float(node), 'incl':float(incl), 'e':float(e), 'a':a, 'compdate':compdate, 'ceu_arcsec':ceu_arcsec,
               'ceurate_arcsecday':ceurate_arcsecday, 'ceudate':ceudate, 'peu_gratest_ceunext10yr':peu_gratest_ceunext10yr,
               'peudate_gratest_ceunext10yr':peudate_gratest_ceunext10yr, 'peu_gratest_peunext10yr':peu_gratest_peunext10yr,
               'peudate_gratest_peunext10yr':peudate_gratest_peunext10yr}
        P = 365.2568984 * a**1.5; dat['n'] = 360.0/P
    return dat
def unpack_epoch(packed_epoch, source='MPCORB'):
    if source=='MPCORB':
        S = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        y = int(str(S.find(packed_epoch[0]))+packed_epoch[1:3])
        m = S.find(packed_epoch[3]); d = S.find(packed_epoch[4])
        return '%04i-%02i-%02i' %(y, m, d)
    if source=='astorb':
        y = packed_epoch[:4]; m = packed_epoch[4:6]; d = packed_epoch[6:]
        return '%s-%s-%s' %(y, m, d)

def download_orbparams_catalogue(service, filename):
    """
    Donwloads the orbital parameter catalogue from the specified service.

    Parameters:
                service: Specifies origin of orbital parameters.
                         Can be 'Lowell' or 'MPC'.
                filename: Name of the saved file.
    """
    if service == 'MPC':
        #url = "https://www.minorplanetcenter.org/iau/MPCORB/MPCORB.DAT"
        url = 'https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz'
    elif service == 'Lowell':
        #url = 'ftp://ftp.lowell.edu/pub/elgb/astorb.dat.gz'
        ve = "Downloading 'astorb.dat' via FTP is not longer supported. "
        ve+= "Download it yourself in https://asteroid.lowell.edu/main/astorb/ and"
        ve+= "provide it via 'filename' parameter."
        raise ValueError(ve)
    else: raise ValueError("'service' must be 'Lowell' or 'MPC'")

    response = urllib.request.urlopen(url)
    fo = open(filename,'wb'); fo.write(gzip.decompress(response.read())); fo.close()


PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))

class OrbParams:
    def __init__ (self, service='MPC', update=False, **kwargs):
        """
        Reads an orbital parameter file for asteroids and computes ephemerides
        
        Parameters:
                    service: Sets origin of orbital parameters. Can be 'Lowell' or 'MPC'.
                             Default to 'Lowell'.
                    update: bool. If True, it downloads the orbital parameter's file
                            from the appropiate servive and saves it in the package's
                            directory. Default to False.
        Other Parameters:
                    filename: Str. Specifies file (and path) with orbital parameters to
                              read. This file must be in Lowell's or MPC's format.
        Attributes:
                    service: 'Lowell' or 'MPC'
                    asts: List of dictionaries containing all asteroids' information.
                    dicc_epochs_M: Dictionary with reference epochs of orbital parameters.
        """
        self.service = service
        if 'filename' in kwargs: self.filename = kwargs['filename']
        else:
            if service=='Lowell': self.filename = PACKAGE_PATH + '/astorb.dat'
            if service=='MPC': self.filename = PACKAGE_PATH + '/MPCORB.DAT'
        
        if update: self.update()
        
        self.asts = self.read(self.filename)
        
        # Getting dicc of epochs
        epochs_M = []
        for ast in self.asts: epochs_M.append(ast['epoch'])
        set_epoch_M = list(set(epochs_M)); epochs_M = Time(set_epoch_M, format='iso', scale='tt')
        self.dicc_epochs_M = {t:Time(t, format='iso', scale='tt') for t in set_epoch_M}
        
    def update(self):
        """
        It downloads the orbital parameter's file from the appropiate servive.
        """
        print('UPDATING',self.filename)
        download_orbparams_catalogue(self.service, self.filename)
        print("UPDATE FINISHED!")
                
    def read(self, filename):
        """
        Reads orbital parameter file and saves its data as a dicctionary in a list.
        """
        if not os.path.exists(filename):
            err = "Orbital parameter file not found. You need to run "
            err+= "'DataBase.download_orbparams_catalogue(service, filename)' or initialize OrbParams with "
            err+= "'update=True' or initialize OrbParams setting 'filename' with an orbital parameter "
            err+= "file in a different directory."
            raise ValueError(err)
        if self.service == 'MPC':
            asts = []
            is_header=True; is_numbered=False; is_unnumbered=False; is_1oposition=False
            i = 0; fi = open(filename,'r')
            for line in fi:
                i+=1
                if is_header: pass
                if '------' in line:
                    is_numbered=True; is_header=False
                    continue
                if is_numbered and len(line)==1:
                    is_numbered=False; is_unnumbered=True
                    continue
                if is_unnumbered and len(line)==1:
                    is_unnumbered=False; is_1oposition=True
                    continue
                if is_numbered or is_unnumbered or is_1oposition:
                    asts.append(get_astdata(line, source='MPCORB'))
            fi.close()
                
        elif self.service == 'Lowell':
            fi = open(filename,'r')
            asts = []
            i = 0
            for line in fi: asts.append(get_astdata(line, source='astorb'))
            fi.close()
                    
        else: raise ValueError("'service' must be 'Lowell' or 'MPC'")
        
        return asts
    
    def get_epoch_state(self, epoch_tdb_jd, ast_indexs=None, njobs=1):
        """
        Returns states of asteroids at specified epoch in heliocentric equatorial coordinates.
        One state consists of a list containing positions and velocities: [x,y,z,vx,vy,z].
        
        Parameters:
                    epoch_tdb_jd: float Epoch. Barycentric Dynamical Time in Julian Day format.
                    ast_indexs: Index list of asteroids in AlerceEphs.asts to use.
                                If not specified, states for all asteroids are calculated.
                    njobs: int Number of parallel processes to run calculations. Default 1.
        """
        if not isinstance(ast_indexs, (list,np.ndarray)): ast_indexs = np.arange(len(self.asts))
        t = epoch_tdb_jd # calculated only once
    
        if njobs==1:
            results = []
            for idx in ast_indexs:
                ast = self.asts[idx]
                r = pycutils.state_equatorial_heliocentric(t, t, ast['a'],ast['e'],ast['incl'],ast['Node'],ast['Arg_Peri'],ast['M_epoch'])
                results.append(r)
        else:
            pool = mp.Pool(processes=njobs)
            arguments = []
            for idx in ast_indexs:
                ast = self.asts[idx]
                arguments.append((t, t, ast['a'],ast['e'],ast['incl'],ast['Node'],ast['Arg_Peri'],ast['M_epoch']))
            results = pool.starmap(pycutils.state_equatorial_heliocentric, arguments)
            pool.close()
        return results
    
    def create_states_database(self, sqlfile=PACKAGE_PATH+'/aleph_states.db', ast_indexs=None, verbose=True): #, njobs=1
        """
        Computes states of all Lowell's asteroids in an interval of 101 days around ASTORB
        reference epoch with a time interval of 0.5 days and saves them in a SQL database.
        Since orbital parameters in MPCORB.DAT can be at different epochs, this utility only
        works with astorb.dat data.
        
        Parameters:
                    sqlfile: File name of SQL database. Default 'aleph_states.db'.
                    ast_indexs: Index list of asteroids in AlerceEphs.asts to use.
                                If not specified, states for all asteroids are calculated.
                    verbose: bool If True, print messages indicating progress. Default True.

        The SQL database contains the following tables:

         * 'asts_info'   : Contains asteroids' information. Columns: 'name', 'number', 'H', 'G'.
         * 'param_T0'    : Contains the orbital parameters of all asteroids as given in 
                           astorb.dat in columns 'name', 'a', 'e', 'i', 'longnode',
                           'argperi' and 'meananom' (heliocentric ecliptic coordinates).
         * 'param_TXXXX' : Same as 'param_T0', but with values in other epochs.
                           Only available after using function 'adding_params_into_db'.
         * 'param_times' : Contains columns 'table_title' (names with tables of states) and
                           'time' with the corresponding tdb time for those states.
         * 'state_T0'    : Contains the states of every asteroid in the astorb's epoch of 
                           osculation (common for all asteroids) in the columns 'name', 'x',
                           'y', 'z', 'vx', 'vy' and 'vz' (heliocentric equatorial coordinates).
         * 'state_TXXXX' : Same as 'state_T0', but with values in other epochs.
                           Only available after using function 'adding_states_into_db'.
         * 'state_times' : Contains columns 'table_title' (names with tables of states) and
                           'time' with the corresponding tdb time for those states.
        """
#                    njobs: int Number of parallel processes to run calculations. Default 1.
        if not isinstance(ast_indexs, (list,np.ndarray)): ast_indexs = np.arange(len(self.asts))
        epoch = Time(self.asts[0]['epoch'], format='iso', scale='tt')
        #t = np.round(epoch.tdb.jd,decimals=2) # calculated only once
        t0 = epoch.tdb.jd
        
        # Creating SQL DB
        #if verbose: print('* Creating SQL Database')
        pkey = [[self.asts[i]['name'].strip()] for i in ast_indexs]
        info = [[self.asts[i]['num'], self.asts[i]['H'], self.asts[i]['G']] for i in ast_indexs]
        
        ## Starting database
        os.system('rm %s' %sqlfile) # Database is constructed again
        conn = sqlite3.connect(sqlfile)
        ## Create param_times table
        conn.execute("CREATE TABLE param_times(table_title, time)")
        ## Create state_times table
        conn.execute("CREATE TABLE state_times(table_title, time)")
        ## Create and fill asts_info table
        conn.execute("CREATE TABLE asts_info(name TEXT, number TEXT, H REAL, G REAL, PRIMARY KEY(name))")
        conn.executemany("INSERT INTO asts_info (name, number, H, G) VALUES (?,?,?,?)", np.hstack((pkey, info)) )

        # Orbital parameters at reference epoch
        if verbose: print('* Saving orbital parameters at reference epoch')
        ti = time.time()
        params = []
        for idx in ast_indexs:
            ast = self.asts[idx]
            params.append([ast['name'].strip(),ast['a'],ast['e'],ast['incl'],ast['Node'],ast['Arg_Peri'],ast['M_epoch']])
        #params = np.array(params, dtype=str)
        ### Create state table
        table_name = "param_T0"
        conn.execute("CREATE TABLE param_T0(name TEXT,a REAL,e REAL,i REAL,longnode REAL,argperi REAL,meananom REAL, PRIMARY KEY(name))")
        ### Fill times table
        conn.execute("INSERT INTO param_times(table_title, time) VALUES (?, ?)", (table_name,t0))
        ### Fill state table
        conn.executemany("INSERT INTO param_T0(name,a,e,i,longnode,argperi,meananom) VALUES (?,?,?,?,?,?,?)", params)
        del params
        T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60
        if verbose: print('  Done in Time: %02i : %02i : %05.2f' %(h,m,s))
    
        # Positions at reference epoch
        if verbose: print('* Calculating positions at reference epoch')
        ti = time.time()
        states0 = self.get_epoch_state(t0, ast_indexs=ast_indexs)#, njobs=njobs)
        T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60
        print('  Done in Time: %02i : %02i : %05.2f' %(h,m,s))
        
        ## Reference epoch
        if verbose: print('* Saving states at reference epoch')
        ti = time.time()
        ### Create state table
        table_name = "state_T0"
        conn.execute("CREATE TABLE state_T0(name TEXT, x REAL,y REAL,z REAL,vx REAL,vy REAL,vz REAL, PRIMARY KEY(name))")
        ### Fill times table
        conn.execute("INSERT INTO state_times(table_title, time) VALUES (?, ?)", (table_name,t0))
        ### Fill state table
        conn.executemany("INSERT INTO state_T0(name,x,y,z,vx,vy,vz) VALUES (?,?,?,?,?,?,?)", np.hstack((pkey, states0)))
        T = time.time()-ti; h=int(T/3600); rest=T%3600; m = int(rest/60); s=rest%60
        if verbose: print('  Done in Time: %02i : %02i : %05.2f' %(h,m,s))

        # Closing db
        conn.commit()
        conn.close()

    def dbconn(self, sqlfile=PACKAGE_PATH+'/aleph_states.db'):
        """
        Returns an SQL connection (sqlite3.connect) to the database 
        (see 'sqlite3' documentation for all options). Remember to
        close the connection when finished or before using any 'aleph'
        function.
        """
        return sqlite3.connect(sqlfile)
        
    def adding_states_into_db (self, tdb_jds, sqlfile=PACKAGE_PATH+'/aleph_states.db', njobs=1, onebodysim=True, verbose=True):
        """
        Integrates states of all bodies in 'sqlfile' in the given 'tdb_jds' epochs
        starting from the states already stored in the database.

        Integration is made using `rebound.Simulation`, wich receives the Sun, all planets
        (and Pluto) and one or more asteroids (determined by 'onebodysim' parameter) starting
        from the states at a near time stored in the database ('sqlfile' parameter).
        
        Parameters:
                    tdb_jds: Array with the requested epochs in JD format and TDB scale.
                    sqlfile: File name of SQL database. Default 'aleph_states.db'.
                    njobs: int Number of parallel processes to run calculations. Default 1.
                    verbose: bool. If True, print messages indicating progress. Default True.
                    onebodysim: bool. If True, each `rebound.Simulation` receives only one asteroid.
                                If False, each `rebound.Simulation` receives several asteroids
                                (determined by 'njobs').
        
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
        pkey = np.reshape(df.index.values, (len(df),1))
        states0 = df.values
        ss_states = get_SS_states(Time(t_start, format='jd', scale='tdb'), frame='helio')

        # Integrating from starting states to requested times
        if dt>0: delta_tdbs = np.array(sorted(t_add)) - t_start
        else: delta_tdbs = np.array(sorted(t_add, reverse=True)) - t_start
        if njobs == 1:
            if verbose: print('* Integrating (can take several hours)')
            ti_ = time.time()
            if onebodysim:
                results = []
                for i in range(len(pkey)):
                    results.append(get_integration([states0[i]], ss_states, dt, delta_tdbs))
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
            l = len(pkey); intervalo=l*1./njobs
            argums = []
            if onebodysim:
                for i in range(l):
                    argums.append(([states0[i]], ss_states, dt, delta_tdbs))
            else:
                for i in range(njobs):
                    i_0 = int(i*intervalo); i_f = int((i+1)*intervalo)
                    if i==njobs-1: i_f = l
                    argums.append((states0[i_0:i_f], ss_states, dt, delta_tdbs))
            pool = mp.Pool(processes=njobs)
            ##results = [pool.apply(get_integration, args=arg) for arg in argums]
            ##res_raw = []
            ##for p in results:
            ##    r = p.get(); res_raw.append(r)
            results = pool.starmap(get_integration, argums)
            pool.close()
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
        ti_ = time.time()
        lt = len(ts_db)+1
        for j, it in enumerate(delta_tdbs):
            tdbs = t_start+it; table_name = "state_T%04i" %(j+lt)
            #if verbose: print('        ', j+lt, it, table_name, tdbs)
            data = np.hstack((pkey, dicc_integ[it]))
            ## Fill times table
            conn.execute("insert into state_times(table_title, time) values (?, ?)", (table_name,tdbs))
            ## Create state table    
            conn.execute("create table %s(name TEXT, x REAL,y REAL,z REAL,vx REAL,vy REAL,vz REAL, primary key(name))" %table_name)
            ## Fill state table
            conn.executemany("insert into "+table_name+"(name, x,y,z,vx,vy,vz) values (?,?,?,?,?,?,?)", data)
        T = time.time()-ti_; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
        if verbose: print('    Done in Time: %02i : %02i : %05.2f' %(h_,m_,s_))
        del dicc_integ; del data

        # Closing db
        conn.commit()
        conn.close()

    def adding_params_into_db (self, tdb_jds=None, sqlfile=PACKAGE_PATH+'/aleph_states.db', verbose=True):
        """
        Function that receives epochs with states already stored in the database and computes their
        corresponding orbital parameters (using the 2-body equations). Take into considerations that
        the states are computed in the heliocentric equatorial frame while the orbital parameters
        are in the heliocentric ecliptic frame (as those in astorb.dat or MPCORB.DAT).
        
        Parameters:
                    tdb_jds: Array with the requested epochs in JD format and TDB scale.
                             These epochs must have stored states. Default to None (meaning that
                             orbital parameters will be computed for states at all apochs).
                    sqlfile: File name of SQL database. Optional.
                    verbose: bool. If True, print messages indicating progress. Default True.
        """
        ttotal_ = time.time()
        # Opening SQL DB
        if verbose: print('* Opening SQL Database')
        conn = sqlite3.connect(sqlfile)
        c = conn.cursor()

        c.execute("SELECT time FROM state_times")
        ts_db = np.array(c.fetchall()).ravel()
        
        if tdb_jds is None: tdb_jds = ts_db
        
        # Checking validity of times
        t_intersec = set(tdb_jds)&set(ts_db)
        if len(t_intersec)==0:
            warnings.warn("None of the requested epochs have states calculated. Ending task.")
            return
        c.execute("SELECT time FROM param_times")
        ts_db = np.array(c.fetchall()).ravel()
        t_add = list(set(tdb_jds)-set(ts_db))
        if len(t_add)==0: 
            warnings.warn("All requested epochs have orbital parameters calculated. Ending task.")
            return

        # Iterating for each epoch
        for t in t_add:
            # Computing orbital parameters
            if verbose: print('    Calculating orbital parameters for epoch =',t)
            ti_ = time.time()
            c.execute("SELECT table_title FROM state_times WHERE time=?", (t,))
            (table,) = c.fetchone()
            
            df = pd.read_sql('SELECT * FROM %s' %table, conn, index_col='name')
            pkey = np.reshape(df.index.values, (len(df),1))
            states = df.values
            
            params = []
            for state in states: params.append( pycutils.orbparams_from_state(*state) )
            T = time.time()-ti_#; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
            if verbose: print('       Done in %05.2f s' %(T))
            
            # Saving orbital parameters
            if verbose: print('    Saving orbital parameters for epoch =',t)
            ti_ = time.time()
            ## Getting table name
            table_name = 'param_'+table.split('_')[1]
            ## Fill times table
            conn.execute("INSERT INTO param_times(table_title, time) VALUES (?, ?)", (table_name,t))
            ## Create state table
            data = np.hstack((pkey, params))
            conn.execute("create table %s(name TEXT,a REAL,e REAL,i REAL,longnode REAL,argperi REAL,meananom REAL, PRIMARY KEY(name))" %table_name)
            ## Fill state table
            conn.executemany("insert into "+table_name+"(name,a,e,i,longnode,argperi,meananom) values (?,?,?,?,?,?,?)", data)
            for state in states: params.append( pycutils.orbparams_from_state(*state) )
            T = time.time()-ti_#; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
            if verbose: print('       Done in %05.2f s' %(T))

        # Closing db
        conn.commit()
        conn.close()

        T = time.time()-ttotal_; h_=int(T/3600); rest=T%3600; m_ = int(rest/60); s_=rest%60
        if verbose: print('    Done in Time: %02i : %02i : %05.2f' %(h_,m_,s_))

    def get_available_db_times(self, sqlfile=PACKAGE_PATH+'/aleph_states.db'):
        """
        Returns an `astropy.time.Time` object with all available times stored in the database. 
        """
        conn = sqlite3.connect(sqlfile)
        c = conn.cursor()
        c.execute("SELECT time FROM state_times")
        ts = np.array(c.fetchall()).ravel()
        conn.close()
        return Time(ts, format='jd', scale='tdb')
