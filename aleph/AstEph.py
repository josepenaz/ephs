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

def getint(s):
    try: return int(s)
    except: return np.nan

PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
IS_SQL_DB = os.path.exists(PACKAGE_PATH+'/aleph_states.db')

class AstEph:
    def __init__(self, service='Lowell', open_db=IS_SQL_DB, dbfile=PACKAGE_PATH+'/aleph_states.db'):
        """
        Initializing 'Query' enviroment. Orbital parameters are loaded
        and the sql database with asteroid's states is connected.

        Parameters:
                    
        """
        if (service=='MPC' and open_db==True):
            w = "You are using the 'MPC' survey AND opening the 'dbfile'='%s'" %dbfile
            w+= ". Only 'query_2b_cat()' and 'query_nb_cat()' will use 'MPCORB.DAT' data."
            warnings.warn(w)
        self.service = service
        alephs = DataBase.OrbParams(service)
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
            
            self.cursor.execute('SELECT name,number,H,G FROM asts_info')
            info = np.array(self.cursor.fetchall(), dtype=str)
            self.numbers_db = np.array([getint(num) for num in info[:,1]])
            self.names_db = np.array([nm.strip() for nm in info[:,0]])
            self.H_db = info[:,2].astype(float); self.G_db = info[:,3].astype(float)