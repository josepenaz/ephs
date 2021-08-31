import requests as req
from ftplib import FTP
import os, time
import numpy as np
from scipy.optimize import root_scalar
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric, get_body_barycentric_posvel, solar_system_ephemeris, EarthLocation, GCRS
import multiprocessing as mp
import rebound, sqlite3
from .aleph_utils import *

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
        num = s[:7]; H = s[8:13]; G = s[14:19]; epoch = s[20:25]; M = s[26:35]; peri = s[37:46]; node = s[48:57]
        incl = s[59:68]; e = s[70:79]; n = s[80:91]; a = s[92:103]; U = s[105:106]; ref = s[107:116]; nobs = s[117:122]
        nops = s[123:126]; arc = s[127:136]; rms = s[137:141]; perts1 = s[142:145]; perts2 = s[146:149]; compnm = s[150:160]
        flags = s[161:165]; readdes = s[166:194]; lastobs = s[194:202]
        dat = {'pack_num':num, 'H':getfl(H), 'G':getfl(G), 'epoch':unpack_epoch(epoch,source), 'M_epoch':float(M),
               'Arg_Peri':float(peri), 'Node':float(node), 'incl':float(incl), 'e':float(e), 'n':float(n), 'a':float(a),
               'U':U, 'ref':ref, 'Nobs':getfl(nobs), 'Nops':float(nops), 'arc':arc, 'rms':getfl(rms), 'perts1':perts1,
               'perts2':perts2, 'compnm':compnm, 'flags':flags, 'readdes':readdes, 'lastobs':lastobs}
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


class AlerceEphs:
    def __init__ (self, service='Lowell', update=False):
        """
        Reads an orbital parameter file for asteroids and computes ephemerides
        
        Parameters:
                    service: Sets origin of orbital parameters. Can be 'Lowell' or 'MPC'
                    update: bool. If True, it downloads the orbital parameter's file
                            from the appropiate servive.
        Attributes:
                    service: 'Lowell' or 'MPC'
                    asts: List of dictionaries containing all asteroids' information.
                    epoch: Time for ephemerides. Default None.
                    observer: Observer's coordinate. Default None.
                    obsx: Observer's heliocentric equatorial x coordinates. Default None.
                    obsy: Observer's heliocentric equatorial y coordinates. Default None.
                    obsz: Observer's heliocentric equatorial z coordinates. Default None.
                    dicc_epochs_M: Dictionary with reference epochs of orbital parameters.
        """
        self.service = service
        
        if update: self.update()
        
        self.asts = self.read()
        self.epoch = None
        self.observer = None
        self.obsx = None; self.obsy = None; self.obsz = None
        
        # Getting dicc of epochs
        epochs_M = []
        for ast in self.asts: epochs_M.append(ast['epoch'])
        set_epoch_M = list(set(epochs_M)); epochs_M = Time(set_epoch_M, format='iso', scale='tt')
        self.dicc_epochs_M = {t.iso:t for t in epochs_M}
        
    def update(self):
        """
        Tt downloads the orbital parameter's file from the appropiate servive.
        """
        print('UPDATING')
        if self.service == 'MPC':
            #url = "https://www.minorplanetcenter.org/iau/MPCORB/MPCORB.DAT"
            url = "https://minorplanetcenter.net//iau/MPCORB/MPCORB.DAT"
            resp = req.get(url)
            fo = open('MPCORB.DAT','w'); fo.write(resp.text); fo.close()
            resp.close()
        elif self.service == 'Lowell':
            ftp = FTP('ftp.lowell.edu')     # connect to host, default port
            ftp.login('anonymous', 'jpena@das.uchile.cl')
            ftp.cwd('pub/elgb')
            flnm = 'astorb.dat.gz'
            ftp.retrbinary("RETR " + flnm, open(flnm, 'wb').write)
            try:ftp.quit()
            except: ftp.close()
            os.system('gunzip astorb.dat.gz')
        else: raise ValueError("'service' must be 'Lowell' or 'MPC'")
        print("UPDATING FINISHED!")
                
    def read(self):
        """
        Reads orbital parameter file and saves its data as a dicctionary in a list.
        """
        if self.service == 'MPC':
            asts = []
            is_header=True; is_numbered=False; is_unnumbered=False; is_1oposition=False
            i = 0; fi = open('MPCORB.DAT','r')
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
            fi = open('astorb.dat','r')
            asts = []
            i = 0
            for line in fi: asts.append(get_astdata(line, source='astorb'))
            fi.close()
                    
        else: raise ValueError("'service' must be 'Lowell' or 'MPC'")
        
        return asts
    
    def set_epoch(self, epoch): self.epoch = epoch
    def set_observer(self, observer):
        """
        Specifies observer's position
        Parameters:
                    observer: astropy.coordinates.SkyCoord containing observation's origin coordinate.
        """
        self.observer = observer
        self.obsx = observer.hcrs.cartesian.x.to('au').value
        self.obsy = observer.hcrs.cartesian.y.to('au').value
        self.obsz = observer.hcrs.cartesian.z.to('au').value