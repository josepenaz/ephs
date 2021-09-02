import ctypes
import os

# Find .so full suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"


PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))

cutils = ctypes.CDLL(PACKAGE_PATH+"/../cutils"+suffix)

cutils.state_equatorial_heliocentric.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                                 ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))

cutils.apparent_coords_equatorial_heliocentric.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                                           ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                                           ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))

cutils.params_from_state.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                     ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                     ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                     ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))

def state_equatorial_heliocentric(epoch, paramsepoch, a, e, i, longnode, argperi, meananom):
    """
    Given an epoch and orbital parameters (in heliocentric ecliptic coordinates)
    returns the state (in heliocentric equatorial coordinates).

    Parameters:
               epoch : float. Epoch where the state is calculated (in days).
               paramsepoch: float. Epoch of the orbital parameters (in days).
               a : float. Semi-major axis (in au).
               e : float. Eccentricity.
               i : float. Inclination (in degrees).
               longnode: float. Longitude of ascending node (in degrees).
               argperi: float. Argument of perihelion (in degrees).
               meananom: floaf. Mean anomaly (at pamsepoch, in degrees).
    
    Returns: tuple of floats: x, y, z, vx, vy, vz. State in au and au/day
    """
    jd_ = ctypes.c_double(epoch); tjdepoch_ = ctypes.c_double(paramsepoch)
    
    a_ = ctypes.c_double(a); e_ = ctypes.c_double(e); i_ = ctypes.c_double(i)
    longnode_ = ctypes.c_double(longnode)
    argperi_ = ctypes.c_double(argperi)
    meananom_ = ctypes.c_double(meananom)

    x = ctypes.c_double(0); y = ctypes.c_double(0); z = ctypes.c_double(0)
    vx = ctypes.c_double(0); vy = ctypes.c_double(0); vz = ctypes.c_double(0)
    
    args = (jd_, tjdepoch_, a_, e_, i_, longnode_, argperi_, meananom_,
            x, y, z, vx, vy, vz)

    cutils.state_equatorial_heliocentric(*args)
    return x.value, y.value, z.value, vx.value, vy.value, vz.value

def apparent_coords_equatorial_heliocentric(epoch, obsx, obsy, obsz, paramsepoch, a, e, i, longnode, argperi, meananom):
    """
    Given an epoch, orbital parameters (in heliocentric ecliptic coordinates)
    and observer's coordinates (in heliocentric equatorial coordinates)
    returns the 'apparent' (as seen by the observer, namely, corrected by 
    light-travel time) coordinates, state and distance.

    Parameters:
               epoch : float. Epoch where the state is calculated (in days).
               obsx : float. Observer's 'x' coordinate (in au).
               obsy : float. Observer's 'y' coordinate (in au).
               obsz : float. Observer's 'z' coordinate (in au).
               paramsepoch: float. Epoch of the orbital parameters (in days).
               a : float. Semi-major axis (in au).
               e : float. Eccentricity.
               i : float. Inclination (in degrees).
               longnode: float. Longitude of ascending node (in degrees).
               argperi: float. Argument of perihelion (in degrees).
               meananom: floaf. Mean anomaly (at pamsepoch, in degrees).
    
    Returns: tuple of floats: ra, dec, x, y, z, vx, vy, vz, dist.
             ra, dec : Equatorial coordinates of body as seen by the observer.
             x, y, z, vx, vy, vz : Heliocentric state of body in au and au/day
                                   when observed light was emited.
             dist : Apparent distance observer-body in au. 
                    Equivalent to sqrt(x*x+y*y+z*z). 
    """
    jd_ = ctypes.c_double(epoch); tjdepoch_ = ctypes.c_double(paramsepoch)
    
    xobs = ctypes.c_double(obsx)
    yobs = ctypes.c_double(obsy); zobs = ctypes.c_double(obsz)
    
    a_ = ctypes.c_double(a); e_ = ctypes.c_double(e); i_ = ctypes.c_double(i)
    longnode_ = ctypes.c_double(longnode)
    argperi_ = ctypes.c_double(argperi)
    meananom_ = ctypes.c_double(meananom)

    ra = ctypes.c_double(0); dec = ctypes.c_double(0)
    x = ctypes.c_double(0); y = ctypes.c_double(0); z = ctypes.c_double(0)
    vx = ctypes.c_double(0); vy = ctypes.c_double(0); vz = ctypes.c_double(0)
    gdist = ctypes.c_double(0)
    
    args = (jd_, xobs, yobs, zobs, tjdepoch_, a_, e_, i_, longnode_, argperi_, meananom_,
            ra, dec, x, y, z, vx, vy, vz, gdist)

    cutils.apparent_coords_equatorial_heliocentric(*args)
    return ra.value, dec.value, x.value, y.value, z.value, vx.value, vy.value, vz.value, gdist.value

def orbparams_from_state(x, y, z, vx, vy, vz):
    """
    Given a body's state (in heliocentric equatorial coordinates), returns 
    orbital parameters (in heliocentric ecliptic coordinates).

    Parameters:
               x : float. Body's 'x' coordinate (in au).
               y : float. Body's 'y' coordinate (in au).
               z : float. Body's 'z' coordinate (in au).
               vx : float. Body's 'vx' velocity component (in au/day).
               vy : float. Body's 'vy' velocity component (in au/day).
               vz : float. Body's 'vz' velocity component (in au/day).
    
    Returns: tuple of floats: a, e, i, longnode, z, argperi, meananom.
             a : Semi-major axis (in au).
             e : Eccentricity.
             i : Inclination (in degrees).
             longnode: Longitude of ascending node (in degrees).
             argperi: Argument of perihelion (in degrees).
             meananom: Mean anomaly (at pams. epoch, in degrees).
    """
    x_ = ctypes.c_double(x); y_ = ctypes.c_double(y); z_ = ctypes.c_double(z)
    vx_ = ctypes.c_double(vx); vy_ = ctypes.c_double(vy); vz_ = ctypes.c_double(vz)
    
    a = ctypes.c_double(0); e = ctypes.c_double(0); i = ctypes.c_double(0)
    longnode = ctypes.c_double(0); argperi = ctypes.c_double(0)
    meananom = ctypes.c_double(0)
    
    args = (x_, y_, z_, vx_, vy_, vz_, a, e, i, longnode, argperi, meananom)

    cutils.params_from_state(*args)
    return a.value, e.value, i.value, longnode.value, argperi.value, meananom.value

PI         = 3.14159265358979323846
TWOPI      = 6.283185307179586476925287
#C          = 173.14463348
G0         = 39.476926421373015
Msun       = 1.0000059769998062e+00
SSMASS     = 1.00134       # Total SS mass
DAY        = (1./365.25)	# Julian day, 86400 s
AXIAL_TILT = (23+26/60.+21.406/3600.)*PI/180. # tilt between equatorial and ecliptic coordinates ec=rot_x(eq,A_T); eq=rot_x(ec,-A_T)
ASEC2RAD   = 4.848136811095359935899141e-6 # Constant conversion from arcsec to rad. From NOVAS.
C_AUDAY    = 173.1446326846693
