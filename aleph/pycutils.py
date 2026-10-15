from ._cutils import ffi, lib


def state_equatorial_heliocentric(
    epoch:float,
    paramsepoch:float,
    a:float,
    e:float,
    i:float,
    longnode:float,
    argperi:float,
    meananom:float,
) -> tuple[float, float, float, float, float, float]:
    """
    Given an epoch and orbital parameters (in heliocentric ecliptic coordinates)
    returns the state (in heliocentric equatorial coordinates).

    Parameters:
        epoch : float
            Epoch where the state is calculated (in days).
        paramsepoch: float
            Epoch of the orbital parameters (in days).
        a : float
            Semi-major axis (in au).
        e : float
            Eccentricity.
        i : float
            Inclination (in degrees).
        longnode: float
            Longitude of ascending node (in degrees).
        argperi: float
            Argument of perihelion (in degrees).
        meananom: float
            Mean anomaly (at paramsepoch, in degrees).
    
    Returns: 
        tuple[float]: 
            - (x, y, z): Position in au
            - (vx, vy, vz): Velocity in au/day
    """
    x = ffi.new("double *")
    y = ffi.new("double *")
    z = ffi.new("double *")
    vx = ffi.new("double *")
    vy = ffi.new("double *")
    vz = ffi.new("double *")

    lib.state_equatorial_heliocentric(
        epoch,
        paramsepoch,
        a,
        e,
        i,
        longnode,
        argperi,
        meananom,
        x,
        y,
        z,
        vx,
        vy,
        vz,
    )

    return (
        x[0],
        y[0],
        z[0],
        vx[0],
        vy[0],
        vz[0],
    )


def apparent_coords_equatorial_heliocentric(
    epoch:float,
    obsx:float,
    obsy:float,
    obsz:float,
    paramsepoch:float,
    a:float,
    e:float,
    i:float,
    longnode:float,
    argperi:float,
    meananom:float,
) -> tuple[float, float, float, float, float, float, float, float, float]:
    """
    Given an epoch, orbital parameters (in heliocentric ecliptic coordinates)
    and observer's coordinates (in heliocentric equatorial coordinates)
    returns the 'apparent' (as seen by the observer, namely, corrected by 
    light-travel time) coordinates, state and distance.

    Parameters:
        epoch : float
            Epoch where the state is calculated (in days).
        obsx : float
            Observer's 'x' coordinate (in au).
        obsy : float
            Observer's 'y' coordinate (in au).
        obsz : float
            Observer's 'z' coordinate (in au).
        paramsepoch: float
            Epoch of the orbital parameters (in days).
        a : float
            Semi-major axis (in au).
        e : float
            Eccentricity.
        i : float
            Inclination (in degrees).
        longnode: float
            Longitude of ascending node (in degrees).
        argperi: float
            Argument of perihelion (in degrees).
        meananom: float
            Mean anomaly (at pamsepoch, in degrees).
    
    Returns:
        tuple[float]:
            - (ra, dec): Equatorial coordinates of body as seen by the observer.
            - (x, y, z): Heliocentric state of body in au when observed light was emited.
            - (vx, vy, vz): Heliocentric velocity of body in au/day when observed light was emited.
            - dist : Apparent distance observer-body in au.
    """
    ra = ffi.new("double *")
    dec = ffi.new("double *")
    x = ffi.new("double *")
    y = ffi.new("double *")
    z = ffi.new("double *")
    vx = ffi.new("double *")
    vy = ffi.new("double *")
    vz = ffi.new("double *")
    gdist = ffi.new("double *")

    lib.apparent_coords_equatorial_heliocentric(
        epoch,
        obsx,
        obsy,
        obsz,
        paramsepoch,
        a,
        e,
        i,
        longnode,
        argperi,
        meananom,
        ra,
        dec,
        x,
        y,
        z,
        vx,
        vy,
        vz,
        gdist,
    )

    return (
        ra[0],
        dec[0],
        x[0],
        y[0],
        z[0],
        vx[0],
        vy[0],
        vz[0],
        gdist[0],
    )


def orbparams_from_state(
    x: float,
    y: float,
    z: float,
    vx: float,
    vy: float,
    vz: float,
) -> tuple[float, float, float, float, float, float]:
    """
    Given a body's state (in heliocentric equatorial coordinates), returns 
    orbital parameters (in heliocentric ecliptic coordinates).

    Parameters:
        x : float
            Body's 'x' coordinate (in au).
        y : float
            Body's 'y' coordinate (in au).
        z : float
            Body's 'z' coordinate (in au).
        vx : float
            Body's 'vx' velocity component (in au/day).
        vy : float
            Body's 'vy' velocity component (in au/day).
        vz : float
            Body's 'vz' velocity component (in au/day).
    
    Returns:
        tuple[float]:
            - a : Semi-major axis (in au).
            - e : Eccentricity.
            - i : Inclination (in degrees).
            - longnode: Longitude of ascending node (in degrees).
            - argperi: Argument of perihelion (in degrees).
            - meananom: Mean anomaly (at pams. epoch, in degrees).
    """
    a = ffi.new("double *")
    e = ffi.new("double *")
    i = ffi.new("double *")
    longnode = ffi.new("double *")
    argperi = ffi.new("double *")
    meananom = ffi.new("double *")

    lib.params_from_state(
        x,
        y,
        z,
        vx,
        vy,
        vz,
        a,
        e,
        i,
        longnode,
        argperi,
        meananom,
    )

    return (
        a[0],
        e[0],
        i[0],
        longnode[0],
        argperi[0],
        meananom[0],
    )

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
