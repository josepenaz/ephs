import astropy.constants as const
from astropy.coordinates import SkyCoord, get_body_barycentric, get_body_barycentric_posvel, solar_system_ephemeris, EarthLocation, GCRS
from astropy import units as u
import rebound
import numpy as np
from scipy.optimize import root_scalar

solar_system_ephemeris.set('de432s')

def zeroE(E,M,e): return M-E+(e*np.sin(E))
def get_heliocentric_xyz(epoch_jd, epoch_ref_jd, ast):
    """Returns equatorial coordinates"""
    # Moving parameters at obs epoch
    Mt = (ast['M_epoch'] + ast['n']*(epoch_jd-epoch_ref_jd))%360
    N = ast['Node']
    E = root_scalar(zeroE, args=(Mt*np.pi/180,ast['e']), bracket=(-2*np.pi, 2*np.pi)).root # In rads
    ecl = 23.4393
    # Getting coordinates in orbit's plane
    xplane = ast['a']*(np.cos(E) - ast['e']); yplane = ast['a']*np.sqrt(1-(ast['e']**2))*np.sin(E)
    # Getting heliocentric ecliptic coordinates
    v = np.arctan2(yplane,xplane) # True anomaly
    dist = np.sqrt(xplane**2+yplane**2)
    i = ast['incl']*np.pi/180; w = ast['Arg_Peri']*np.pi/180; N *= np.pi/180; ecl *= np.pi/180 # Parameters in radians
    xrot = np.cos(N) * np.cos(v+w) - np.sin(N) * np.sin(v+w) * np.cos(i)
    yrot = np.sin(N) * np.cos(v+w) + np.cos(N) * np.sin(v+w) * np.cos(i)
    zrot = np.sin(v+w) * np.sin(i)
    xh = dist * xrot; yh = dist * yrot; zh = dist * zrot
    # Velocities
    vx = -ast['a'] * (ast['n']*np.pi/180) * np.sin(E) / (1.0 - ast['e'] * np.cos(E)); vz = 0.0
    vy = np.sqrt(1.0 - ast['e']**2) * ast['a'] * (ast['n']*np.pi/180) * np.cos(E) / (1.0 - ast['e'] * np.cos(E))
    v = np.sqrt(vx**2+vy**2)
    ## Rotate by argument of perihelion in orbit plane
    cosw = np.cos(w); sinw = np.sin(w)
    xdp = vx * cosw - vy * sinw; ydp = vx * sinw + vy * cosw; zdp = vz;
    ## Rotate by inclination about x axis
    cosi = np.cos(i); sini = np.sin(i)
    xd = xdp; yd = ydp * cosi - zdp * sini; zd = ydp * sini + zdp * cosi;
    ## Rotate by longitude of node about z axis
    cosnode = np.cos(N); sinnode = np.sin(N)
    vxh = xd * cosnode - yd * sinnode; vyh = xd * sinnode + yd * cosnode; vzh = zd
    # Getting heliocentric equatorial coordinates
    xh, yh, zh = rot_x(xh, yh, zh, -ecl); vxh, vyh, vzh = rot_x(vxh, vyh, vzh, -ecl)
    return xh, yh, zh, vxh, vyh, vzh

au_to_d = (1*u.au/const.c.to('au/s')).to('d').value
def apparent_coords_equatorial_heliocentric(epoch_jd, obsx, obsy, obsz, epoch_ref_jd, ast):
    # Asteroid's coordinates
    xh, yh, zh, vxh, vyh, vzh = get_heliocentric_xyz(epoch_jd, epoch_ref_jd, ast)
    # Correcting asteroid's position by light-travel time
    gdist = np.sqrt((xh-obsx)**2+(yh-obsy)**2+(zh-obsz)**2)
    lighttime = gdist*au_to_d
    t3 = epoch_jd-lighttime; niter = 0
    while niter < 31:
        niter +=1; t2 = t3
        xh, yh, zh, vxh, vyh, vzh = get_heliocentric_xyz(t2, epoch_ref_jd, ast)
        gdist = np.sqrt((xh-obsx)**2+(yh-obsy)**2+(zh-obsz)**2)
        lighttime = gdist*au_to_d
        t3 = epoch_jd-lighttime
        #print('dt = ', (t3.jd - t2.jd))
        if abs(t3 - t2) < 1.0e-10: break
    # Spherical coords from Earth
    xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
    ra_geo  = np.arctan2( ye, xe ) * 180/np.pi; dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) ) * 180/np.pi
    return ra_geo,dec_geo,xh,yh,zh,vxh,vyh,vzh,gdist

def state_equatorial_heliocentric(epoch_jd, epoch_ref_jd, ast): return get_heliocentric_xyz(epoch_jd, epoch_ref_jd, ast)
    
################################################################
def get_apparent_coordinates(state_ast, states_ss, obsx, obsy, obsz, delta_tdb):
    # Creating simulation
    sim = get_SS_sim(states_ss); sim.N_active = sim.N # Solar system
    sim.add(x=state_ast[0], y=state_ast[1], z=state_ast[2], vx=state_ast[3], vy=state_ast[4], vz=state_ast[5]) # Asteroid
    sun = sim.particles[0]; ast = sim.particles[-1]; #sim.dt = -0.1
    # Integration until observation time
    sim.integrate(delta_tdb)
    # Correcting by light's travel time
    gdist = np.sqrt((ast.x-sun.x-obsx)**2+(ast.y-sun.y-obsy)**2+(ast.z-sun.z-obsz)**2)
    lighttime = gdist*au_to_d
    t3 = delta_tdb - lighttime
    niter=0
    while niter<31:
        niter += 1; t2=t3
        sim.integrate(t2)
        gdist = np.sqrt((ast.x-sun.x-obsx)**2+(ast.y-sun.y-obsy)**2+(ast.z-sun.z-obsz)**2)
        lighttime = gdist*au_to_d
        t3 = delta_tdb - lighttime
        if abs (t3 - t2) < 1.0e-10: break
    xh , yh, zh =  ast.x-sun.x , ast.y-sun.y , ast.z-sun.z
    vxh,vyh,vzh = ast.vx-sun.vx,ast.vy-sun.vy,ast.vz-sun.vz
    xe = xh-obsx; ye = yh-obsy; ze = zh-obsz
    ra_geo  = np.arctan2( ye, xe ); dec_geo = np.arctan2( ze, np.sqrt(xe*xe+ye*ye) )
    return ra_geo, dec_geo, xh,yh,zh, vxh,vyh,vzh, lighttime/au_to_d

# X-axis rotation
def rot_x(x,y,z, theta):
    ct = np.cos(theta); st = np.sin(theta);
    xr = x; yr = ct*y + st*z; zr =-st*y + ct*z;
    return xr, yr, zr
################################################################
def get_SS_states(epoch, planets=True, frame='bary'):
    """
    Computes states (positions and velocities) for the Solar System in a given epoch.
    Parameters:
                epoch: astropy.time.Time value.
                planets: Bool. If True, returns states of all planets (and Pluto).
                         If False, returns only Sun's state. Default True.
                frame: Which coordinate frame to use: 'bary' for barycentric and
                       'helio' for heliocentric. Default 'bary'.
    Returns: numpy.array with the solar system states
    """
    names = ['mercury', 'venus', 'earth-moon-barycenter', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto']
    states = []
    # Sun
    sun_pos_bary, sun_vel_bary = get_body_barycentric_posvel('sun', epoch)
    x = sun_pos_bary.x.to('au').value; y = sun_pos_bary.y.to('au').value; z = sun_pos_bary.z.to('au').value
    vx = sun_vel_bary.x.to('au/d').value; vy = sun_vel_bary.y.to('au/d').value; vz = sun_vel_bary.z.to('au/d').value
    if frame == 'bary': states.append([x,y,z,vy,vy,vz])
    elif frame == 'helio': states.append([0,0,0,0,0,0])
    else: raise ValueError("'frame' must be 'bary' or 'helio'")
    # Planets
    if planets:
        for i, name in enumerate(names):
            pos_bary, vel_bary = get_body_barycentric_posvel(name, epoch)
            if frame=='bary':
                pos = pos_bary; vel = vel_bary
            else:
                pos = pos_bary-sun_pos_bary; vel = vel_bary-sun_vel_bary
            x = pos.x.to('au').value; y = pos.y.to('au').value; z = pos.z.to('au').value
            vx = vel.x.to('au/d').value; vy = vel.y.to('au/d').value; vz = vel.z.to('au/d').value
            states.append([x,y,z,vx,vy,vz])
    return np.array(states)
def get_SS_sim(ss_states, planets=True):
    """
    Starts a Solar System simulation for integration.
    Parameters: 
                ss_states: numpy.array with solar system states
                           tu start simulation.
                planets: Bool. If False, considers only the Sun.
                         Default True.
    Returns: rebound.Simulation
    """
    masses = np.array([1, 1.6601141530543488e-07, 2.4478382877847715e-06, 3.040432648022642e-06, 3.2271560375549977e-07,
              0.0009547919152112404, 0.0002858856727222417, 4.36624373583127e-05, 5.151383772628674e-05, 7.361781606089469e-09])
    if not planets: masses = masses[:1]
    sim = rebound.Simulation()
    sim.units = ('d', 'AU', 'Msun')
    for i in range(len(masses)):
        st = ss_states[i]
        sim.add(m=masses[i],x=st[0], y=st[1], z=st[2], vx=st[3], vy=st[4], vz=st[5])
    return sim

#def get_integration(states, epoch0, dt, itimes):
def get_integration(states, ss_states, dt, itimes):
    """
    From initial states of asteroids at epoch0, returns dictionary with integrated
    positions of those asteroids at specified relative times from epoch0.
    Integration is performed using rebound
    
    Parameters:
                states: List of equatorial heliocentric states of all bodies to integrate from.
                epoch: astropy.time.Time Epoch to initialize Solar System simulation.
                dt: Maximum time-step for integration.
                relative_times: List with times (in days) after epoch0 to get integrated states.
    Returns: dicctionary with states, with the epoch (tdb.jd) as key.
    """
    # Starting Solar System
    #sim = start_SS(epoch0, ref='helio', planets=1); sim.move_to_hel()
    sim = get_SS_sim(ss_states)
    sim.N_active = sim.N
    # Adding asteroids
    for state in states:
        xh, yh, zh, vxh, vyh, vzh = state
        sim.add(x=xh, y=yh, z=zh, vx=vxh, vy=vyh, vz=vzh)
    # Integration for each time
    sim.dt = dt; sun = sim.particles[0]
    dicc_states = {}
    for t in itimes:
        sim.integrate(t); sim.move_to_hel()
        statesi = np.zeros((sim.N,6),dtype="float64")
        sim.serialize_particle_data(xyzvxvyvz=statesi)
        dicc_states[t] = statesi[10:]
    return dicc_states