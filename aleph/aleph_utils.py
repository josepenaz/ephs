import astropy.constants as const
from astropy.coordinates import SkyCoord, get_body_barycentric, get_body_barycentric_posvel, solar_system_ephemeris, EarthLocation, GCRS
from astropy import units as u
from . import pycutils
import rebound
import numpy as np

solar_system_ephemeris.set('de432s')

au_to_d = (1*u.au/const.c.to('au/s')).to('d').value
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
    if ra_geo < 0: ra_geo += pycutils.TWOPI
    return ra_geo, dec_geo, xh,yh,zh, vxh,vyh,vzh, lighttime/au_to_d

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
    names = ['mercury', 'venus', 'earth-moon-barycenter', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']#, 'pluto']
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
              0.0009547919152112404, 0.0002858856727222417, 4.36624373583127e-05, 5.151383772628674e-05])#, 7.361781606089469e-09])
    if not planets: masses = masses[:1]
    sim = rebound.Simulation()
    sim.units = ('d', 'AU', 'Msun')
    for i in range(len(masses)):
        st = ss_states[i]
        sim.add(m=masses[i],x=st[0], y=st[1], z=st[2], vx=st[3], vy=st[4], vz=st[5])
    return sim

def get_integration(states, ss_states, dt, itimes):
    """
    From initial states of asteroids at epoch0, returns dictionary with integrated
    positions of those asteroids at specified relative times from epoch0.
    Integration is performed using rebound.
    
    Parameters:
                states: List of first equatorial heliocentric states of all bodies to integrate from.
                ss_states: List of first equatorial heliocentric states of all 8 planets.
                dt: Maximum time-step for integration.
                relative_times: List with times (in days) relative to first states' time to integrate.
    Returns: dicctionary with states, with the relative epoch (tdb.jd) as key.
    """
    # Starting Solar System
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
        dicc_states[t] = statesi[9:] #[10:] #[9:] without Pluto; [10:] with Pluto
    return dicc_states
################################################################
def angle_between_vectors(x1, x2):
    """
    Computes angle between vectors x1 and x2 using the dot and the cross product.
    
    Parameters:
                x1: array_like with _x_, _y_ and _z_ coordinates.
                x2: array_like with _x_, _y_ and _z_ coordinates.
                
    Returns:
             array with angle between x1 and x2 in radians.
    """
    #x1mod = np.sqrt((x1[0]*x1[0]) + (x1[1]*x1[1]) + (x1[2]*x1[2]))
    #x2mod = np.sqrt((x2[0]*x2[0]) + (x2[1]*x2[1]) + (x2[2]*x2[2]))
    x1_dot_x2 = (x1[0]*x2[0]) + (x1[1]*x2[1]) + (x1[2]*x2[2])
    x1_cross_x2_i = x1[1]*x2[2] - x1[2]*x2[1]
    x1_cross_x2_j = x1[2]*x2[0] - x1[0]*x2[2]
    x1_cross_x2_k = x1[0]*x2[1] - x1[1]*x2[0]
    x1_croos_x2_mod = np.sqrt(x1_cross_x2_i*x1_cross_x2_i+x1_cross_x2_j*x1_cross_x2_j+x1_cross_x2_k*x1_cross_x2_k)
    #return abs( np.arccos(x1x2/(x1mod*x2mod)) )
    return abs( np.arctan2( x1_croos_x2_mod , x1_dot_x2 ) )
    #return abs( np.arctan( x1_croos_x2_mod / x1_dot_x2) )
