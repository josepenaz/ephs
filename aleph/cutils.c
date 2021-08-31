/* gcc cutils.c -o cutils.o -fPIC -shared */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#include <time.h>*/
/*#include "novas.h"*/
/*#include "novascon.h"*/
#include "cutils.h"

/*#define TJDEPOCH 2452665.5*/

double gconst;
double machine_epsilon;

double determine_machine_epsilon()
{
  double u, den;
  
  u = 1.0;
  do {
    u /= 2.0;
    den = 1.0 + u;
  } while(den>1.0);

  return(10.0 * u);
}

/*-----------------  Functions from NOVAS   -----------------*/
/********bary2obs */
/* Receives a pos[3] (double) vector, a 'observer' vector pos_obs[3] (double), 
   another 'relative to observer' vector pos2[3] (double) and a lightime (double). 
   Function computes pos2 = pos - pos_obs and lightime = || pos2 || / C_AUDAY */
void bary2obs (double *pos, double *pos_obs,

               double *pos2, double *lighttime)
{
   short int j;

/* Translate vector to geocentric coordinates. */

   for (j = 0; j < 3; j++)
   {
      pos2[j] = pos[j] - pos_obs[j];
   }

/* Calculate length of vector in terms of light time. */

   *lighttime = sqrt (pos2[0] * pos2[0] + pos2[1] * pos2[1] +
       pos2[2] * pos2[2]) / C_AUDAY;

   return;
}

/********vector2radec */
/* Receives a position vector pos[3] (double) and computes its equatorial coordinates (ra, dec angles). }
   RETURNED
   VALUE:
      (short int)
         = 0 ... Everything OK.
         = 1 ... All vector components are zero; 'ra' and 'dec' are
                 indeterminate.
         = 2 ... Both pos[0] and pos[1] are zero, but pos[2] is nonzero;
                 'ra' is indeterminate.
*/
short int vector2radec (double *pos, double *ra, double *dec)
{
   double xyproj;

   xyproj = sqrt (pos[0] * pos[0] + pos[1] * pos[1]);
   if ((xyproj == 0.0) && (pos[2] == 0))
   {
      *ra = 0.0;
      *dec = 0.0;
      return 1;
   }
    else if (xyproj == 0.0)
   {
      *ra = 0.0;
      if (pos[2] < 0.0)
         *dec = -90.0;
       else
         *dec = 90.0;
      return 2;
   }
    else
   {
      *ra = atan2 (pos[1], pos[0]) / ASEC2RAD / 54000.0;
      *dec = atan2 (pos[2], xyproj) / ASEC2RAD / 3600.0;

      if (*ra < 0.0)
         *ra += 24.0;
   }
   return 0;
}
/*-----------------------------------------------------------*/

void rot_x(Vector v_in, Vector *v_out, double theta)
{
  
  double ct, st;

  ct = cos(theta);
  st = sin(theta);

  v_out->x = v_in.x;
  v_out->y = ct*v_in.y + st*v_in.z;
  v_out->z =-st*v_in.y + ct*v_in.z;
  
}

void rot_y(Vector v_in, Vector *v_out, double theta)
{
  
  double ct, st;

  ct = cos(theta);
  st = sin(theta);

  v_out->x = ct*v_in.x - st*v_in.z;
  v_out->y = v_in.y;
  v_out->z = st*v_in.x + ct*v_in.z;

  
}

void rot_z(Vector v_in, Vector *v_out, double theta)
{
  
  double ct, st;

  ct = cos(theta);
  st = sin(theta);

  v_out->x = ct*v_in.x + st*v_in.y;
  v_out->y =-st*v_in.x + ct*v_in.y;
  v_out->z = v_in.z;
  
}

/*-----------------  Functions for Orbital Elements   -----------------*/

/* Computes orbital elements from a State and 'gm' = 'Gravitational constant' * 'System Mass'*/
void keplerian(double gm, State state, 
	  double *a, double *e, double *i, double *longnode, double *argperi, double *meananom)
{
  double rxv_x, rxv_y, rxv_z, hs, h, parameter;
  double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
  double rcosu, rsinu, u, trueanom, eccanom;

  /* find direction of angular momentum vector */
  rxv_x = state.y * state.zd - state.z * state.yd;
  rxv_y = state.z * state.xd - state.x * state.zd;
  rxv_z = state.x * state.yd - state.y * state.xd;
  hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z;
  h = sqrt(hs);

  r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
  vs = state.xd * state.xd + state.yd * state.yd + state.zd * state.zd;
  /* v = sqrt(vs);  unnecessary */
  rdotv = state.x * state.xd + state.y * state.yd + state.z * state.zd;
  rdot = rdotv / r;
  parameter = hs / gm;

  *i = acos(rxv_z / h);

  if(rxv_x!=0.0 || rxv_y!=0.0) {
    *longnode = atan2(rxv_x, -rxv_y);
  } else {
    *longnode = 0.0;
  }

  *a = 1.0 / (2.0 / r - vs / gm);

  ecostrueanom = parameter / r - 1.0;
  esintrueanom = rdot * h / gm;
  *e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom);

  if(esintrueanom!=0.0 || ecostrueanom!=0.0) {
    trueanom = atan2(esintrueanom, ecostrueanom);
  } else {
    trueanom = 0.0;
  }

  cosnode = cos(*longnode);
  sinnode = sin(*longnode);

  /* u is the argument of latitude */
  rcosu = state.x * cosnode + state.y * sinnode;
  rsinu = (state.y * cosnode - state.x * sinnode)/cos(*i);

  if(rsinu!=0.0 || rcosu!=0.0) {
    u = atan2(rsinu, rcosu);
  } else {
    u = 0.0;
  }

  *argperi = u - trueanom;

  eccanom = 2.0 * atan(sqrt((1.0 - *e)/(1.0 + *e)) * tan(trueanom/2.0));
  *meananom = eccanom - *e * sin(eccanom);
  if (*meananom < 0) *meananom += TWOPI;

  return;
}

/* Computes a State from orbital elements and 'gm' = 'Gravitational constant' * 'System Mass'*/
void cartesian(double gm, 
	       double a, double e, double i, double longnode, double argperi, double meananom, 
	       State *state)
{
  double meanmotion, cosE, sinE, foo;
  double x, y, z, xd, yd, zd;
  double xp, yp, zp, xdp, ydp, zdp;
  double cosw, sinw, cosi, sini, cosnode, sinnode;
  double E0, E1, E2, den;
  int niter=0;

  /* first compute eccentric anomaly */
  E0 = meananom; 
  do {
    niter++;
    E1 = meananom + e * sin(E0);
    E2 = meananom + e * sin(E1);

    den = E2 - 2.0*E1 + E0;
    if(fabs(den) > machine_epsilon) {
      E0 = E0 - (E1-E0)*(E1-E0)/den;
    }
    else {
      E0 = E2;
      E2 = E1;
    }
    if(niter>3){
      E0 = 0.5*(E2+E1);
      niter=0;
    }
  } while(fabs(E0-E2) > machine_epsilon);

  cosE = cos(E0);
  sinE = sin(E0);

  /* compute unrotated positions and velocities */
  foo = sqrt(1.0 - e*e);
  meanmotion = sqrt(gm/(a*a*a));
  x = a * (cosE - e);
  y = foo * a * sinE;
  z = 0.0;
  xd = -a * meanmotion * sinE / (1.0 - e * cosE);
  yd = foo * a * meanmotion * cosE / (1.0 - e * cosE);
  zd = 0.0;

  /* rotate by argument of perihelion in orbit plane*/
  cosw = cos(argperi);
  sinw = sin(argperi);
  xp = x * cosw - y * sinw;
  yp = x * sinw + y * cosw;
  zp = z;
  xdp = xd * cosw - yd * sinw;
  ydp = xd * sinw + yd * cosw;
  zdp = zd;

  /* rotate by inclination about x axis */
  cosi = cos(i);
  sini = sin(i);
  x = xp;
  y = yp * cosi - zp * sini;
  z = yp * sini + zp * cosi;
  xd = xdp;
  yd = ydp * cosi - zdp * sini;
  zd = ydp * sini + zdp * cosi;

  /* rotate by longitude of node about z axis */
  cosnode = cos(longnode);
  sinnode = sin(longnode);
  state->x = x * cosnode - y * sinnode;
  state->y = x * sinnode + y * cosnode;
  state->z = z;
  state->xd = xd * cosnode - yd * sinnode;
  state->yd = xd * sinnode + yd * cosnode;
  state->zd = zd;

  return;
}
/*---------------------------------------------------------------------*/


/* state_equatorial_heliocentric: computes the position amd velocity of an orbit
**                                and returns it in equatorial cartesian coordinate
*/
void state_equatorial_heliocentric(double jd, double tjdepoch, double a, double e, double incl, double longnode, double argperi, double meananom, double *x, double *y, double *z, double *vx, double *vy, double *vz)
{
  double meananom_epoch;
  //double argperi_epoch;
  double mean_motion;
  double determine_machine_epsilon();
  void cartesian(double gm, double a, double e, double i, double longnode, double argperi, double meananom, State *state);
  State p;
  Vector vec_tmp, vec_tmp1;
  double tjd0;


  gconst = G0*DAY*DAY; /* heliocentric constant */

  incl           *= DEG2RAD; /*PI/180.*/
  longnode       *= DEG2RAD; /*PI/180.*/
  argperi        *= DEG2RAD; /*PI/180.*/
  meananom_epoch =  meananom*DEG2RAD; /*meananom*PI/180.*/
  
  mean_motion = sqrt(gconst/(a*a*a));

  tjd0 = jd;

  /* since the coordinates are given in ecliptic coordinates I should move later to equatorial */
  machine_epsilon = determine_machine_epsilon();
  cartesian(gconst, a, e, incl, longnode, argperi, meananom_epoch + (tjd0 - tjdepoch)*mean_motion, &p);

  /* before continuing I switch to equatorial coordinates as the vectors already computed */
  vec_tmp.x=p.x; vec_tmp.y=p.y; vec_tmp.z=p.z;
  rot_x(vec_tmp, &vec_tmp1, -AXIAL_TILT);
  p.x=vec_tmp1.x; p.y=vec_tmp1.y; p.z=vec_tmp1.z;
  
  vec_tmp.x=p.xd; vec_tmp.y=p.yd; vec_tmp.z=p.zd;
  rot_x(vec_tmp, &vec_tmp1, -AXIAL_TILT);
  p.xd=vec_tmp1.x; p.yd=vec_tmp1.y; p.zd=vec_tmp1.z;
  
  *x = p.x;
  *y = p.y;
  *z = p.z;
  *vx = p.xd;
  *vy = p.yd;
  *vz = p.zd;
}

/* Given state's values (in the equatorial frame) computes orbital parameters (in the ecliptic frame). */
void params_from_state(double x, double y, double z, double vx, double vy, double vz, double *a, double *e, double *i, double *longnode, double *argperi, double *meananom)
{
  State p;
  double a_, e_, i_, longnode_, argperi_, meananom_;
  void keplerian(double gm, State state, double *a, double *e, double *i, double *longnode, double *argperi, double *meananom);
  Vector pos_ec, vel_ec, pos_eq, vel_eq;

  /* Vectors for equatorial position and velocity */
  pos_eq.x = x; pos_eq.y = y; pos_eq.z = z;
  vel_eq.x =vx; vel_eq.y =vy; vel_eq.z =vz; 

  /* Equatorial vectors are rotated to get the ecliptic vectors */
  rot_x(pos_eq, &pos_ec, AXIAL_TILT);
  rot_x(vel_eq, &vel_ec, AXIAL_TILT);

  /* Constructions of ecliptic 'State' */
  p.x = pos_ec.x; p.y = pos_ec.y; p.z = pos_ec.z;
  p.xd= vel_ec.x; p.yd= vel_ec.y; p.zd= vel_ec.z; 

  gconst = G0*DAY*DAY; // G0*Msun*SSMASS*DAY*DAY;
  
  /* Getting orbital parameters */
  keplerian(gconst, p, &a_, &e_, &i_, &longnode_, &argperi_, &meananom_);

  *a = a_;
  *e = e_;
  *i = i_*RAD2DEG;
  *longnode = longnode_*RAD2DEG;
  *argperi = argperi_*RAD2DEG;
  *meananom = meananom_*RAD2DEG;
}

/* apparent_coords_equatorial_heliocentric: Given orbital parameters, observer's position and epoch, returns body's 
                                            apparent coordinates, its cartesian coordinates at light emission time
                                            and distance traveled by light.
*/
void apparent_coords_equatorial_heliocentric(double jd, double pobsx,  double pobsy,  double pobsz, double tjdepoch, double a, double e, double incl, double longnode, double argperi, double meananom, double *ra, double *dec, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *gdist)
{
  double meananom_epoch;
  double mean_motion;
  double determine_machine_epsilon();
  void cartesian(double gm, double a, double e, double i, double longnode, double argperi, double meananom, State *state);
  State p;
  Vector pobs, vec_tmp, vec_tmp1;
  double pob[3], pos1[3], pos2[3];
  double tjd0, lighttime, t2, t3;

  int niter=0;
  double 

  gconst = G0*DAY*DAY; // G0*Msun*SSMASS*DAY*DAY;

  incl           *= DEG2RAD; /*PI/180.*/
  longnode       *= DEG2RAD; /*PI/180.*/
  argperi        *= DEG2RAD; /*PI/180.*/
  meananom_epoch =  meananom*DEG2RAD; /*meananom*PI/180.*/
  
  mean_motion = sqrt(gconst/(a*a*a));
  tjd0 = jd;

  pobs.x = pobsx;
  pobs.y = pobsy;
  pobs.z = pobsz;

  pob[0] = pobs.x;
  pob[1] = pobs.y;
  pob[2] = pobs.z;

  /* since the coordinates are given in ecliptic coordinates I should move later to equatorial */
  machine_epsilon = determine_machine_epsilon();
  cartesian(gconst, a, e, incl, longnode, argperi, meananom_epoch + (tjd0 - tjdepoch)*mean_motion, &p);
  /* before continuing I switch to equatorial coordinates as the vectors already computed */
  // x,y,z
  vec_tmp.x=p.x; vec_tmp.y=p.y; vec_tmp.z=p.z;
  rot_x(vec_tmp, &vec_tmp1, -AXIAL_TILT);
  p.x=vec_tmp1.x; p.y=vec_tmp1.y; p.z=vec_tmp1.z;
  // xd,yd,zd
  vec_tmp.x=p.xd; vec_tmp.y=p.yd; vec_tmp.z=p.zd;
  rot_x(vec_tmp, &vec_tmp1, -AXIAL_TILT);
  p.xd=vec_tmp1.x; p.yd=vec_tmp1.y; p.zd=vec_tmp1.z;

  pos1[0] = p.x;
  pos1[1] = p.y;
  pos1[2] = p.z;



  bary2obs (pos1,pob, pos2,&lighttime); /*bary_to_geo (pos1,pob, pos2,&lighttime);*/
  t3 = tjd0 - lighttime;

  do{
	niter++;
	t2 = t3;

	cartesian(gconst, a, e, incl, longnode, argperi, meananom_epoch + (t2-tjdepoch)*mean_motion, &p);
	/* before continuing I switch to equatorial coordinates as the vectors already computed */
	vec_tmp.x=p.x; vec_tmp.y=p.y; vec_tmp.z=p.z;
	rot_x(vec_tmp, &vec_tmp1, -AXIAL_TILT);
	p.x=vec_tmp1.x; p.y=vec_tmp1.y; p.z=vec_tmp1.z;

	vec_tmp.x=p.xd; vec_tmp.y=p.yd; vec_tmp.z=p.zd;
	rot_x(vec_tmp, &vec_tmp1, -AXIAL_TILT);
	p.xd=vec_tmp1.x; p.yd=vec_tmp1.y; p.zd=vec_tmp1.z;

	pos1[0] = p.x;
	pos1[1] = p.y;
	pos1[2] = p.z;

	bary2obs (pos1,pob, pos2,&lighttime); /* bary_to_geo (pos1,pob, pos2,&lighttime); */
	t3 = tjd0 - lighttime;

	if(niter>31){
	  break;
	}

  } while (fabs (t3 - t2) > 1.0e-10);

  vector2radec (pos2, ra, dec);

  *ra *= 15.0;
  *dec *= 1.0;
  *x = p.x;
  *y = p.y;
  *z = p.z;
  *vx = p.xd;
  *vy = p.yd;
  *vz = p.zd;
  *gdist = lighttime*C_AUDAY;

}
