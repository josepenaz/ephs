#include <string.h>
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TWOPI 6.283185307179586476925287;
/*#define C 173.14463348*/
#define G0 39.476926421373015
#define Msun 1.0000059769998062e+00
#define SSMASS  1.00134       /*Total SS mass*/
#define DAY	(1./365.25)	/*Julian day, 86400 s*/
#define AXIAL_TILT (23+26/60.+21.406/3600.)*PI/180. /* tilt between equatorial and ecliptic coordinates ec=rot_x(eq,A_T); eq=rot_x(ec,-A_T) */
#define ASEC2RAD 4.848136811095359935899141e-6 /* Constant conversion from arcsec to rad. From NOVAS. */
#define C_AUDAY 173.1446326846693 /* Speed of light in AU/day. From NOVAS. */
#define DEG2RAD 0.017453292519943296
#define RAD2DEG 57.295779513082321

typedef struct {
  double x, y, z, xd, yd, zd;
} State;

typedef struct {
  double x, y, z;
} Vector;

void state_equatorial_heliocentric(double jd, double tjdepoch, double a, double e, double incl, double longnode, double argperi, double meananom, double *x, double *y, double *z, double *vx, double *vy, double *vz);
void apparent_coords_equatorial_heliocentric(double jd, double pobsx,  double pobsy,  double pobsz, double tjdepoch, double a, double e, double incl, double longnode, double argperi, double meananom, double *ra, double *dec, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *gdist);
void params_from_state(double x, double y, double z, double vx, double vy, double vz, double *a, double *e, double *i, double *longnode, double *argperi, double *meananom);
