#ifndef CUTILS_H
#define CUTILS_H

#include <math.h>

#define PI 3.14159265358979323846
#define TWOPI 6.283185307179586476925287
#define G0 39.476926421373015
#define Msun 1.0000059769998062e+00
#define SSMASS 1.00134
#define DAY (1.0 / 365.25)
#define AXIAL_TILT ((23.0 + 26.0 / 60.0 + 21.406 / 3600.0) * PI / 180.0)
#define ASEC2RAD 4.848136811095359935899141e-6
#define C_AUDAY 173.1446326846693
#define DEG2RAD 0.017453292519943296
#define RAD2DEG 57.295779513082321

typedef struct {
    double x, y, z, xd, yd, zd;
} State;

typedef struct {
    double x, y, z;
} Vector;

void state_equatorial_heliocentric(
    double jd,
    double tjdepoch,
    double a,
    double e,
    double incl,
    double longnode,
    double argperi,
    double meananom,
    double *x,
    double *y,
    double *z,
    double *vx,
    double *vy,
    double *vz
);

void apparent_coords_equatorial_heliocentric(
    double jd,
    double pobsx,
    double pobsy,
    double pobsz,
    double tjdepoch,
    double a,
    double e,
    double incl,
    double longnode,
    double argperi,
    double meananom,
    double *ra,
    double *dec,
    double *x,
    double *y,
    double *z,
    double *vx,
    double *vy,
    double *vz,
    double *gdist
);

void params_from_state(
    double x,
    double y,
    double z,
    double vx,
    double vy,
    double vz,
    double *a,
    double *e,
    double *i,
    double *longnode,
    double *argperi,
    double *meananom
);

#endif /* CUTILS_H */