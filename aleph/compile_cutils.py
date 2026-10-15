from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(
    """
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
    """
)

ffibuilder.set_source(
    "aleph._cutils",
    """
    #include "cutils.h"
    """,
    sources=["aleph/cutils.c"],
    include_dirs=["aleph"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)