# Aleph

Aleph is a fast ephemeris producer joining [`astropy`](https://www.astropy.org/) and [`rebound`](https://rebound.readthedocs.io/en/latest/) with orbital parameters from [MPC](https://minorplanetcenter.net/)'s [MPCORB.DAT](https://minorplanetcenter.net/iau/MPCORB.html) or Lowell's [astorb.dat](https://asteroid.lowell.edu/main/astorb/). It allows you to get multi-epoch ephemeris of any asteroid in the catalogue or single-epoch ephemeris of _all_ asteroids in a given field of view.

You can choose if ephemeris are calculated using the simple 2-Body equations (fast but not very accurate) or integrating the N-Body problem including the Sun and the planets with their masses and each asteroid as massless particles (much more accurate but more time expensive). To speed up the ephemeris calculation of so many bodies, almost all functions allows parallel integrations.

## Instalation

You can clone (or download) this repository and install it from there:
``` {.sourceCode .bash}
git clone https://github.com/josepenaz/ephs.git
cd ephs
python3 setup.py install
```
You can also install it using `pip` from the cloned repository:
``` {.sourceCode .bash}
git clone https://github.com/josepenaz/ephs.git
cd ephs
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m build
python3 -m pip install dists aleph-X.X.tar.gz
```

## Getting started - DataBase class

The most basic class in Aleph is `DataBase.OrbParams` which reads and saves an orbital parameter catalogue. By default, it will try to read [MPCORB.DAT](https://minorplanetcenter.net/iau/MPCORB.html). The first time you start an `OrbParams` object, it will raise an error since the catalogue hasn't been downloaded. To avoid this, you can set `update=True` and it will download the catalogue and save it in the Aleph's directory:
``` {.sourceCode .python}
from aleph import DataBase
orbs = DataBase.OrbParams(update=True)
```
Use `update=True` only when you want to download the catalogue. Once it is download and saved, it will automatically read it each time you initialize an `OrbParams` object. Do not use `update=True` unless you want to download it again.
If you want to use another version of MPCORB.DAT, you can give give the path to the file you want to open:
``` {.sourceCode .python}
orbs = DataBase.OrbParams(filename='path/MPCORB.DAT')
```

`OrbParams` can also read Lowell's [astorb.dat](https://asteroid.lowell.edu/main/astorb/) by setting `service='Lowell'`. Since FTP access to astorb.dat is not allowed, you must download it manually and specify the path to the file:
``` {.sourceCode .python}
orbs = DataBase.OrbParams(service='Lowell', filename='path/astorb.dat')
```

Once you have loaded the catalogue with `OrbParams` you can access to all asteroid's orbital parameters with the `asts` attribute, which is a list of dictionaries. Each dictionary saves the orbital parameters of each asteroid. For example, you can access to Ceres (the first body of the catalogue):
``` {.sourceCode .python}
ceres = orbs.asts[0]
```
You also can access to the epochs of the asteroids' mean anomaly with the attribute `dicc_epochs_M`. That atribute is important because it tells you how much the orbits need to be integrated to get the ephemeris (the integration starts at that epoch and ends at the epoch you want the ephemeris to be).

## Getting Ephemeris

To look for all the asteroids that are visible in a given field of view you first have to create a `Query` instance:
``` {.sourceCode .python}
from aleph.Query import Query
q = Query()
```
By default, it will load the MPCORB.DAT database. If you want to load the astorb.dat catalogue, you need to set `service='Lowell'` and provide the catalogue itself with `filename='path/astorb.dat'`:
``` {.sourceCode .python}
from aleph.Query import Query
q = Query(service='Lowell', filename='path/astorb.dat')
```
### Virtual Observatory: Using the catalogues
With the `Query` instance initialized you can _query_ for the ephemeris you want. In the next example, you ask for the asteroids in a field of 0.5 degrees of radius centered in the equatorial coordinates (0.5, 0.5) degrees, observed in July 25th, 2020 at Palomar Observatory. And to speed up things, we tell `Query` to compute the ephemeris in 5 parallel process:
``` {.sourceCode .python}
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy import units as u
 
palomar = EarthLocation.of_site('Palomar')
field_center = SkyCoord(0.5*u.deg, 0.5*u.deg)
field_radius = 0.5*u.deg
epoch = Time('2020-07-25')

ephs = q.query_mixed_cat(field_center, field_radius, epoch=epoch, observer=palomar, njobs=5)
```
The function `query_mixed_cat` starts integrating using the equations of the 2-Body problem, it takes all bodies fallen in a similar but wider area than the field of view and integrates them along with the 8 planets of our Solar System (this integration is performed with the `IAS15` integrator of `rebound`).

You can use the 2-Body or the N-Body integrations just by using the `query_2b_cat` and `query_nb_cat` functions respectively.

### Creating your own Data Base of integrated orbits
Lets say you need many ephemeris far from the mean anomaly epoch of your database. This means that for every time you run a `q.query_mixed_cat` or a `q.query_nb_cat` all asteroids will be integrated again and again the same period of time till the ephemeris epoch each time you call those methods. So the idea is you integrate all asteroids only once to an epoch close to the epoch you need to query, so every time you make a query the integrations will be much shorter. This can only be done using the astorb.dat catalogue. In the next example you will create the database (which is done just by saving the orbital parameters of the catalogue) and then you will integrate all bodies and save their _states_ (positions and velocities) at a far epoch. Finally, you can compute and save the orbital parameters that corresponds to the states you just compute.
``` {.sourceCode .python}
from astropy.time import Time
from aleph import DataBase

aleph = DataBase.OrbParams(service='Lowell', filename='path/astorb.dat')
aleph.create_states_database() # You can provide 'sqlfile' with
                               # the name of the file where you
                               # want to save the database.

# Lets say you want to integrate till two epochs:
integ_epochs = Time(['2022-01-02','2022-06-02']).jd # You need to provide the Julian Day
aleph.adding_states_into_db(integ_epochs, njobs=5)  # Remember to provide 'sqlfile' if
                                                    # you use it before
aleph.adding_params_into_db()    # Again, remember 'sqlfile' if you use it before
```
That process can take a long time, but it will greatly increase the execution time of any ephemeride queried around the integrated epochs.

### Using your own Data Base of integrated orbits
Once the database is created, you can query for ephemeris using the `query_2b_db`, `query_nb_db` and `query_mixed_db` functions (that work in the same manner than the functions above). Remember that those functions work only with astorb.dat:
``` {.sourceCode .python}
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy import units as u
from aleph.Query import Query

q = Query(service='Lowell', filename='path/astorb.dat')
 
palomar = EarthLocation.of_site('Palomar')
field_center = SkyCoord(0.5*u.deg, 0.5*u.deg)
field_radius = 0.5*u.deg
epochs = Time(['2022-01-01','2022-01-02','2022-01-03',
               '2022-06-01','2022-06-02','2022-06-03'])

all_ephs = []
for epoch in epochs:
    ephs = q.query_mixed_cat(field_center, field_radius, epoch=epoch, observer=palomar, njobs=5)
    all_ephs.append(ephs)
```

### Getting multi-epoch ephemeris of a single body
You can get ephemeris of a single body at as many epochs you want. This ephemeris are computed integrating the N-Body problem (Sun, planets and the body).
``` {.sourceCode .python}
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy import units as u
from aleph.Query import Query

q = Query(service='Lowell', filename='path/astorb.dat')
 
palomar = EarthLocation.of_site('Palomar')
epochs = Time(['2022-01-01','2022-01-02','2022-01-03',
               '2022-06-01','2022-06-02','2022-06-03'])

ephs = q.asteph(0, epochs=epochs, observer=palomar)
```
To use this function, you need to create an Aleph database (see previous section), although you do not need to integrate all the asteroids to any epoch (running `aleph.create_states_database()` is enough).
