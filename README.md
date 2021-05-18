
  (BeAST) solenoid magnetic field map
  ===================================

  Was originally written to import BeAST magnetic field map. As long as 
the new 3T solenoid field maps are provided as a 4-column ASCII file 
(r, z, br, bz) on a 2cm NZxNR grid, will work with those ones as well. 

  The repository contains ASCII files with the field maps and a 
C++ class to handle it.


Installation
------------

```
  git clone https://github.com/eic/BeastMagneticField.git
```

Compiling
---------

```
  cd BeastMagneticField && mkdir build && cd build
  # Assume installation in the current build directory; replace "=." by a different path if needed;
  cmake -DCMAKE_INSTALL_PREFIX=. -Wno-dev ..
  make install
```

Usage (standalone executable)
-----------------------------

```
  # It is assumed that libbmf.so installation directory is in the LD_LIBRARY_PATH (Linux); 
  # the command below will import the ASCII field map located in the ../data directory (check the path!) 
  # and return the field value (radial and longitudinal components) at [R=3.0cm , Z=100.0cm];

  # assume we are still in the "build" directory;
  ./bmf-main ../data/EIC_Magnetic_Field_Map_2021_05_07_radial_coords_[cm]_[T].120000.lines.Bmap 3.0 100.0
```

Usage (shared library)
----------------------

```
  # Refer to the source/main.cc file; the essential lines are 
 
    auto bmf = new BeastMagneticField(<ASCII-map-file-location>);
    bool ret = bmf->GetFieldValue(r, z, br, bz);

  # include/BeastMagneticField.h has a couple of methods to control the GetFieldValue() call behaviour;
```

Usage (ROOT)
------------

```
  # Make sure ROOTSYS is defined before running cmake;
  #
  # Refer to the visualization scripts/viewer.C as an example; assume we are still in the "build" directory;

  root -l '../scripts/viewer.C("../data/EIC_Magnetic_Field_Map_2021_05_07_radial_coords_[cm]_[T].120000.lines.Bmap")' 

```

  Original BeAST field map
  ========================

  The field map was calculated for a coil configuration tuned to 
provide a relatively constant 3T field inside the central part of 
the bore (+/-1.0 m along the beam line and up 0.8 m in radius)
, as well as to align the field lines in the forward gaseous RICH 
volume with the charged particle trajectories originating from 
the IP. 

  Open source Elmer calculator was used. Field map is azimuthally
symmetric and provided in the range +/-5.0 m along the beam 
direction and up to a radius of 2.5 m (however only the area up 
to 137 cm in radius is considered to be the "inner bore").