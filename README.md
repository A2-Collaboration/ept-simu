# End-point tagger simulation

This repository contains somewhat structured resources to simulate the
electron flight paths according inside the measured magnetic field of
the end-point tagger. It also has [some documentation](doc) in German.
All this was written by Juergen Ahrens, and put into git and
cmake'ified by Andreas Neiser.

## Installation and Usage

To build this project, you need
  * GNU compiler version >4.7.2 with gfortran
  * cmake version >2.6

Then do a standard out-of-source build by executing the following
inside the projects directory:

    mkdir build
    cd build
    cmake ..
    make

You should find various executables inside the build directory then.
They probably don't work out of the box.
