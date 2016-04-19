# neoParSA: a scalable parallel simulated annealing library

`neoParSA` is a highly modularized parallel and serial simulated annealing library in C++.
It is designed to be both easy to use and easy to extend.

## 1. Dependencies
### MPI
`neoParSA` library proper requires an implementation of MPI version 2.0 or above.
Common ones are

  - [MPICH](https://www.mpich.org)
  - [Open MPI](https://www.open-mpi.org)

[Intel MPI Library] or other
vender specific MPI implementations also work.

[Intel MPI Library]: https://software.intel.com/en-us/intel-mpi-library

### Other dependencies
In addition to MPI, `neoParSA` library propery also requires
  - [CMake] >=2.8.0
  - [Boost] >=1.56
[CMake]: https://cmake.org/
[Boost]: http://www.boost.org/
The pattern formation model test problem `fly` requires an extra library
  - [GNU Scientific Library][gsl]
[gsl]: http://www.gnu.org/software/gsl/

## 2. Build
`neoParSA` uses CMake as a build system to allows cross platform, out of source
builds with flexible build options. It is always recommended to compile the library
in a separate directory. Suppose the source files are in the directory `neoParSA`,
it is recommended to create a subdirectory `build` inside.

    $ cd neoParSA
    $ mkdir build

If all dependencies are in their usual location, the build can be configured
without extra options.

    $ cd build
    $ cmake ..


You can specify options to for configuration:

    $ CC=icc CXX=icpc cmake -D CMAKE_BUILD_TYPE=Release -D BOOST_ROOT=$HOME/boost_1_59 ..

To compile just the library, do

    $ make parsa

To make the library and all the test problems, do

    $ make

## 3. Usage
