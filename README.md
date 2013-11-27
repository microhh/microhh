MicroHH
=======
MicroHH is a computational fluid dynamics code made for Direct Numerical Simulation (DNS) and Large-Eddy Simulation of turbulent flows in the atmospheric boundary layer.

MicroHH is created by Chiel van Heerwaarden (chielvanheerwaarden@gmail.com) and Thijs Heus (thijsheus@gmail.com).

Requirements
============
In order to compile MicroHH you need:
[x] C++ compiler
[x] MPI2 implementation (optional, but recommended)
[x] FFTW3 libraries
[x] NetCDF4
[x] CMake
[x] Doxygen (optional)
[x] Python + python-netcdf4 + numpy + matplotlib (optional for running example cases)

Compilation of the code
=======================
First, enter the config directory: 

cd config

Here, you find a potential series of settings with the extension .cmake for different systems. Check whether your system is there. If not, create a file with the correct compiler settings and the proper location for all libraries. Then, copy your system file to default.cmake. Let us assume your system is Ubuntu:

cp ubuntu.cmake default.cmake

Then, go back to the main directory and create a subdirectory with an arbitrary name in which you will compile the code. Let us assume this directory is called "build":

mkdir build
cd build

From this directory, run cmake with the suffix .. to point to the parent directory:

cmake ..

This should trigger the build system and create the make files, if the default.cmake file contains the correct settings. In case this works correctly, you can start the compilation of the code and create the "microhh" executable:

make -j

Running an example case
=======================
To start one of the included test cases, open the directory "cases". Here, a collection of test cases has been included. In this example, we start the drycblles case, a simple large-eddy simulation of a dry convective boundary layer.

cd cases/drycblles

First, we have to create the vertical profiles for our prognostic variables:

python drycbllesprof.py

Then, we have to copy or link the microhh executable to the current directory. Here we assume the executable is in the build directory that we have created before.

cp ../../build/microhh .

Now, we can start microhh in init mode to create the initial 3d fields.

./microhh init drycblles

If everything works out properly, a series of files has been created. The model can be started now following:

./microhh run drycblles

This will take some time. Now, a statistics file called drycblles.0000000.nc has been created. You can open this file with your favorite plotting tool, or run some example plots using the provided plotting script that uses python and matplotlib. This is most easily done in interactive python:

ipython
run drycbllesstats

Happy MicroHHing!
