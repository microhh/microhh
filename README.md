MicroHH
-------
[![Travis](https://api.travis-ci.org/microhh/microhh.svg?branch=master)](https://travis-ci.org/microhh/microhh) [![Documentation Status](https://readthedocs.org/projects/microhh/badge/?version=latest)](https://microhh.readthedocs.io/en/latest/?badge=latest)

MicroHH is a computational fluid dynamics code made for Direct Numerical Simulation (DNS) and Large-Eddy Simulation of turbulent flows in the atmospheric boundary layer. The code is written in C++.

MicroHH is hosted on GitHub (http://github.com/microhh). Here, the latest version of the source code can be found, as well as all releases. Bug notifications and fixes are always welcome.

MicroHH is described in detail in [Van Heerwaarden et al. (2017)](http://www.geosci-model-dev-discuss.net/gmd-2017-41/#discussion). In case you decide to use MicroHH for your own research, the developers would appreciate to be notified and kindly request to cite their reference paper. The version described in the reference paper has been assigned a DOI via [Zenodo](https://zenodo.org).

[![DOI](https://zenodo.org/badge/14754940.svg)](https://zenodo.org/badge/latestdoi/14754940)

Requirements
------------
In order to compile MicroHH you need:
* C++ compiler
* FFTW3 libraries
* Boost libraries
* NetCDF4
* CMake
* MPI2/3 implementation (optional for MPI support)
* CUDA (optional for GPU support)
* Python + numpy + python-netcdf4 (optional for running example cases)
* Ipython + python-netcdf4 + matplotlib (optional for plotting results example cases)

Compilation of the code
-----------------------
First, enter the config directory: 

    cd config

Here, you find a potential series of settings with the extension .cmake for different systems. Check whether your system is there. If not, create a file with the correct compiler settings and the proper location for all libraries. Then, copy your system file to default.cmake. Let us assume your system is Ubuntu:

    cp ubuntu.cmake default.cmake

Then, go back to the main directory and create a subdirectory with an arbitrary name in which you will compile the code. Let us assume this directory is called "build":

    cd ..  
    mkdir build  
    cd build   

From this directory, run cmake with the suffix .. to point to the parent directory where the CMakeLists.txt is found. This builds the model without Message Passing Interface (MPI) and CUDA support.

    cmake ..

In case you prefer to enable either MPI or CUDA, run INSTEAD of the previous command:
    
    cmake .. -DUSEMPI=TRUE

or

    cmake .. -DUSECUDA=TRUE

(Note that once the build has been configured and you wish to change the `USECUDA` or `USEMPI` setting, you must delete the build directory or create an additional empty directory from which `cmake` is run.)

With the previous command you have triggered the build system and created the make files, if the `default.cmake` file contains the correct settings. Now, you can start the compilation of the code and create the microhh executable with:

    make -j

Your directory should contain a file named `microhh` now. This is the main executable.

Running an example case
-----------------------
To start one of the included test cases, go back to the main directory and  open the directory `cases`. Here, a collection of test cases has been included. In this example, we start the `drycblles` case, a simple large-eddy simulation of a dry convective boundary layer.

    cd cases/drycblles

First, we have to create the vertical profiles for our prognostic variables:

    python drycblles_input.py

Then, we have to copy or link the `microhh` executable to the current directory. Here we assume the executable is in the build directory that we have created before.

    cp ../../build/microhh .

Now, we can start `microhh` in initialization mode to create the initial fields:

    ./microhh init drycblles

If everything works out properly, a series of files has been created. The model can be started now following:

    ./microhh run drycblles

This will take some time. Now, a statistics file called `drycblles.default.0000000.nc` has been created. You can open this file with your favorite plotting tool, or run some example plots using the provided plotting script that uses Python and matplotlib. This is most easily done in interactive python:

    ipython  
    run drycbllesstats

This should show you a set of basic plots. Congratulations, you have just completed your first run of MicroHH.

Happy MicroHHing!

