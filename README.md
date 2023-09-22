# BLseparation

This repo contains the fortran code to solve an incompressible boundary layer 
with a suction/blowing boundary condition at the top wall.

The code solves the Navier-Stokes equations of an incompressible flow
using finite differences in a staggered grid with a RK3 temporal scheme.

The code has been written by Adrian Lozano-Duran.


## Download

To download this repo, type: 

```bash  
git clone git@github.mit.edu:Computational-Turbulence-Group/BLseparation.git <path>/<name>
```
where `<path>` is the path where you want to install the code, and `<name>` is the name of the repo. If `<name>` is not given, the default *BLseparation* is used.


## Pre-requisites

`BLseparation` requires the [FFTW library](https://www.fftw.org) and the intel
compilers `ifort` and `mpiifort` [intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html#gs.6296q9).


## Compilation
Provided that the `build` folder of fftw is found in `<path-to_fftw>` and that `ifort`
and `mpiifort` are available commands, follow these steps:

1. Add the include and lib paths of your local installation of fftw to your `~/.bashrc`:
    ```bash
    export FFTW_INCLUDE_DIR=<path-to-fftw>/include
    export FFTW_LIBRARY_DIR=<path-to-fftw>/lib
    ```
2. Load the new variables:
    ```bash
    source ~/.bashrc
    ```
3. Go to `<path>/<name>` and type:
    ```bash
    make
    ```
4. (optional) For debugging purposes, you can also compile the code using
    ```bash
    make clean
    make debug
    ```

    NOTE: If you are in supercloud, intel compilers can be loaded as:
    ```bash
    module load intel-oneapi/2023.1
    ```

