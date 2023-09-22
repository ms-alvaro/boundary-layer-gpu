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

Provided that FFTW is installed under `<path-to_fftw>` and that `ifort`
and `mpiifort` are available commands, follow these steps:

1. Add the include and lib paths of your local installation of FFTW to your `~/.bashrc`:
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

### Advanced compilation

In case you want to make some modification to the source code. A clean way to 
do it is by creating and empty directory (`my_modifications` in this example)
elsewhere.

Then, copy [Makefile](Makefile) into `my_modifications` and edit the line
```Makefile
SRC = pwd
```
replacing `pwd` by `<path>/<name>`


These steps will generate an executable file `boundary_layer_<version>`
in `<path>/<name>` directory.

> NOTE: If you are in supercloud.mit.edu, intel compilers can be loaded as:
> ```bash
> module load intel-oneapi/2023.1
> ```

## Setting-up a simulation

To set up a simulation, `BLseparation` needs and input file. A template of the
input file can be found in [input_parameters.turbb](input_parameters.turbb). 
You should copy this file and edit the local copy.

> NOTE: If `inflow_flag` > 1, a binary input file has to be specified 
> as `inflow_file`. For that, you can use [generate_inflow_file.py](pyfiles/generate_inflow_file.py).


## Run a simulation

Once you have defined your input file, `<myinput_file>.turbb`, 
the simulation is launched as
```bash
mpirun -np <np> ./boundary_layer_<version> -i <myinput_file>.turbb > output
```
where `<np>` is the number of processors.


## Post-processing

After executing the code, `BLseparation` generates two different output files: the instantaneous
flow field, and a `txt` file with statistics of the flow.

Python scripts to post-process these files can be found in [tests/pyfiles](tests/pyfiles).


