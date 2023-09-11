# FV1D
Small 1D finite volume code for experimentation purposes. Header-only + 1 main file.

## Build

Assuming you have `cmake` (3.16+), a C++ compiler and hdf5 (1.8+) installed on your computer, at the root of the project :

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
```

Should do the configuration. Then build the project with `make`. Once built, the `build`folder should have an executable `fv1d`

## Running and analysing

To run a simulation, just type `./fv1d INIFILE` zith `INIFILE` a `.ini` file for configuration. You can find some examples in the `settings` folder at the root of the project. Once the program has finished running, you can analyse the results using one of the python scripts in the `python` folder at the root of the project.
