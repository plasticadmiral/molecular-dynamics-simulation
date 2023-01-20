


This project was created as a part of the [Molecular Dynamics Course](https://imtek-simulation.github.io/MolecularDynamics/) at Uni Freiburg. An explanation of the project structure can be found on the report. The link to the large data files is the following.

https://mega.nz/folder/CeBj1SBK#RaQxGveHUBE82Nj6cY-f2Q

The visualization of the project data can be found in `results/plots_and_figures.ipynb` .
A part of the source files used in this project were provided by the lecture content. The complete explanation and bibliography can be found in the report. This project has been build using [Ubuntu on WSL](https://ubuntu.com/wsl), therefore, the possibility of running all sections of this project using alternate configurations/operating systems is not guaranteed. The data for the parallel section ([MPI](https://www.open-mpi.org/) operations) of this project was generated using [BwUniCluster 2.0](https://wiki.bwhpc.de/e/BwUniCluster2.0).

## Tools used:
The datastructures are provided by the [Eigen](https://eigen.tuxfamily.org/) library and tests are run on the implemented features using [Googletest](https://github.com/google/googletest).  [Ovito](https://www.ovito.org/about/) is used to visualize the simulation. The plots are made with python using [Matplotlib](https://matplotlib.org/) and are managed using a [Jupyter notebook](https://jupyter.org/).

# Simulation Videos:
## Equilibriated Molecule
https://user-images.githubusercontent.com/49732300/212583771-2ec5e089-9911-46ab-99ff-7cbaf0bb3288.mov

## Phase change of an icosahedral gold structure
[![Phase change](https://vumbnail.com/790627243.jpg)](https://vimeo.com/790627243)

## Stretching of a gold nanowire
https://user-images.githubusercontent.com/49732300/212583811-12fa6379-0030-45c9-be37-b26629d84d59.mov

# Build 
To build the project, do the following set of commands in the project folder. 
```bash
mkdir build
cmake -DCMAKE_BUILD_TYPE=Release ../code
make
```
The executable should now be present in the build directory. The executable for the tests should be under `build/tests/` .


The same process follows to build on the unicluster. Except `-DEIGEN3_INCLUDE_DIR=<path>` and `-DMPI_BUILD=True` flags need to be set.

This project uses a modified version of [this](http://www.pas.rochester.edu/~wangyt/algorithms/ih/) code to create a [Makay Icosahedron Structure](http://doye.chem.ox.ac.uk/jon/structures/Morse/paper/node6.html). The files for this is present under the `data` folder. This project can be built in it's directory the following way. 

```bash
g++ -o output.out ih.C
```
and run using 
```bash
./output.out
```
if it does not run, execute the following command before running
```bash
chmod+x output.out
```
