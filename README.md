# Surface Multigrid via Intrinsic Prolongation
<img src="./assets/teaser.jpg" width="100%">

Public code release for [Surface Multigrid via Intrinsic Prolongation](https://www.dgp.toronto.edu/projects/intrinsic-prolongation/). For more details, please refer to:

**Surface Multigrid via Intrinsic Prolongation**<br>
[Hsueh-Ti Derek Liu](https://www.dgp.toronto.edu/~hsuehtil/), [Jiayi Eris Zhang](https://eriszhang.github.io/), [Mirela Ben-Chen](https://mirela.net.technion.ac.il/), and [Alec Jacobson](https://www.cs.toronto.edu/~jacobson/)<br>
ACM Transaction on Graphics (Proceedings of SIGGRAPH 2021)<br>
**[[Paper](http://www.dgp.toronto.edu/~hsuehtil/pdf/surfMG_65mb.pdf)]** **[[ArXiv](https://arxiv.org/abs/2104.13755)]** **[[Project Page](https://www.dgp.toronto.edu/projects/intrinsic-prolongation/)]**

## Installation
To get started, clone this repository *recursively*
```
git clone --recursive https://github.com/HTDerekLiu/surface_multigrid_code.git
```
On all platforms, we assume you have installed cmake and a modern c++ compiler on Mac OS X, Linux, or Windows.

## Layout
The main folder contains 6 separate examples that demonstrate some core functionalities and typical usage of our code. All of them have a similar directory and file layout:
```
cmake/
  CMakeLists.txt
README.md
main.cpp
```
+ `01_single_SSP/`: visualize our prolongation operator by mapping the fine mesh vertices onto the coarse mesh.
+ `02_mg_hierarchy/`: show the construction of our multigrid hierarchy and visualize the corresponding prolongation operators between different levels.
+ `03_mg_solver/`: show the usage of our multigrid solver on surface meshes with boundaries.
+ `04_mg_solver_nobd/`: show the usage of our multigrid solver on surface meshes without boundaries.
+ `05_example_mean_curvature_flow/`: show the usage of our multigrid solver in mean curvature flow.
+ `06_example_balloon_sim/`: show the usage of our multigrid solver in balloon simulation.

And they share a common `src` folder for source code and a `meshes` folder for input meshes.

## Compilation
Inside each subfolder, for example `01_single_SSP`, compile in release mode using the following typical cmake/make build routine:
```
cd 01_single_SSP
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```
If all goes well, you should be able to find and run the executable `main_bin` directly with no arguments.

## Usage
We aim to make our multigrid solver as easy to use as possible. To use our multigrid solver for linear systems defined on a fixed surface mesh, first call `mg_precompute` to construct the multigrid hierarchy. Then, similar to [libigl](https://github.com/libigl/libigl/)'s style of calling direct solvers [min_quad_with_fixed](https://github.com/libigl/libigl/blob/main/include/igl/min_quad_with_fixed.h), our multigrid solver can be used by calling `min_quad_with_fixed_mg_precompute` given a left-hand-side system matrix `A` and then calling `min_quad_with_fixed_mg_solve` to solve the system with the right-hand-side `b`. Please refer to the above examples for more details.