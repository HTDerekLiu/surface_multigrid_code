# Multigrid Hierarchy
## Compilation
To run this example, compile in release mode using the following typical cmake/make build routine:
```
cd 02_mg_hierarchy
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```
If all goes well, you should be able to find and run the executable `main_bin` directly with no arguments.

## Demo
<img src="./02.gif" width="100%">

Show the construction of our multigrid hierarchy and visualize the corresponding prolongation operators between different levels. Press `0` `1` `2` to toggle between different levels.

