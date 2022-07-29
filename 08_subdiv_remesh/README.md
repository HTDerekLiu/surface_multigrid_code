# Subdivision Remeshing
## Compilation
To run this example, compile in release mode using the following typical cmake/make build routine:
```
cd 08_subdiv_remesh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```
If all goes well, you should be able to find and run the executable `main_bin` directly with no arguments and you will see some `output_s?.obj` files.

## Demo
<img src="../assets/08.png" width="100%">

This demo shows how one can use our successive self-parameterization to do subdivision remeshing by upsampling a base mesh (left) and the map to the input mesh.