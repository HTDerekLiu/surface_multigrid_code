clc; clear all; close all;

[V,F] = readOBJ("hilbert_cube_remesh.obj");
[T,N,B,I] = per_vertex_frames(V,F);
avgE = avgedge(V,F);
nV = size(V,1);
V = V + 0.1 * avgE * (T+B);

writeOBJ('test.obj',V,F)