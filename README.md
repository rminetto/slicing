# An Optimal Algorithm for 3D Triangle Mesh Slicing 

This repository has the source codes (versions in c++, python and opengl) of an algorithm for slicing an unstructured triangular mesh model by a series of parallel planes. We prove that the algorithm is asymptotically optimal: its time complexity is O(n log k + k + m) for irregularly spaced slicing planes, where n is the number of triangles, k is the number of slicing planes, and m is the number of triangle-plane intersections segments. The time complexity reduces to O(n + k + m) if the planes are uniformly spaced or the triangles of the mesh are given in the proper order. It also contains an asymptotically expected lineartime algorithm for constructing a set of polygons from the unsorted lists of line segments produced by the slicing step. 
