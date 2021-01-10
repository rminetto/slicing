# An Optimal Algorithm for 3D Triangle Mesh Slicing 

This repository has the source codes (versions in c++, python and opengl) of an algorithm for slicing an unstructured triangular mesh model by a series of parallel planes. We prove that the algorithm is asymptotically optimal: its time complexity is O(n log k + k + m) for irregularly spaced slicing planes, where n is the number of triangles, k is the number of slicing planes, and m is the number of triangle-plane intersections segments. The time complexity reduces to O(n + k + m) if the planes are uniformly spaced or the triangles of the mesh are given in the proper order. It also contains an asymptotically expected lineartime algorithm for constructing a set of polygons from the unsorted lists of line segments produced by the slicing step. 

![](/assets/slicing.jpg)

# Download STL models

[Drive] https://drive.google.com/file/d/1cF7Yokc5vJmBHYp1YO1NI0Ny2kfE76mR/view?usp=sharing

# Citing

```
@ARTICLE{9234742,
  author={R. {Minetto} and N. {Volpato} and J. {Stolfi} and R.M.M.H. {Gregori} and M.V.G. da {Silva}},
  journal={Computer-Aided Design (Elsevier)}, 
  title={An Optimal Algorithm for 3D Triangle Mesh Slicing}, 
  year={2017},
  volume={92},
  number={1},
  pages={1-10},
  doi={http://dx.doi.org/10.1016/j.cad.2017.07.001},
  url = {http://www.sciencedirect.com/science/article/pii/S0010448517301215}
}
'''
