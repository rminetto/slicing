import numpy as np
from stl import mesh as mesh_data

## Triangle class for STL.
#
#  Class for a triangle.
class stltriangle:
    
    ## Constructor.
    # @param v0 First vertex.
    # @param v1 Second vertex.
    # @param v2 Third vertex.
    # @param normal Normal vector.
    def __init__(self, v0, v1, v2, normal):

        ## First vertex of the triangle.
        self.v0     = v0
        ## Second vertex of the triangle.
        self.v1     = v1
        ## Third vertex of the triangle.
        self.v2     = v2
        ## Normal of the triangle.
        self.normal = normal

    ## Minimum z coordinate.
    #
    # Return minimum z coordinate.
    #
    # @return Minimum z.
    def min_z(self):

        min_coordinate = self.v0[2]

        if min_coordinate > self.v1[2]:
            min_coordinate = self.v1[2]
        if min_coordinate > self.v2[2]:
            min_coordinate = self.v2[2]

        return min_coordinate
    
    ## Maximum z coordinate.
    #
    # Return maximum z coordinate.
    #
    # @return Maximum z.
    def max_z(self):

        max_coordinate = self.v0[2]

        if max_coordinate < self.v1[2]:
            max_coordinate = self.v1[2]
        if max_coordinate < self.v2[2]:
            max_coordinate = self.v2[2]

        return max_coordinate



## Mesh class for STL files.
#
#  Class for a mesh read from an STL file.
class stlmesh:
    
    ## Constructor.
    # @param stl_file STL file name.
    def __init__(self, stl_file):

        ## Number of tringles in the mesh.
        self.ntriangles = 0 
        ## List of triangles.
        self.triangles = []

        # Read STL file.
        data = mesh_data.Mesh.from_file(stl_file)

        # Create list of triangles.
        for v0, v1, v2, normal in zip(data.v0, data.v1, data.v2, data.normals):
            self.triangles.append(stltriangle(v0,v1,v2,normal))

        # Set number of triangles.
        self.ntriangles = len(self.triangles)

    ## Return a vector of maximum values for coordinates in the mesh.
    # @return Array of maxima.
    def max_coordinates(self):
        
        Mval = np.ones(3, dtype='float32')*np.finfo(np.float32).min
        for i in range(self.ntriangles): 
            Mval = np.maximum(Mval,self.triangles[i].v0)
            Mval = np.maximum(Mval,self.triangles[i].v1)
            Mval = np.maximum(Mval,self.triangles[i].v2)

        return Mval
    
    ## Return a vector of minimum values for coordinates in the mesh.
    # @return Array of minima.
    def min_coordinates(self):
        
        Mval = np.ones(3, dtype='float32')*np.finfo(np.float32).max
        for i in range(self.ntriangles): 
            Mval = np.minimum(Mval,self.triangles[i].v0)
            Mval = np.minimum(Mval,self.triangles[i].v1)
            Mval = np.minimum(Mval,self.triangles[i].v2)

        return Mval

    ## Create an array of vertex and normal data for drawing triangles in OpenGL.
    #  
    #  @param bmin Minimim bounding box value.
    #  @param bmax Maximum bounding box value.
    #  @return Data Array.
    #  @return Number of vertices.
    #
    #  The returned array is contains the vertices and normals of the triangles structures as:
    #  
    #  | triangle 1                  | triangle 2                  |...| triangle n                  |
    #
    #  such that the respective vertices and normals are intertwined as :
    #
    #  |v0|normal|v1|normal|v2|normal|v0|normal|v1|normal|v2|normal|...|v0|normal|v1|normal|v2|normal|
    #
    #  The coordinates of each vertex is normalized to the range [bmin,bmax].
    #  OpenGL draws coordinates in the range [-1,1]. Therefore, \f$-1 \leq bmin <
    #  bmax \leq 1\f$.
    #
    # @return Array.
    def OpenGLData(self,bmin,bmax):

        nvertices = 3*self.ntriangles

        # Data array with triangles stored as: |v0|normal|v1|normal|v2|normal|
        # so that there are 3 vertices + 3 normals for n triangles of 3
        # coordinates triangles: (3+3)*3*n.
        data = np.zeros((3+3)*3*self.ntriangles, dtype='float32')

        # Maxima and minima for normalization.
        mesh_max = self.max_coordinates()
        mesh_min = self.min_coordinates()

        # Fill normalized data.
        j = 0
        for i in range(self.ntriangles):
            data[j:j+3] = (self.triangles[i].v0-mesh_min)/(mesh_max-mesh_min)*(bmax-bmin)+bmin
            j += 3
            data[j:j+3] = self.triangles[i].normal
            j += 3
            data[j:j+3] = (self.triangles[i].v1-mesh_min)/(mesh_max-mesh_min)*(bmax-bmin)+bmin
            j += 3
            data[j:j+3] = self.triangles[i].normal
            j += 3
            data[j:j+3] = (self.triangles[i].v2-mesh_min)/(mesh_max-mesh_min)*(bmax-bmin)+bmin
            j += 3
            data[j:j+3] = self.triangles[i].normal
            j += 3

        return data, nvertices

