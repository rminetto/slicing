## @file Slicer of mesh of triangles.
#
# Implements the algorithms from:
#
#  [1] Rodrigo Minetto, Neri Volpato, Jorge Stolfi, Rodrigo M.M.H. Gregori, Murilo
#  V.G. da Silva. An optimal algorithm for 3D triangle mesh slicing,
#  Computer-Aided Design, Volume 92, 2017.


import math
import numpy as np
import geotypes as gt

## Slicer class.
#
#  Class for slicing a mesh of triangles.
class slicer:
    
    ## Constructor.
    # @param T List of triangles in mesh.
    # @param P List of planes for slices (adaptive slicing).
    # @param delta Offset between planes for slices (regular slicing). 
    # @param srt True if planes are sorted, False otherwise.
    # @return List of planes/slices with polygons.
    def __init__(self, T, P, delta, srt):
        ## Final planes.
        self.planes = []
        ## Segments.
        self.segments = []
        ## Mesh triangles.
        self.T = T
        ## Slicing planes.
        self.P = P
        ## delta
        self.delta = delta
        ## Sorted.
        self.srt = srt
        ## Rounding eps.
        self.eps = 0.0004

        # Round coordinates of triangles to even multiples of eps.
        self.roundTrianglesEven()
        
        # Round delta to even multiple of eps.
        self.delta = self.mround(self.delta, self.eps, 2, 0)

        # Create planes.
        if delta > 0:
            # Create uniform planes to odd multiples of eps.
            self.uniformPlanes()
        

    ## Incremental slicing
    #
    # Algorithm 1 in [1].
    def incremental_slicing(self):

        # Split triangle list.
        L = self.build_triangle_list()

        # Plane sweeping
        A = []
        for i in range(len(self.P)):

            self.segments.append([])
            A.extend(L[i])

            t = 0
            while t < len(A):
                if self.T[A[t]].max_z() < self.P[i]:
                    A.pop(t)
                else:
                    I = self.compute_intersection(self.T[A[t]],self.P[i])
                    if I[0] != None and I[1] != None:
                        self.segments[i].append(I)
                    t += 1

        # Buil oriented polygons from segments.
        self.contour_construction()    


    ## List of triangles.
    #
    # Build list of triangles for each plane. Algorithm 2 in [1].
    #
    # @todo Implement pre-sorted part of algorithm.
    def build_triangle_list(self):
    
        L = []
        for i in range(len(self.P)+1):
            L.append([])

        # Uniform slicing.
        if self.delta > 0.0:
            for t in range(len(self.T)):
                if self.T[t].min_z() < self.P[0]: i = 0
                elif self.T[t].min_z() > self.P[-1]: i = len(self.P)
                else: i = math.floor((self.T[t].min_z()-self.P[0])/self.delta)+1
                L[i].append(t)
        # Pre-sorted triangles.
        #elif self.srt:
        # General case.
        else:
            for t in range(len(self.T)):
                i = self.binary_search(t)
                L[i].append(t)

        return L

    ## Binary seach
    #
    # Search for the plane of a triangle.
    def binary_search(self, t):

        min_z = self.T[t].min_z()

        if min_z > self.P[len(self.P)-1]:
            return len(self.P)
        
        if min_z < self.P[0]:
            return 0

        l = 0
        r = len(self.P)-1
        while r-l > 1:
            m = math.floor((l+r)/2.0)
        
            if min_z > self.P[m]:
                l = m
            else:
                r = m

        return r

    ## Compute intersection.
    #
    # Compute the segment intersection of a triangle and a plane.
    #
    # @param t Triangle.
    # @param p Plane.
    # @return endpoints of the intersection segment.
    def compute_intersection(self,t,p):
       
        intersection = [None,None]
        n = 0 # Index intersection positions to copy vertices.
        
        # First edge v0-v1.
        v = t.v1-t.v0
        if math.fabs(v[2]) > 0.0:
            alpha = (p-t.v0[2])/v[2]
            if alpha >= 0.0 and alpha <= 1.0:
                intersection[n] = gt.vertex3f(t.v0[0]+alpha*v[0], t.v0[1]+alpha*v[1], p)
                n += 1
         
        # Second edge v0-v2.
        v = t.v2-t.v0
        if math.fabs(v[2]) > 0.0:
            alpha = (p-t.v0[2])/v[2]
            if alpha >= 0.0 and alpha <= 1.0:
                intersection[n] = gt.vertex3f(t.v0[0]+alpha*v[0], t.v0[1]+alpha*v[1], p)
                n += 1

        # Third edge v1-v2.
        v = t.v2-t.v1
        if math.fabs(v[2]) > 0.0:
            alpha = (p-t.v1[2])/v[2]
            if alpha >= 0.0 and alpha <= 1.0:
                intersection[n] = gt.vertex3f(t.v1[0]+alpha*v[0], t.v1[1]+alpha*v[1], p)
                n += 1
                
        return intersection


    ## Create contours.
    #
    # Set oriented group of segments to form a polygon.
    def contour_construction(self):

        # Iterate planes.
        for i in range(len(self.segments)):

            # New plane to fill with polygons.
            self.planes.append([])
        
            # Round x and y coordinates (z is already rounded).
            #eps = 1.0/128.0
            for seg in self.segments[i]:
            #    seg[0].coord[0:2] = np.round(seg[0].coord[0:2]/eps)*eps
            #    seg[1].coord[0:2] = np.round(seg[1].coord[0:2]/eps)*eps
                seg[0].coord[0:2] = self.mround(seg[0].coord[0:2], self.eps, 2, 1)
                seg[1].coord[0:2] = self.mround(seg[1].coord[0:2], self.eps, 2, 1)

            # Fill hash table (dictionary).
            H = {}
            for seg in self.segments[i]:
                if np.linalg.norm(seg[0].coord-seg[1].coord) > 0.0001:
                    if seg[0].coord.tostring() in H:
                        H[seg[0].coord.tostring()].append(seg[1])
                    else:
                        H[seg[0].coord.tostring()] = [seg[0]] # save for future reference as current vertex
                        H[seg[0].coord.tostring()].append(seg[1])
                    if seg[1].coord.tostring() in H:
                        H[seg[1].coord.tostring()].append(seg[0])
                    else:
                        H[seg[1].coord.tostring()] = [seg[1]] # save for future reference as current vertex
                        H[seg[1].coord.tostring()].append(seg[0])

            while len(H) > 0:
                # Init polygon.
                polygon = gt.polygon3f()

                # Choose first returned key.
                u_key    = list(H.keys())[0]
                u_vertex = H[u_key][0]

                # Set final and next (v) keys and vertices.
                final_key    = H[u_key][2].coord.tostring() 
                final_vertex = H[u_key][2] 
                v_key        = H[u_key][1].coord.tostring()
                v_vertex     = H[u_key][1] 

                # Not needed anymore.
                H.pop(u_key)

                # Insert vertex to polygon.
                polygon.vertices.append(u_vertex)

                while u_key != final_key:
                    
                    # Get next.
                    last_key    = u_key
                    last_vertex = u_vertex

                    u_key    = v_key
                    u_vertex = H[u_key][0]

                    v_key   = H[u_key][1].coord.tostring()
                    if v_key != last_key:
                        v_vertex = H[u_key][1] 
                    else:
                        v_key    = H[u_key][2].coord.tostring()
                        v_vertex = H[u_key][2] 
                        
                    H.pop(u_key)
                
                    # Insert vertex to polygon.
                    polygon.vertices.append(u_vertex)

                # Insert polygon to plane.
                self.planes[i].append(polygon)

        # Compute proper orientation for polygons.
        for plane in self.planes:
            for p1 in range(len(plane)):
                orientation = gt.COUNTERCLOCKWISE
                for p2 in range(len(plane)):
                    if p1 != p2:
                        if plane[p2].is_inside(plane[p1].vertices[0]):
                            orientation = gt.CLOCKWISE
                
                if orientation != plane[p1].orientation():
                    plane[p1].invert_orientation()


    ## OpenGL data for planes.
    #
    # Create data array to draw planes in OpenGL.
    def OpenGLPlanesData(self,mmin,mmax,bmin,bmax):
        
        # Normalize planes.
        planes = (self.P-mmin)/(mmax-mmin)*(bmax-bmin)+bmin

        step = 2
        # Lines * columns * 2 triangles * planes 
        ntriangles = (2.0/step)*(2.0/step)*2.0*len(planes)
        nvertices  = int(ntriangles*3.0)
        # Data elements = nvertices*(3 coordinates + 3 normal coordinates)
        data = np.zeros(nvertices*6, dtype='float32')

        l = 0
        for p in planes:
            for j in range(-1,1,step):
                for k in range(-1,1,step):
                    data[l:l+6] = [j, k, p, 0, 1, 0]
                    l += 6
                    data[l:l+6] = [j, k+step, p, 0, 1, 0]
                    l += 6
                    data[l:l+6] = [j+step, k, p, 0, 1, 0]
                    l += 6
                    data[l:l+6] = [j, k+step, p, 0, 1, 0]
                    l += 6
                    data[l:l+6] = [j+step, k+step, p, 0, 1, 0]
                    l += 6
                    data[l:l+6] = [j+step, k, p, 0, 1, 0]
                    l += 6

        return data,nvertices

    ## OpenGL data for polygons.
    #
    # Create data array to draw polygons in OpenGL.
    def OpenGLPolygonsData(self,mmin,mmax,bmin,bmax):

        # Number of segments
        nsegments = 0
        for plane in self.planes:
            for polygon in plane:
                nsegments += 2*len(polygon.vertices)
       
        # Number of vertices
        nvertices = nsegments*2
        
        # Create data array of size nsegments*2 vertices*(3 coordinates + 3 colors)
        data = np.zeros(nsegments*2*6, dtype='float32')

        # Copy segments to data array.
        j = 0
        for plane in self.planes:
            for polygon in plane:

                if polygon.orientation() == gt.COUNTERCLOCKWISE:
                    color = np.array([0.0, 1.0, 0.0], dtype='float32')
                else:
                    color = np.array([0.0, 0.0, 1.0], dtype='float32')

                for v in range(len(polygon.vertices)):
                    data[j:j+3] = (polygon.vertices[v].coord-mmin)/(mmax-mmin)*(bmax-bmin)+bmin
                    j += 3
                    data[j:j+3] = color
                    j += 3
                    data[j:j+3] = (polygon.vertices[(v+1)%len(polygon.vertices)].coord-mmin)/(mmax-mmin)*(bmax-bmin)+bmin
                    j += 3
                    data[j:j+3] = color
                    j += 3

        return data,nvertices

       
    ## OpenGL data.
    #
    # Create data array to draw segments in OpenGL.
    def OpenGLData(self,mmin,mmax,bmin,bmax):
       
        # Number of segments
        nsegments = 0
        for l in self.segments:
            nsegments += len(l)

        # Number of vertices
        nvertices = nsegments*2
        
        # Create data array of size nsegments*2 vertices*3coordinates
        data = np.zeros(nsegments*2*3, dtype='float32')

        # Copy segments to data array.
        j = 0
        for l in self.segments:
            for s in l:
                data[j:j+3] = (s[0:3]-mmin)/(mmax-mmin)*(bmax-bmin)+bmin
                j += 3
                data[j:j+3] = (s[3:6]-mmin)/(mmax-mmin)*(bmax-bmin)+bmin
                j += 3

        return data,nvertices

    ## Round
    #
    #
    def mround(self, v, eps, mod, rem):
        y = np.round(v/(mod*eps))
        z = (y*mod+rem)*eps;
        return z

    ## Round coordinates of triangles.
    #
    # round the coordinates of triangles to even multiples of eps.
    def roundTrianglesEven(self):
        
        for t in self.T:
            t.v0 = self.mround(t.v0, self.eps, 2, 0)
            t.v1 = self.mround(t.v1, self.eps, 2, 0)
            t.v2 = self.mround(t.v2, self.eps, 2, 0)
            
            
    ## Return a vector of minimum values for coordinates in the mesh.
    def min_coordinates(self):
        
        Mval = np.ones(3, dtype='float32')*np.finfo(np.float32).max
        for t in self.T: 
            Mval = np.minimum(Mval,t.v0)
            Mval = np.minimum(Mval,t.v1)
            Mval = np.minimum(Mval,t.v2)

        return Mval
    
    ## Return a vector of maximum values for coordinates in the mesh.
    def max_coordinates(self):
        
        Mval = np.ones(3, dtype='float32')*np.finfo(np.float32).min
        for t in self.T: 
            Mval = np.maximum(Mval,t.v0)
            Mval = np.maximum(Mval,t.v1)
            Mval = np.maximum(Mval,t.v2)

        return Mval

    ## Create uniform planes.
    #
    # Create planes at odd multiples of delta.
    def uniformPlanes(self):

        z_min = self.min_coordinates()[2]
        z_max = self.max_coordinates()[2]
        
        # First plane.
        P0 = self.mround(z_min - self.delta, self.eps, 2, 1)

        # Compute planes.
        self.P = []
        i = 1
        Pi = P0 + i * self.delta
        while Pi < z_max:
            if Pi > z_min:
                self.P.append(Pi)
            i += 1
            Pi = P0 + i * self.delta
