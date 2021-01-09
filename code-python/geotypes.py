## @file Basic geometrical types.
#
# Defines a set of basis geometrical types and operations.

import numpy as np

## Constant for X position.
X = 0
## Constant for Y position.
Y = 1
## Constant for Z position.
Z = 2

## Constant for polygon orientatio.
CLOCKWISE        = -1
## Constant for polygon orientatio.
COUNTERCLOCKWISE =  1

## Constant for point in polygon test.
INSIDE_POLYGON  = 1
## Constant for point in polygon test.
OUTSIDE_POLYGON = 0


## 3D vertex.
#
# A vertex of three coordinates (x,y,z).
class vertex3f:
    
    ## Constructor.
    #
    # Create an object with a vector of coordinates (x,y,z).
    #
    # @param x The x coordinate.
    # @param y The y coordinate.
    # @param z The z coordinate.
    def __init__(self, x, y, z):
        ## Array of coordinates.
        self.coord = np.array([x, y, z], dtype='float32')

    ## Subtract.
    #
    # Return vector from subtracting two vertices.
    #
    # @param Subtract vertices.
    def __sub__(self, v):

        # Subtract coordinates.
        coord = self.coord - v.coord

        return vertex3f(coord[X], coord[Y], coord[Z])
    
    ## Normalize.
    #
    # Normalize a vector.
    def normalize(self):
        self.coord = self.coord/np.linalg.norm(self.coord)

        
## 3D polygon.
#
# A polygon represented as a sequence of vertices.
class polygon3f:
    
    ## Constructor.
    #
    # Create an object with an empty polygon.
    def __init__(self):
        ## List of vertices.
        self.vertices = []

    ## Signed area.
    #
    # Return the signed area \f$ A \f$ of the polygon: \f$ A < 0 \f$ if
    # clockwise oriented polygon or \f$ A > 0 \f$ if counterclockwise oriented polygon.
    def signed_area(self):

        # Fix first vertex as anchor vertex.
        v = self.vertices[0]

        # Compute the area of all triangles fixed with v.
        A = 0
        for i in range(1,len(self.vertices)-1):
            A += triangle_signed_area(v, self.vertices[i], self.vertices[i+1])

        return A


    ## Orientation.
    #
    # Return the orientation (CLOCKWISE or COUNTERCLOCKWISE) of the polygon.
    def orientation(self):
        
        return COUNTERCLOCKWISE if self.signed_area() > 0.0 else CLOCKWISE
    
    ## Invert orientation.
    #
    # Invert the vertex sequence of the polygon.
    def invert_orientation(self):

        self.vertices = np.flip(self.vertices)

    
    ## Point in polygon.
    #
    # Check whether point p lies inside or outside the polygon.
    # @param p Point.
    # @return INSIDE_POLYGON or OUTSIDE_POLYGONS constants.
    def is_inside(self, p):
        
        winding_number = 0
        for i in range(0,len(self.vertices)):

            winding_number += vectors_angle(p, self.vertices[i], self.vertices[(i+1)%len(self.vertices)])

        return (INSIDE_POLYGON if np.fabs(winding_number-2.0*np.pi) < np.fabs(winding_number) else OUTSIDE_POLYGON)

## 3D line segment.
#
# A line segment defined by two vertices. 
class segment3f:
    
    ## Constructor.
    #
    # Create an object with a vector of coordinates (x,y,z).
    def __init__(self, x, y, z):
        ## Array of coordinates.
        self.coord = np.array([x, y, z], dtype='float32')

## Triangle signed area.
#
# Return the signed area of a triangle.
# @param v0 First vertex of the triangle.
# @param v1 Second vertex of the triangle.
# @param v2 Third vertex of the triangle.
# @return Signed area.
def triangle_signed_area(v0, v1, v2):

    x0 = v0.coord[X]
    y0 = v0.coord[Y]
    x1 = v1.coord[X]
    y1 = v1.coord[Y]
    x2 = v2.coord[X]
    y2 = v2.coord[Y]

    return (((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0))/2.0)


## Angle between vertices.
#
# Return the angle of v1ov2.
# @param o Firts vertex (at vertex).
# @param v1 Second vertex of the triangle.
# @param v2 Third vertex of the triangle.
# @return Angle.
def vectors_angle(o, v1, v2):

    # Vectors with o.
    u = v1-o
    w = v2-o

    # Normalize vectors.
    u.normalize()
    w.normalize()

    # Return angle.
    cos_theta = np.minimum(np.dot(u.coord,w.coord),1.0) 
    cos_theta = np.maximum(cos_theta,-1.0) 
    return (np.arccos(cos_theta) if triangle_signed_area(o, v1, v2) >= 0.0 else -np.arccos(cos_theta))
