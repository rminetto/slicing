## @file Basic functions for OpenGL.
#
# Defines a set of util functions for OpenGL.

import os
import sys
import math
import numpy as np
import OpenGL.GL as gl

## Read shader.
#
# Read a shader from a file.
#
# @param File name.
# @return String with shader code.
def readShaderFile(shader_file):

    # Check if file exists.
    if not os.path.isfile(shader_file):
        print("Could not open shader file: " + shader_file)
        sys.exit()

    # Read file.
    shader_code = None
    with open(shader_file, 'r') as f:
        shader_code = f.read()

    return shader_code


## Create program.
#
# Creates a program from given shader codes.
#
# @param vertex_code String with code for vertex shader.
# @param fragment_code String with code for fragment shader.
# @return Compiled program.
def createShaderProgram(vertex_code, fragment_code):
    
    # Request a program and shader slots from GPU
    program  = gl.glCreateProgram()
    vertex   = gl.glCreateShader(gl.GL_VERTEX_SHADER)
    fragment = gl.glCreateShader(gl.GL_FRAGMENT_SHADER)
    
    # Set shaders source
    gl.glShaderSource(vertex, vertex_code)
    gl.glShaderSource(fragment, fragment_code)

    # Compile shaders
    gl.glCompileShader(vertex)
    if not gl.glGetShaderiv(vertex, gl.GL_COMPILE_STATUS):
        error = gl.glGetShaderInfoLog(vertex).decode()
        print(error)
        raise RuntimeError("Shader compilation error")
                
    gl.glCompileShader(fragment)
    if not gl.glGetShaderiv(fragment, gl.GL_COMPILE_STATUS):
        error = gl.glGetShaderInfoLog(fragment).decode()
        print(error)
        raise RuntimeError("Shader compilation error")                

    # Attach shader objects to the program
    gl.glAttachShader(program, vertex)
    gl.glAttachShader(program, fragment)

    # Build program
    gl.glLinkProgram(program)
    if not gl.glGetProgramiv(program, gl.GL_LINK_STATUS):
        print(gl.glGetProgramInfoLog(program))
        raise RuntimeError('Linking error')

    # Get rid of shaders (no more needed)
    gl.glDetachShader(program, vertex)
    gl.glDetachShader(program, fragment)

    return program

## Zeros matrix..
#
# Create a 4x4 matrix with all elements set to zero.
#
# @return The matrix.
def matZeros():
    Z = np.zeros((4,4), dtype='float32')
    
    return Z


## Identity matrix.
#
# Create a 4x4 identity matrix.
#
# @param .
# @return .
def matIdentity():
    I = np.identity(4, dtype='float32')
    
    return I

## Translate.
#
# Create a 4x4 translation matrix.
#
# @param x Displacement in the x-axis.
# @param y Displacement in the y-axis.
# @param z Displacement in the z-axis.
# @return Translation matrix.
def matTranslate(x, y, z):
    T =  np.identity(4, dtype='float32')

    T[0,3] = x
    T[1,3] = y
    T[2,3] = z

    return T

## Scale.
#
# Create a 4x4 scale matrix.
#
# @param x Scale factor for x-axis.
# @param y Scale factor for y-axis.
# @param z Scale factor for z-axis.
# @return Scale matrix.
def matScale(x, y, z):
    T =  np.identity(4, dtype='float32')

    T[0,0] = x
    T[1,1] = y
    T[2,2] = z

    return T

## X rotation.
#
# Create a 4x4 rotation matrix for x-axis.
#
# @param angle Angle in radians.
# @return Rotation matrix.
def matRotateX(angle):
    R =  np.identity(4, dtype='float32')

    acos = math.cos(angle)
    asin = math.sin(angle)

    R[1,1] =  acos
    R[1,2] = -asin
    R[2,1] =  asin
    R[2,2] =  acos

    return R

## Y rotation.
#
# Create a 4x4 rotation matrix for z-axis.
#
# @param angle Angle in radians.
# @return Rotation matrix.
def matRotateY(angle):
    R =  np.identity(4, dtype='float32')

    acos = math.cos(angle)
    asin = math.sin(angle)

    R[0,0] =  acos
    R[0,2] =  asin
    R[2,0] = -asin
    R[2,2] =  acos

    return R

## Z rotation.
#
# Create a 4x4 rotation matrix for z-axis.
#
# @param angle Angle in radians.
# @return Rotation matrix.
def matRotateZ(angle):
    R =  np.identity(4, dtype='float32')

    acos = math.cos(angle)
    asin = math.sin(angle)

    R[0,0] =  acos
    R[0,1] = -asin
    R[1,0] =  asin
    R[1,1] =  acos

    return R

## Perspective.
#
# Create a 4x4 perspective projection matrix.
#
# @param fovy Field of view.
# @param aspect Aspect value.
# @param n Near value.
# @param f Far value.
# @return Perspective matrix.
def matPerspective(fovy, aspect, n, f):

    P = np.zeros((4,4), dtype='float32')
    
    rad = fovy
    tan = math.tan(rad/2.0)

    P[0,0] = 1.0/(aspect*tan)
    P[1,1] = 1.0/tan
    P[2,2] = -(f+n)/(f-n);
    P[2,3] = - (2.0*f*n)/(f-n);
    P[3,2] = -1.0;

    return P

## Frustrum.
#
# Create a 4x4 perspective projection matrix.
#
# @param l Left value.
# @param r Right value.
# @param b Bottom value.
# @param t Top value.
# @param n Near value.
# @param f Far value.
# @return Frustrum matrix.
def matFrustum(l, r, b, t, n, f):
    
    F = np.zeros((4,4), dtype='float32')

    F[0,0] = (2.0*n)/(r-l)
    F[0,2] = (r+l)/(r-l)
    F[1,1] = (2.0*n)/(t-b)
    F[1,2] = (t+b)/(t-b)
    F[2,2] = (-(f+n))/(f-n)
    F[2,3] = (-2.0*f*n)/(f-n)
    F[3,2] = -1.0

    return F

## Ortho matrix.
#
# Create a 4x4 orthogonal projection matrix.
#
# @param l Left value.
# @param r Right value.
# @param b Bottom value.
# @param t Top value.
# @param n Near value.
# @param f Far value.
# @return Orthogonal matrix.
def matOrtho(l, r, b, t, n, f):

    F = np.zeros((4,4), dtype='float32')

    F[0,0] = 2.0/(r-l)
    F[0,3] = -(r+l)/(r-l)
    F[1,1] = 2.0/(t-b)
    F[1,3] = -(t+b)/(t-b)
    F[2,2] = -2.0/(f-n)
    F[2,3] = -(f+n)/(f-n)
    F[3,3] = 1.0

    return F

## Normalize.
#
# Normalize a vector.
#
# @param v Vector.
# @return Normalized vector.
def vecNormalize(v):

    return (v/np.linalg.norm(v))

## View matrix.
#
# Create view matrix.
#
# @param px Camera position x coordinate.
# @param py Camera position y coordinate.
# @param pz Camera position z coordinate.
# @param tx Target position x coordinate.
# @param ty Target position y coordinate.
# @param tz Target position z coordinate.
# @param ux Camera up vector x coordinate.
# @param uy Camera up vector y coordinate.
# @param uz Camera up vector z coordinate.
# @return View matrix.
def matLookAt(px,py,pz,tx,ty,tz,ux,uy,uz):
    
    cameraPosition  = np.array([px, py, pz], dtype='float32')
    cameraTarget    = np.array([tx, ty, tz], dtype='float32')
    cameraDirection = vecNormalize(cameraPosition-cameraTarget)
    up              = np.array([ux, uy, uz], dtype='float32')
    cameraRight     = vecNormalize(np.cross(up,cameraDirection))
    cameraUp        = np.cross(cameraDirection, cameraRight);

    L = np.identity(4, dtype='float32')

    R = np.identity(4, dtype='float32')
    R[0,0] = cameraRight[0]
    R[0,1] = cameraRight[1]
    R[0,2] = cameraRight[2]
    R[1,0] = cameraUp[0]
    R[1,1] = cameraUp[1]
    R[1,2] = cameraUp[2]
    R[2,0] = cameraDirection[0]
    R[2,1] = cameraDirection[1]
    R[2,2] = cameraDirection[2]
    
    T = np.identity(4, dtype='float32')
    T[0,3] = -cameraPosition[0]
    T[1,3] = -cameraPosition[1]
    T[2,3] = -cameraPosition[2]

    L = np.matmul(R,T)

    return L
