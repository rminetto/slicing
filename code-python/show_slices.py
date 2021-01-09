#!/usr/bin/python3

import sys
import os
import argparse
from ctypes import c_void_p
import numpy as np
import math
import OpenGL.GL as gl
import OpenGL.GLUT as glut
import glutils as ut
import stlmesh
import slicer

## Screen width.
win_width  = 800
## Screen height.
win_height = 600

## Shaders for drawing the mesh and the planes.
program      = None
## Shaders for drawing the contours.
program2     = None
## Vertex array for the mesh.
VAO_mesh     = None
## Vertex buffer for the mesh.
VBO_mesh     = None
## Vertex array for the planes.
VAO_planes   = None
## Vertex buffer for the planes.
VBO_planes   = None
## Vertex array for the segments/contours.
VAO_segments = None
## Vertex buffer for the segments/contours.
VBO_segments = None

## Auxiliar to track mouse x movement.
lastX = win_width/2.0
## Auxiliar to track mouse y movement.
lastY = win_height/2.0
## To check when mouse is clicked.
firstMouse = True

## Set camera position.
camPosition = np.array([0.0, 0.0, -5.0], dtype='float32')
## Set camera target position.
camTarget = np.array([0.0, 0.0, 0.0], dtype='float32')
## Set camera up.
camUp = np.array([0.0, 1.0, 0.0], dtype='float32')
## Set field of view.
fov = 45.0

## Set light color.
lightColor    = np.array([1.0, 1.0, 1.0], dtype='float32')
## Set light position.
lightPosition = np.array([0.0, 0.0, -10.0], dtype='float32')

## STL file name.
stl_file = None
## Number of vertices to be drawn in mesh.
num_vertices_mesh = 0
## Set mesh color.
meshColor = np.array([0.8, 0.8, 0.8, 0.3], dtype='float32')

## Displacement for planes.
delta = None

## Number of vertices to be drawn in planes.
num_vertices_planes = 0
## Set planes color.
planeColor = np.array([1.0, 1.0, 0.0, 0.1], dtype='float32')

## Number of vertices to be drawn in segments/contours.
num_vertices_segments = 0

## Flag to draw mesh.
show_mesh     = True
## Flag to draw planes.
show_planes   = False
## Flag to draw segments.
show_segments = True
## Normalize object to view cube.
view_max =  1.0
## Normalize object to view cube.
view_min = -1.0
## Set model matrix.
model = ut.matRotateX(math.radians(-90.0))

## Mesh and planes vertex shader.
mesh_vs = """
#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 vNormal;
out vec3 fragPosition;

void main()
{
    gl_Position = projection * view * model * vec4(position, 1.0);
    fragPosition = vec3(model * vec4(position, 1.0));
    vNormal = normal;
}
"""

## Mesh and planes fragment shader.
mesh_fs = """
#version 330 core
in vec3 vNormal;
in vec3 fragPosition;

out vec4 fragColor;

uniform vec4 objectColor;
uniform vec3 lightColor;
uniform vec3 lightPosition;
uniform vec3 cameraPosition;

void main()
{
    float ambientStrength = 0.1;
    vec3 ambient = ambientStrength * lightColor;

    vec3 norm = normalize(vNormal);
    vec3 lightDirection = normalize(lightPosition - fragPosition);
    float diff = max(dot(norm, lightDirection), 0.0);
    vec3 diffuse = diff * lightColor;

    float specularStrength = 0.5;
    vec3 cameraDirection = normalize(cameraPosition - fragPosition);
    vec3 reflectDirection = reflect(-lightDirection, norm);
    float spec = pow(max(dot(cameraDirection, reflectDirection), 0.0), 1);
    vec3 specular = specularStrength * spec * lightColor;

    vec4 result = (vec4(ambient,1.0) + vec4(diffuse,1.0) + vec4(specular,1.0)) * objectColor;
    fragColor = result;
} 
"""


## Segments vertex shader.
segment_vs = """
#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec3 color;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec4 vertexColor;

void main()
{
    gl_Position = projection * view * model * vec4(position, 1.0);
    vertexColor = vec4(color, 1.0);
}
"""

## Segments fragment shader.
segment_fs = """
#version 330 core

in vec4 vertexColor;
out vec4 fragColor;

void main()
{
    fragColor = vertexColor;
} 
"""

## Display function.
#
# Draws mesh, planes and segments.
def display():

    gl.glClearColor(0.0, 0.0, 0.0, 1.0)
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    
    if show_segments: draw_segments()
    if show_mesh: draw_mesh()
    if show_planes: draw_planes()
        
    glut.glutSwapBuffers()

## Draw segments.
#
# Set how to draw segments.
def draw_segments():
    # Enable shader program for segments
    gl.glUseProgram(program2)
    
    # Set view matrix
    view = ut.matLookAt(camPosition[0], camPosition[1], camPosition[2], camTarget[0], camTarget[1], camTarget[2], camUp[0], camUp[1], camUp[2])
    loc = gl.glGetUniformLocation(program2, "view")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, view.transpose())

    # Set projection matrix
    projection = ut.matPerspective(np.radians(fov), win_width/win_height, 0.1, 100.0)
    loc = gl.glGetUniformLocation(program2, "projection")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, projection.transpose())

    # Set current model matrix (global)
    loc = gl.glGetUniformLocation(program2, "model")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, model.transpose())

    # Draw segments
    gl.glLineWidth(1.5)
    gl.glBindVertexArray(VAO_segments)
    gl.glDrawArrays(gl.GL_LINES, 0, num_vertices_segments)
    gl.glBindVertexArray(0)
    gl.glLineWidth(1.0)

## Draw mesh.
#
# Set how to draw mesh.
def draw_mesh():
    # Enable shader program
    gl.glUseProgram(program)

    # Set view matrix
    view = ut.matLookAt(camPosition[0], camPosition[1], camPosition[2], camTarget[0], camTarget[1], camTarget[2], camUp[0], camUp[1], camUp[2])
    loc = gl.glGetUniformLocation(program, "view")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, view.transpose())

    # Set projection matrix
    projection = ut.matPerspective(np.radians(fov), win_width/win_height, 0.1, 100.0)
    loc = gl.glGetUniformLocation(program, "projection")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, projection.transpose())

    # Set current model matrix (global)
    loc = gl.glGetUniformLocation(program, "model")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, model.transpose())

    # Set light attributes
    loc = gl.glGetUniformLocation(program, "lightColor")
    gl.glUniform3f(loc, lightColor[0], lightColor[1], lightColor[2])
    loc = gl.glGetUniformLocation(program, "lightPosition")
    gl.glUniform3f(loc, lightPosition[0], lightPosition[1], lightPosition[2])
    
    # Set color attributes for mesh
    loc = gl.glGetUniformLocation(program, "objectColor")
    gl.glUniform4f(loc, meshColor[0], meshColor[1], meshColor[2], meshColor[3])
   
    # Draw mesh
    gl.glBindVertexArray(VAO_mesh)
    gl.glDrawArrays(gl.GL_TRIANGLES, 0, num_vertices_mesh)
    gl.glBindVertexArray(0)

## Draw planes.
#
# Set how to draw planes.
def draw_planes():
    # Enable shader program
    gl.glUseProgram(program)
    
    # Set view matrix
    view = ut.matLookAt(camPosition[0], camPosition[1], camPosition[2], camTarget[0], camTarget[1], camTarget[2], camUp[0], camUp[1], camUp[2])
    loc = gl.glGetUniformLocation(program, "view")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, view.transpose())

    # Set projection matrix
    projection = ut.matPerspective(np.radians(fov), win_width/win_height, 0.1, 100.0)
    loc = gl.glGetUniformLocation(program, "projection")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, projection.transpose())

    # Set current model matrix (global)
    loc = gl.glGetUniformLocation(program, "model")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, model.transpose())

    # Set light attributes
    loc = gl.glGetUniformLocation(program, "lightColor")
    gl.glUniform3f(loc, lightColor[0], lightColor[1], lightColor[2])
    loc = gl.glGetUniformLocation(program, "lightPosition")
    gl.glUniform3f(loc, lightPosition[0], lightPosition[1], lightPosition[2])

    # Set color attributes for mesh
    loc = gl.glGetUniformLocation(program, "objectColor")
    gl.glUniform4f(loc, planeColor[0], planeColor[1], planeColor[2], planeColor[3])
    
    # Set current model matrix (global)
    S = ut.matScale(1.5,1.5,1.0)
    loc = gl.glGetUniformLocation(program, "model")
    gl.glUniformMatrix4fv(loc, 1, gl.GL_FALSE, np.matmul(model,S).transpose())

    # Draw planes
    gl.glBindVertexArray(VAO_planes)
    gl.glDrawArrays(gl.GL_TRIANGLES, 0, num_vertices_planes)
    gl.glBindVertexArray(0)

## Reshape function.
#
# Set how to reshape screen.
def reshape(width,height):
    win_width  = width
    win_height = height
    gl.glViewport(0, 0, width, height)

## Keyboard callback.
#
# Set actions for keyboard keys.
def keyboard( key, x, y ):
    global show_mesh    
    global show_planes  
    global show_segments

    if key == b'\x1b'or key == b'q':
        sys.exit( )
    if key == b'm':
        show_mesh = not(show_mesh)
    if key == b'p':
        show_planes = not(show_planes)
    if key == b's':
        show_segments = not(show_segments)

    glut.glutPostRedisplay()


## Mouse function.
#
# Set actions for mouse events.
def mouse(btn,state,x,y):
    global fov
    global firstMouse

    if btn == glut.GLUT_LEFT_BUTTON and state == glut.GLUT_UP:
        firstMouse = True
    if btn == 3:
        camPosition[2] += 0.1
    if btn == 4:
        camPosition[2] -= 0.1
    
    glut.glutPostRedisplay()

## Movement function.
#
# Set actions for mouse movement.
def motion(x,y):
    global lastX
    global lastY
    global firstMouse
    global model

    if firstMouse:
        lastX = x
        lastY = y
        firstMouse = False
    else:
        xoffset = x - lastX
        yoffset = lastY - y 
        lastX = x
        lastY = y

        angle_inc = 0.5
        if xoffset > 0:
            angleY = xoffset*angle_inc
            Ry = ut.matRotateY(math.radians(angleY))
            model = np.matmul(Ry,model)
        if xoffset < 0:
            angleY = xoffset*angle_inc
            Ry = ut.matRotateY(math.radians(angleY))
            model = np.matmul(Ry,model)
        if yoffset > 0:
            angleX = yoffset*angle_inc
            Rx = ut.matRotateX(math.radians(angleX))
            model = np.matmul(Rx,model)
        if yoffset < 0:
            angleX = yoffset*angle_inc
            Rx = ut.matRotateX(math.radians(angleX))
            model = np.matmul(Rx,model)

    glut.glutPostRedisplay()

## Init program.
#
# All initial stuff, reading file, reading shaders, creating buffers.
def init():
    global program
    global program2
    global VAO_mesh
    global VBO_mesh
    global VAO_planes
    global VBO_planes
    global VAO_segments
    global VBO_segments
    global num_vertices_mesh
    global num_vertices_planes
    global num_vertices_segments

    # Build programs (shaders).
    program  = ut.createShaderProgram(mesh_vs, mesh_fs)
    program2 = ut.createShaderProgram(segment_vs, segment_fs)

    # Load mesh.
    mesh = stlmesh.stlmesh(stl_file)
    mesh_min = mesh.min_coordinates()[2]
    mesh_max = mesh.max_coordinates()[2]

    # Compute slices.
    P = None
    srt   = False
    mesh_slicer = slicer.slicer(mesh.triangles,P,delta,srt)
    mesh_slicer.incremental_slicing()

    # Create mesh data for GPU.
    data,num_vertices_mesh = mesh.OpenGLData(view_min, view_max)

    # Copy mesh data to GPU
    VAO_mesh = gl.glGenVertexArrays(1)
    VBO_mesh = gl.glGenBuffers(1)
    
    gl.glBindVertexArray(VAO_mesh)
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, VBO_mesh)
    gl.glBufferData(gl.GL_ARRAY_BUFFER, data.nbytes, data, gl.GL_STATIC_DRAW)
    gl.glVertexAttribPointer(0, 3, gl.GL_FLOAT, gl.GL_FALSE, 6*data.itemsize, None)
    gl.glEnableVertexAttribArray(0)
    gl.glVertexAttribPointer(1, 3, gl.GL_FLOAT, gl.GL_FALSE, 6*data.itemsize, c_void_p(3*data.itemsize))
    gl.glEnableVertexAttribArray(1)
    
    # Unbind
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
    gl.glBindVertexArray(0)

    # Create data for planes.
    data,num_vertices_planes = mesh_slicer.OpenGLPlanesData(mesh_min,mesh_max,view_min,view_max)
    
    # Copy planes data to GPU
    VAO_planes = gl.glGenVertexArrays(1)
    VBO_planes = gl.glGenBuffers(1)
    
    gl.glBindVertexArray(VAO_planes)
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, VBO_planes)
    gl.glBufferData(gl.GL_ARRAY_BUFFER, data.nbytes, data, gl.GL_STATIC_DRAW)
    gl.glVertexAttribPointer(0, 3, gl.GL_FLOAT, gl.GL_FALSE, 6*data.itemsize, None)
    gl.glEnableVertexAttribArray(0)
    gl.glVertexAttribPointer(1, 3, gl.GL_FLOAT, gl.GL_FALSE, 6*data.itemsize, c_void_p(3*data.itemsize))
    gl.glEnableVertexAttribArray(1)
    
    # Unbind
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
    gl.glBindVertexArray(0)

    # Create slices data for GPU.
    mesh_max = mesh.max_coordinates()
    mesh_min = mesh.min_coordinates()
    data,num_vertices_segments = mesh_slicer.OpenGLPolygonsData(mesh_min,mesh_max,view_min,view_max)
    
    # Copy segments data to GPU
    VAO_segments = gl.glGenVertexArrays(1)
    VBO_segments = gl.glGenBuffers(1)
    
    gl.glBindVertexArray(VAO_segments)
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, VBO_segments)
    gl.glBufferData(gl.GL_ARRAY_BUFFER, data.nbytes, data, gl.GL_STATIC_DRAW)
    gl.glVertexAttribPointer(0, 3, gl.GL_FLOAT, gl.GL_FALSE, 6*data.itemsize, None)
    gl.glEnableVertexAttribArray(0)
    gl.glVertexAttribPointer(1, 3, gl.GL_FLOAT, gl.GL_FALSE, 6*data.itemsize, c_void_p(3*data.itemsize))
    gl.glEnableVertexAttribArray(1)
    
    # Unbind
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
    gl.glBindVertexArray(0)

    # Set depth 
    gl.glEnable(gl.GL_DEPTH_TEST)

    gl.glEnable(gl.GL_BLEND)
    gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)


## Parsing.
#
# Parse input arguments.
def parse_input():
    global stl_file
    global delta

    parser = argparse.ArgumentParser(description='Visualize slicer. User can toggle visualized objects using keyboard keys: m for mesh, p for planes and s for segments.')
    parser.add_argument(dest='stl', metavar='STL', help='STL file.')
    parser.add_argument('-d', dest='delta', metavar='delta', type=float, default=2.0, help='Define the displacement between planes.')
    args = parser.parse_args()

    if not(os.path.isfile(args.stl)):
        print ("File " + args.stl + " does not exist.")
        sys.exit()

    stl_file = args.stl
    delta = args.delta

## Main.
#
# Main function.
def main():

    glut.glutInit()
    
    parse_input()

    glut.glutInitDisplayMode(glut.GLUT_DOUBLE | glut.GLUT_RGBA | glut.GLUT_DEPTH)
    glut.glutInitWindowSize(win_width,win_height)
    glut.glutInitContextVersion(3, 3);
    glut.glutInitContextProfile(glut.GLUT_CORE_PROFILE);
    glut.glutCreateWindow('Mesh slices')
    
    init()

    glut.glutReshapeFunc(reshape)
    glut.glutDisplayFunc(display)
    glut.glutKeyboardFunc(keyboard)
    glut.glutMouseFunc(mouse)
    glut.glutMotionFunc(motion)

    glut.glutMainLoop()

if __name__ == '__main__':
    main()
