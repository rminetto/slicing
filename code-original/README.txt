To compile:

g++ -O3 -o slicing slicing.cpp -frounding-math

To run:

./slicing -Incremental -model ./stl_models/02.femur.stl -thickness 2.0 -adaptive false -2D -rotate false -orienting_contours true

where

the options for slicing are: -Trivial, -Park, -Incremental.

the options for adaptive are: true or false.

the options for the output are: -2D, -3D, video, No.

the options for rotation are: true or false.

the options for contour orientation (clockwise or counter-clockwise) are: true or false.

The thickness is given in millimeters.

Author: Rodrigo Minetto.
