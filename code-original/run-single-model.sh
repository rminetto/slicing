#! /bin/bash

video=video;

rm *.svg;

mkdir ${video};

rm ${video}/*;

g++ -O3 -o slicing slicing.cpp -frounding-math;

#The options for the output are: [-2D, -3D, -video, -No]
#The options for slicing are: [-Trivial, -Park, -Incremental]
#The options for adaptive are: [true, false]. 
#The thickness is given in millimeters. 
#The options for rotation are: [true, false] (to change de Z cutting plane).
#The options for contour orientation (clockwise or counter-clockwise) are: [true, false].

./slicing -Incremental -model ./stl_models/02.femur.stl -thickness 2.0 -adaptive false -3D -rotate false -orienting_contours true

#./make-video.sh;
