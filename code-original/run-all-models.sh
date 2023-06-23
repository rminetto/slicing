#! /bin/bash

#The options for the output are: [-2D, -3D, -video, -No]
#The options for slicing are: [-Trivial, -Park, -Incremental]
#The options for adaptive are: [true, false]. 
#The thickness is given in millimeters. 
#The options for rotation are: [true, false] (to change de Z cutting plane).
#The options for contour orientation (clockwise or counter-clockwise) are: [true, false].

g++ -O3 -o slicing slicing.cpp -frounding-math;

foutput=RUNNING_TIME.txt;

for model in $(ls ./stl_models/); do
   echo $model;
   time ./slicing -Trivial -model ./stl_models/${model} -thickness 0.032 -adaptive false -no -rotate false -orienting_contours false >> ${foutput};
   time ./slicing -Incremental -model ./stl_models/${model} -thickness 0.032 -adaptive false -no -rotate false -orienting_contours false >> ${foutput};
   time ./slicing -Park -model ./stl_models/${model} -thickness 0.032 -adaptive false -no -rotate false -orienting_contours false >> ${foutput};
   echo "-------------------------------------" >> ${foutput};
done
