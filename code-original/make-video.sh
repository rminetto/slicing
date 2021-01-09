#! /bin/bash

video=video;

for svg in $(ls ${video}/*.svg); do
   frame=`echo ${svg} | cut -f1 -d"."`;
   inkscape ${frame}.svg --export-png=${frame}.png;
   rm ${frame}.svg;
   echo $frame;
done

mencoder mf://${video}/?????_?????.png -mf w=1024:h=768:fps=90:type=png -ovc lavc -lavcopts vcodec=msmpeg4v2 -nosound -o fatiamento.mpg;
