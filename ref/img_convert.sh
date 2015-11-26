#!/bin/bash

cd huntington_1280/

for i in $( ls *.pgm);
do 
  convert -resize 640x360 $i ../huntington_640/$i;
done
