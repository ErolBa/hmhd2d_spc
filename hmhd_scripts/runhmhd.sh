#!/bin/bash

folder_name="run_"$1
rm -rf $folder_name
mkdir $folder_name
cd $folder_name
cp ../i$1 .

mpirun -n 9 ~/HMHD2D/build/hmhd2d $1
