#!/bin/bash

rootname=$1
rootname="${rootname:1}"

if grep -q "rs$rootname" "i$rootname"; then
  echo "Starting from restart file"
else
  folder_name="run_"$rootname
  rm -rf $folder_name
  mkdir $folder_name
  cd $folder_name
  cp ../i$rootname .
fi

mpirun -n 9 ~/HMHD2D/build/hmhd2d $rootname
