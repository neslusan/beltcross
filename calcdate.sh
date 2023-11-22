#!/bin/bash
#
#$ -cwd
#

echo "Calculation of MOID... (wait for a while until done, please)";
./moid.exes
echo "Calculation of relative velocity and angle between";
echo "  the velocity vectors...";
./relvel.exes
echo "Calculation of basic date of passage...";
./basicdate.exes
echo "Creation of the list of passages in year under interest...";
./datelist.exes
echo "Arranging the passages in order of increasing date...";
./arrange.exes

jyear=$(< year.d);

echo " ";
echo "For the input parameters given in file <object.d>,";
echo "the prediction of the passages of object through the meteoroid";
echo "stream in year "`printf "%04d" $jyear`" is written in file <datelist"`printf "%04d" $jyear`".dat>.";
echo " ";
echo "Program successfully terminated.";
