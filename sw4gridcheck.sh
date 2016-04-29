#!/bin/bash

if [ "$#" -ne 3 ];then
	echo "Wrong Input!"
	echo "Input: Vsmin(km/s) h(km) fmax(Hz)"
	exit 1
fi
args=("$@")
Vsmin=${args[0]} 
h=${args[1]}
fmax=${args[2]}
P=$(echo "$Vsmin/$h/$fmax" | bc -l)
echo "P value is: "$P
echo "Recommended: 6<=P<=10"
