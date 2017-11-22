#!/bin/bash
# Usage: ./script.sh PHI RHO 30
#                    ini end cpu
#str1=$1
#str2=$2
nCPUs=32
#for structure in $(ls ZIF_structures_Compendium/*.cif | sed 's/_ZIF_TOBUNPOROUS_Optimized\.cif//g' | sed 's/ZIF_structures_Compendium\///g' | sed -n "/${str1}/,/${str2}/p") ; do
for CIFFile in $(ls ZIF_structures_Compendium/*.cif) ; do
  structure=$(echo $CIFFile | sed 's/\.cif//g' | sed 's/ZIF_structures_Compendium\///g')
  file=$(echo $CIFFile | sed 's/ZIF_structures_Compendium\///g')
  cp ZIF_structures_Compendium/$file .
  bash main.sh $file $nCPUs > salida.${structure}.txt 
done
exit 0
