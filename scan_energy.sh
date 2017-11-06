#!/bin/bash
#cat *_ABW_*_emmd/logs/minimization_postMD.txt
for structure in GIS ; do
 energy0=$(tail -n1 00_${structure}_*_emmd/logs/minimization_postMD.txt | awk '{print $5}')
 for molecule in empty Ar Xe ; do
  for file in *_${structure}_*_${molecule}_*_emmd/logs/minimization_postMD.txt ; do
   n_molecules=$(echo $file | sed 's/_/ /g' | awk '{print $7}')
   step=$(tail -n1 $file | awk '{print $2}')
   energy=$(tail -n1 $file | awk '{print $5}')
   echo $step $energy0 $energy $(echo "scale=4; ($energy - $energy0)/(${n_molecules})" | bc -l) $n_molecules $file 
  done
 done
done
