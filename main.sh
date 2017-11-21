#!/bin/bash
CIFFile=$1
nCPU=$2
n_cycles=2
temperature=85.0
pressure=0.0
filling_mode="RASPA" # Rabdel_Code
CyclesEvery=5000
InitCycles=$(echo "$CyclesEvery * 0.1" | bc -l | sed 's/\./ /g' | awk '{print $1}')
MoviesEvery=$((CyclesEvery - 1))
#
structure=$(echo $CIFFile | sed 's/\.cif//g')
seed=$(od --read-bytes=3 --address-radix=n --format=uL /dev/urandom | tr --delete " ")
while [ $(echo "$seed > 900000000" | bc -l) == 1 ] || [ $(echo "$seed < 0" | bc -l) == 1 ] ; do
 seed=$(od --read-bytes=3 --address-radix=n --format=uL /dev/urandom | tr --delete " ")
 sleep 0.5
done
CIFTemporallyFile=${structure}_${seed}.cif
# Files:
loc=$(pwd)
raspa_files_folder=$loc/lib/fff_raspa
lammps_files_folder=$loc/lib/fff_lammps
lib_folder=$loc/lib
src_files_folder=$loc/src
##############################################################
# Functions:
function make_binaries {
# Make binaries
 cp ${src_files_folder}/*.f90 .
 cp ${src_files_folder}/*.c .
 cp ${src_files_folder}/Makefile .
 make install
}
function clean_binaries {
 make clean
}
function update_name {
 NewNameFile=${structure}_${guest}_${n_Ar}
 CyclesNameFile=${cycle_name}_${NewNameFile}_${seed}
}
function check_supercell {
# Check cell size for correct calculation of energies considering cutoff
 dos2unix $CIFTemporallyFile
 cutoff=10.0
 ua=1
 ub=1
 uc=1
 a_cell=$(grep "_cell_length_a" $CIFTemporallyFile | awk '{print $2}')
 b_cell=$(grep "_cell_length_b" $CIFTemporallyFile | awk '{print $2}')
 c_cell=$(grep "_cell_length_c" $CIFTemporallyFile | awk '{print $2}')
 while [ $(echo "$a_cell < 2*$cutoff" | bc -l) == 1 ] ; do
  let "ua++"
  a_cell=$(echo "$ua*$a_cell" | bc -l)
 done
 while [ $(echo "$b_cell < 2*$cutoff" | bc -l) == 1 ] ; do
  let "ub++"
  b_cell=$(echo "$ub*$b_cell" | bc -l)
 done
 while [ $(echo "$c_cell < 2*$cutoff" | bc -l) == 1 ] ; do
  let "uc++"
  c_cell=$(echo "$uc*$c_cell" | bc -l)
 done
 echo "$CIFTemporallyFile"
 echo "Make Supercell: ${ua}x${ub}x${uc} > $a_cell $b_cell $c_cell"
 echo "==========================================================="
}
function makeCIFTopology {
 ./cif2lammps -c ${CIFTemporallyFile} -wq -S
 mv ${structure}_${seed}_topol.cif ${CIFTemporallyFile}
 rm ${structure}_${seed}.gin ${structure}_${seed}.pdb ${structure}_${seed}.data
}
function mc_muVT_raspa {
 # mu V T, ensemble
 # RASPA
 cp ${raspa_files_folder}/*.def .
 cp ${raspa_files_folder}/INPUT .
 cp ../${CIFTemporallyFile} .
 cp ../cif2lammps .
 cp $lib_folder/forcefield.lib .
 sed "s/STRUCTURE/${structure}_${seed}/g" INPUT > simulation.input
 sed -i "s/RANDOMSEED/${seed}/g"         simulation.input
 sed -i "s/TEMPERATURE/50.0/g" simulation.input
 sed -i "s/PRESSURE/1.0e9/g" simulation.input
 sed -i "s/GUEST/${guest}/g" simulation.input
 sed -i "s/CYCLESEVERY/${CyclesEvery}/g" simulation.input
 sed -i "s/INITCYCLES/${InitCycles}/g" simulation.input
 sed -i "s/MOVIESEVERY/${MoviesEvery}/g" simulation.input
 check_supercell
 sed -i "s/SUPERCELL/$ua $ub $uc/g" simulation.input
 makeCIFTopology
 echo "Run RASPA"
 go_raspa
 #wait_for_raspa
 sed '/MODEL    2/,$d' Movies/System_0/Movie*allcomponents.pdb > c
 sed 's/MODEL    2/MODEL    1/g' c > input.pdb
 rm c
 rm -rf VTK Movies/System_0/Framework* Movies/System_0/Movie_*_component_*.pdb Movies/System_0/Movie_*_frameworks.pdb
 cp ../pdb2cif .
 ./pdb2cif
}
function fill_with_guest {
 update_name
 folder=${CyclesNameFile}_filling
 if [ ! -d $folder ] ; then
  mkdir $folder
  case "${filling_mode}" in
   Rabdel_Code)
    cd $folder
     cp ../${CIFTemporallyFile} .
     cp ../adsorption_fast_atom_saturation_01 .
     ./adsorption_fast_atom_saturation_01 ${CIFTemporallyFile} > out.txt
     n_Ar=$(grep 'Ar atoms' out.txt | awk '{print $3}')
     Arname=$(echo $CIFTemporallyFile | sed 's/\.cif//g')_Ar.cif
     update_name
     cp ${Arname} ../${CIFTemporallyFile}
     cp ${Arname} ../${CyclesNameFile}.cif
    cd ..
    ;;
   RASPA)
    cd $folder
     mc_muVT_raspa
     n_Ar=$(grep "Ar " p1.cif | wc -l | awk '{print $1}')
     #n_Ar=$(grep 'absolute adsorption:' Output/System_0/output_*.data | grep 'avg' | awk '{print $3}' | tail -n1 )
     #n_Ar=${n_Ar%.*}
     Arname=$(echo $CIFTemporallyFile | sed 's/\.cif//g')_Ar.cif
     update_name
     cp p1.cif ../${CIFTemporallyFile}
     cp p1.cif ../${CyclesNameFile}.cif
    cd ..
    ;;
   *) echo "Invalid option"
      exit 0
    ;;
  esac
 fi
}
function interface_adsorption_lammps {
 lammps_file_lib="in.lmp"
 cp ${lib_folder}/forcefield.lib .
 flags="-S"
 case "${flags_cif2lammps}" in
  post-loading)
   case "${filling_mode}" in
    Rabdel_Code)
    flags="-f -wq -S"
    ;;
    RASPA)
    flags="-wq -S"
    ;;
    *)
    echo "WARNING: filling_mode NOT taken"
    flags="-wq -S"
    ;;
   esac
  ;;
  initialisation)
   flags="-S"
  ;;
  post-Xe-Ar-exchange)
   flags="-wq -S"
  ;;
  *)
   echo "WARNING: flags_cif2lammps NOT taken"
   flags="-wq -S"
  ;;
 esac 
 ./cif2lammps -c ${CyclesNameFile}.cif ${flags}
}
function first_optimisation {
 update_name
 CIF111toSupercell
 cp ${CIFTemporallyFile} ${CyclesNameFile}.cif
 ./cif2lammps -c ${CyclesNameFile}.cif -R -S 
 lammps_file_lib="in.lmp.initialitation"
 em_md_lammps
}
function em_md_lammps {
 folder=${CyclesNameFile}_emmd
 if [ ! -d $folder ] ; then
  mkdir $folder
  mv ${CyclesNameFile}.data $folder/.
  mv ${CyclesNameFile}.gin $folder/.
  mv ${CyclesNameFile}.pdb $folder/.
  mv ${CyclesNameFile}.cif $folder/.
  mv ${CyclesNameFile}_topol.cif $folder/.
  cp ${lammps_files_folder}/${lammps_file_lib} $folder/in.lmp
  mv atom_types_for_dump.txt $folder/.
  cd $folder
   sed -i "s/TEMPERATURE/${temperature}/g" in.lmp
   sed -i "s/RANDOMSEED/${seed}/g"   in.lmp
   sed -i "s/PRESSURE/${pressure}/g" in.lmp
   sed -i "s/FILENAME/${CyclesNameFile}/g" in.lmp
   elements=$(cat atom_types_for_dump.txt | sed 's/[0-9]//g' | sed 's/  / /g')
   sed -i "s/ELEMENTS/${elements}/g" in.lmp
   go_lammps
   lammps_raspa
  cd ..
 fi
}
function go_raspa {
 count_used_CPUs
 while [ $(echo "${n_used} >= ${nCPU}" | bc -l) == 1 ] ; do
  sleep 30
  count_used_CPUs
 done
 raspa
 sleep 1
}
function go_lammps {
 count_used_CPUs
 while [ $(echo "${n_used} >= ${nCPU}" | bc -l) == 1 ] ; do
  sleep 30
  count_used_CPUs
 done
 lammps
 sleep 1
}
function count_used_CPUs {
# 
 n_used=0
 for process in "simulate" "lmp_mpi" "lmp_fftw" "gulp" ; do
  n=$(ps aux | grep ${process} | sed '/grep/d' | wc -l | awk '{print $1}')
  n_used=$((${n}+${n_used}))
 done
}
function lammps {
 time > salida.lammps
 nohup mpirun --np ${nCPU} ~/lammps/src/lmp_fftw -in in.lmp -sf opt >> salida.lammps 
 time >>  salida.lammps
}
function raspa {
 time > salida.raspa
 nohup simulate >> salida.raspa
 time >> salida.raspa
}
function lammps_raspa {
 cp ../lammpstrj2pdb .
 cp ../pdb2cif .
 ./lammpstrj2pdb < movs/opti.lammpstrj
 n_lines=$(wc -l out.pdb |awk '{print $1}')
 line=$(sed -n '/MODEL/{=;p}' out.pdb | sed '{N;s/\n/ /}' | tail -n1 | awk '{print $1}')
 end=$((n_lines - line + 1))
 tail -n$end out.pdb > input.pdb
 # Remove guest!!!
 if [ "${remove_guest}" == "true" ] ; then
  sed -i '/ Ar /d' input.pdb
  sed -i '/ Xe /d' input.pdb
  sed -i '/ Kr /d' input.pdb
 fi
 #
 ./pdb2cif
 mv p1.cif ${CyclesNameFile}.cif
 cp ${CyclesNameFile}.cif ../${CIFTemporallyFile}
 cp ${CyclesNameFile}.cif ../${CyclesNameFile}.cif
 rm input.pdb pdb2cif lammpstrj2pdb out.pdb
}
function CIF111toSupercell {
 check_supercell
 echo "SimulationType  MC
NumberOfCycles               0
NumberOfInitializationCycles 0
PrintEvery                   1
Forcefield  GenericMOFs
Framework 0
FrameworkName ${structure}_${seed}
UnitCells $ua $ub $uc" > simulation.input
 ${HOME}/RASPA/simulations/bin/simulate
 mv Movies/System_0/Framework_0_final_${ua}_${ub}_${uc}_P1.cif ${CIFTemporallyFile}
 sed -i '/^$/d' ${CIFTemporallyFile}
}
function distance_angle_measure {
 echo "${CIFTemporallyFile} 5
${CyclesNameFile}.Geometrical_Analysis.txt
1.7 2.6
1.25 1.60" > main.txt
 ./zif_dist_angle_02 main.txt > salida.zif_dist_angle_02.txt
 rm -rf Movies Output VTK Restart main.txt input.cif 
 max_dist_ZnN=$(cat ${CyclesNameFile}.Geometrical_Analysis.txt | grep "ZnN_distances" | awk '{print $4}')
 ave_dist_ZnN=$(cat ${CyclesNameFile}.Geometrical_Analysis.txt | grep "ZnN_distances" | awk '{print $3}')
 overall_goodness=$(cat ${CyclesNameFile}.Geometrical_Analysis.txt | grep "Overall_goodness" | awk '{print $2}')
 resumen
}
function resumen {
 echo "============================="
 echo "Cycle: ${cycle_name}"
 echo " # Loading ${guest}: $n_Ar molecules [$(echo "$n_Ar > $rrr" | bc -l)]"
 echo " # Energy (Total - Empty Structure / Nmolecules): $statu [eV]"
 echo " # ( Total: $energy [eV], Empty: $energy0 [eV] )"
 echo " # Geometry:"
 echo "   1. max. Zn...N distance: $max_dist_ZnN [A]"
 echo "   2. ave. Zn...N distance: $ave_dist_ZnN [A]"
 echo "   3. Overall Goodness:     $overall_goodness [-]"
 echo "============================="
}
##############################################################
# main program:
main_folder=${structure}_${temperature}_${pressure}_${seed}
mkdir ${main_folder}
cd ${main_folder}
 cp ${lib_folder}/forcefield.lib .
 cp ../${CIFFile} ${CIFTemporallyFile}
 cycle=0
 n_Ar=0
 cycle_name=$(echo $cycle | awk '{ printf("%02d\n", $1) }')
 make_binaries
 guest='empty'
 first_optimisation
 distance_angle_measure
 energy0=$(tail -n1 ${CyclesNameFile}_emmd/logs/minimization.txt | awk '{print $5}')
 for i in $(seq 1 ${n_cycles}) ; do
  # 
  guest='argon'
  let cycle++
  cycle_name=$(echo $cycle | awk '{ printf("%02d\n", $1) }')
  fill_with_guest
  # 
  flags_cif2lammps="post-loading"
  interface_adsorption_lammps
  remove_guest="false"
  em_md_lammps
  energy=$(tail -n1 ${CyclesNameFile}_emmd/logs/minimization.txt | awk '{print $5}')
  statu=$(echo "scale=4; ($energy - $energy0)/(${n_Ar})" | bc -l)
  distance_angle_measure
  previous_name=${CyclesNameFile}
  guest='xenon'
  cycle_name=$(echo $cycle | awk '{ printf("%02d\n", $1) }')
  update_name
  sed "s/Ar /Xe /g" ${previous_name}.cif > ${CyclesNameFile}.cif
  #
  remove_guest="true"
  flags_cif2lammps="post-Xe-Ar-exchange"
  interface_adsorption_lammps
  em_md_lammps
  energy=$(tail -n1 ${CyclesNameFile}_emmd/logs/minimization.txt | awk '{print $5}')
  statu=$(echo "scale=4; ($energy - $energy0)/(${n_Ar})" | bc -l)
  distance_angle_measure
 done
 clean_binaries
cd ..
exit 0
