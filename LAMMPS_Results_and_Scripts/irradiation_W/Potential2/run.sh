#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


inputfile='in.W_hybrid'
minfile='in.W_hybrid_min'
filename='velocities.txt'

halfL=0
cp ../${inputfile} .
cp ../${minfile} .
counter=0
while read velocities;
do
  line=(${velocities})
  halfL=${line[0]}
  vx=${line[1]}
  vy=${line[2]}
  vz=${line[3]}

  mpirun /home/jinmiao/projects/lammps-16Feb16/src/lmp_mpi -var boxL ${halfL} -var pkaVx ${vx} -var pkaVy ${vy} -var pkaVz ${vz} <  ${inputfile}  
  mv log.lammps log.lammps_${counter}
  mv dump_irra_minimized.data dump_irra_minimized_${counter}.data
  mv realtime.txt realtime_${counter}.txt
  counter=$((counter+1))

  ##################################################
  #mpirun -np 16 /home/mmjin/projects/lammps-5Sep14/src/lmp_linux -var xpos $xpos -var ypos ${ypos[i]} < in.radiation_second
done < ${filename} 

mpirun /home/jinmiao/projects/lammps-16Feb16/src/lmp_mpi -var boxL ${halfL}  <  ${minfile}  
