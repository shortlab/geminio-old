#!/bin/bash


: << 'END'
#######################create folders for each energy############
#energies=(0.1 0.15 0.2 0.3 0.5 0.75 1  1.5 2  3  5  7.5 10 15 20 30 40 50 60 75 100 120 150)
#halfL=(   15  15   15  15  15  15   15 20  20 20 30 30  30 30 30 30 40 40 50 50 60  65  65)
energies=(0.1 0.14 0.22 0.35 0.55 0.86  1.34 2.1 3.3 5.13 8.02 12.5  19.5 30 46.9 72.7 107 171.9)
halfL=(   15  15   15   15   15   15    20   20  20  30   30   30    40   40 50   60   60  65)
#energies=(0.1 0.13 0.22 0.35 0.65 0.86 1.27  2  3.3  4  6  8  10 12.5 20 23 45 55 75 107)
#halfL=(   15  15   15   15   15   20   20    20 30   30 30 30 30 30   30 40 40 50 50 60)
#energies=(0.15 0.2 0.35 0.63 0.75 1.34 2.1 3  4  6.3 8.02 10 12.6 19.5 23 30.4 45 55 72.7 100)
#halfL=(   15   15  15   15   15   20   20  20 30 30  30   30 30   40   40 40   40 50 60   65)
counter=0
for energy in ${energies[@]};
do
   datafolder=30K${energy}keV
   if [ -d "${datafolder}" ]; then
     counter=$((counter+1))
     continue
   fi
   mkdir ${datafolder}
   cd ${datafolder}
   cp ../run.sh run_${energy}.sh
   chmod +x run_${energy}.sh
  
   python ../generate_random_directions.py ${halfL[counter]} ${energy}
   
   cp ../run.pbs run_${energy}.pbs
   sed -i "s/run.sh/run_${energy}.sh/g" run_${energy}.pbs
   sed -i "s/RUN/RUN_${energy}/g" run_${energy}.pbs
   qsub run_${energy}.pbs
   
   counter=$((counter+1))
   cd ..
done
END

: << 'END'
#######################create folders for each energy############
#energies=(250 400)
#halfL=(   70  80)
energies=(12.6)
halfL=(   30 )
counter=0
for energy in ${energies[@]};
do
   datafolder=${energy}keV
   if [ ! -d "${datafolder}" ]; then
     mkdir ${datafolder}
   fi
   cd ${datafolder}
   cp ../run.sh run_${energy}.sh
   chmod +x run_${energy}.sh
  
   python ../generate_random_directions.py ${halfL[counter]} ${energy}
   
   cp ../run.pbs run_${energy}.pbs
   sed -i "s/run.sh/run_${energy}.sh/g" run_${energy}.pbs
   sed -i "s/RUN/RUN_${energy}/g" run_${energy}.pbs
   qsub run_${energy}.pbs
   
   counter=$((counter+1))
   cd ..
done
END

: << 'END'
#######################create folders for each energy############
energies=(0.16 0.25 0.4 0.63  1.58 2.51 4  6.3 15.8  25.1  63)
halfL=(   15   15   15  15    20   20   20 30  30    40    50)
counter=0
for energy in ${energies[@]};
do
   datafolder=${energy}keV
   if [ ! -d "${datafolder}" ]; then
     mkdir ${datafolder}
   fi
   cd ${datafolder}
   cp ../run.sh run_${energy}.sh
   chmod +x run_${energy}.sh
  
   python ../generate_random_directions.py ${halfL[counter]} ${energy}
   
   cp ../run.pbs run_${energy}.pbs
   sed -i "s/run.sh/run_${energy}.sh/g" run_${energy}.pbs
   sed -i "s/RUN/RUN_${energy}/g" run_${energy}.pbs
   qsub run_${energy}.pbs
   
   counter=$((counter+1))
   cd ..
done
END


: << 'END'
#######################detect clusters for each energy############
#energies=(0.1 0.15 0.2 0.3 0.5 0.75 1  1.5 2  3  5  7.5 10 15 20 30 40 50 60 75 100 120 150)
#halfL=(   15  15   15  15  15  15   15 20  20 20 30 30  30 30 30 30 40 40 50 50 60  65  65)
#energies=(0.16 0.25 0.4 0.63  1.58 2.51 4  6.3 15.8  25.1  63)
#halfL=(   15   15   15  15    20   20   20 30  30    40    50)
#energies=(0.1 0.16 0.25 0.4 0.63 1  1.58 2.51 4  6.3 10 15.8  25.1 40 63 100 150 250 400)
#halfL=(   15  15   15   15  15   15 20   20   20 30  30 30    40   40 50 60  65  64  70)
#energies=(250 400)
#halfL=(   64  70)
#energies=(12.6)
#halfL=(   30 )
energies=(0.1 0.15 0.2 0.35 0.63 0.75 1.34 2.1 3  4  6.3 8.02 10 12.6 19.5 23 30.4 45 55 72.7 100)
halfL=(   15  15   15  15   15   15   20   20  20 30 30  30   30 30   40   40 40   40 50 60   60)
#energies=(0.1 0.14 0.22 0.35 0.55 0.86  1.34 2.1 3.3 5.13 8.02 12.5  19.5 30.4 46.9 72.7 106.6 171.9)
#halfL=(   15  15   15   15   15   15    20   20  20  30   30   30    40   40   50   60   65    65)
#energies=(0.1 0.13 0.22 0.35 0.65 0.86 1.27  2  3.3  4  6  8  10 12.5 20 23 45 55 75 107)
#halfL=(   15  15   15   15   15   20   20    20 30   30 30 30 30 30   30 40 40 50 50 60)
#energies=(0.1 0.14 0.22 0.35 0.55 0.86  1.34 2.1 3.3 5.13 8.02 12.5  19.5 30 46.9 72.7 107 171.9)
#halfL=(   15  15   15   15   15   15    20   20  20  30   30   30    40   40 50   60   60  65)
counter=0
for energy in ${energies[@]};
do
   echo START energy: ${energy}
   datafolder=${energy}keV
   cd ${datafolder}
   ~/projects/lammps-16Feb16/mytools/ovito-2.8.1-x86_64/bin/ovitos ../detect_clusters.py ${halfL[counter]}
   counter=$((counter+1))
   cd ..
   echo END energy: ${energy}
done
END

: << 'END'
#######################move files for each energy############
energies=(0.1 0.15 0.2 0.3 0.5 0.75 1  1.5 2  3  5  7.5 10 15 20 30 40 50 60 75 100 120 150)
halfL=(   15  15   15  15  15  15   15 20  20 20 30 30  30 30 30 30 40 40 50 50 60  65  65)
counter=0
cd download/
for energy in ${energies[@]};
do
   echo START energy: ${energy}
   datafolder=${energy}keV
   mkdir ${datafolder}
   cp ../${datafolder}/cluster_frequency.txt ${datafolder}/
   
   counter=$((counter+1))
   echo END energy: ${energy}
done
END

#: << 'END'
#######################move files for each energy############
#energies=(0.1 0.16 0.25 0.4 0.63 1  1.58 2.51 4  6.3 10 15.8  25.1 40 63 100 150 250 400)
#halfL=(   15  15   15   15  15   15 20   20   20 30  30 30    40   40 50 60  65  64  70)
#energies=(12.6)
#halfL=(   30 )
#energies=(0.1 0.15 0.2 0.3 0.5 0.75 1  1.5 2  3  5  7.5 10 12.6 15 20 30 40 50 60 75 100 120 150)
#halfL=(   15  15   15  15  15  15   15 20  20 20 30 30  30 30   30 30 30 40 40 50 50 60  65  65)
#energies=(0.15 0.2 0.35 0.63 0.75 1.34 2.1 3  4  6.3 8.02 10 12.6 19.5 23 30.4 45 55 72.7 100)
#halfL=(   15   15  15   15   15   20   20  20 30 30  30   30 30   40   40 40   40 50 60   60)
#energies=(0.1 0.14 0.22 0.35 0.55 0.86  1.34 2.1 3.3 5.13 8.02 12.5  19.5 30.4 46.9 72.7 106.6 171.9)
#halfL=(   15  15   15   15   15   15    20   20  20  30   30   30    40   40   50   60   65    65)
#energies=(0.1 0.13 0.22 0.35 0.65 0.86 1.27  2  3.3  4  6  8  10 12.5 20 23 45 55 75 107)
#halfL=(   15  15   15   15   15   20   20    20 30   30 30 30 30 30   30 40 40 50 50 60)
energies=(0.1 0.14 0.22 0.35 0.55 0.86  1.34 2.1 3.3 5.13 8.02 12.5  19.5 30 46.9 72.7 107 171.9)
halfL=(   15  15   15   15   15   15    20   20  20  30   30   30    40   40 50   60   60  65)
#energies=(0.1 0.15 0.2 0.35 0.63 0.75 1.34 2.1 3  4  6.3 8.02 10 12.6 19.5 23 30.4 45 55 72.7 100)
#halfL=(   15  15   15  15   15   15   20   20  20 30 30  30   30 30   40   40 40   40 50 60   60)
counter=0
cd energies_result/
for energy in ${energies[@]};
do
   echo START energy: ${energy}
   datafolder=30K${energy}keV
   if [ ! -d ${datafolder} ] ; then
     mkdir ${datafolder}
   fi
   cp ../${datafolder}/*cluster_frequency.txt ${datafolder}/
   counter=$((counter+1))
   echo END energy: ${energy}
done
#END
