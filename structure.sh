#!/bin/R

# "population_file" is a file which each row is indiviual. First column is the Individual ID and the second coulmn is location ID.

for((h=1; h<=3; h++))
do
for i in harded5_00101
do
for j in population_file
do
for k in harded5_00101
do
Rscript strc_prep.R ${i}

sed -e '1d' gt.contra > headcut
rm gt.contra
less headcut | cut -d " " -f 2- > GTTemp
rm headcut
perl -pe 's/0/0 0/g; s/1/1 0/g; s/2/1 1/g; s/-9/-9 -9/g;' GTTemp > GTTempst
paste -d " " ${j} GTTempst > GTTempinput

structure -K 1 -i GTTempinput -o ${i}.1 ${k}.${h}
structure -K 2 -i GTTempinput -o ${i}.2 ${k}.${h}
structure -K 3 -i GTTempinput -o ${i}.3 ${k}.${h}
structure -K 4 -i GTTempinput -o ${i}.4 ${k}.${h}
structure -K 5 -i GTTempinput -o ${i}.5 ${k}.${h}
structure -K 6 -i GTTempinput -o ${i}.6 ${k}.${h}
structure -K 7 -i GTTempinput -o ${i}.7 ${k}.${h}
structure -K 8 -i GTTempinput -o ${i}.8 ${k}.${h}
structure -K 9 -i GTTempinput -o ${i}.9 ${k}.${h}
structure -K 10 -i GTTempinput -o ${i}.10 ${k}.${h} 

structureHarvester.py --dir=path/to/dir/ --out=path/to/dir/ --evanno

done
done
done
done