#!/bin/bash
function calc()
{
	python -c "print($*)"
}

#source /home/apps/ECPM2020/V3.7/bin/GMXRC
mgbin="/home/apps/matrixGenerator/itptoolsV4.0 "
tprfile="./topol.tpr" ## position of tpr

Nele=1360 ##Number of electrode atoms(charge)
Nzero=2720 ##Number of electrode atoms(no charge)

halfN=$(calc "$Nele"/2 )
$mgbin  matrix -s $tprfile  -n $Nele   #get matrix file
#$mgbin  getGCD  -f allMatrixA.bin   -o allGcdMatrixA.bin -n $Nele

Vs=(0 2) ##The voltage between two electrodes
for V in ${Vs[@]} ; do 
	halfV=$( calc "$V"*0.5 )
	$mgbin newCPM_File -pot  [  -"$halfV" "$halfV" ] -npot [ $halfN $halfN ] -nzero "$Nzero" -outQ 20 -bulk [ 24.0 26.0 ]
	mv CPM_ControlFile.dat  CPM_ControlFile.dat_"$V"V
done ##get control file 



