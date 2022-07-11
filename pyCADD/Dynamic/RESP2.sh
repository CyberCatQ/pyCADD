# A script to calculate RESP2 charges by invoking Gaussian and Multiwfn
# Written by Tian Lu (sobereva@sina.com)
# Last update: 2020-Feb-1
# RESP2(0.5) for singlet neutral molecule with water solvent: ./RESP2.sh maki.pdb
# RESP2(0.5) for triplet neutral molecule with water solvent: ./RESP2.sh nozomi.xyz 0 3
# RESP2(0.5) for singlet anion with ethanol solvent: ./RESP2.sh nico.mol -1 1 ethanol


#Edited by YH. W 
#Edit Date: 2022/07/08

#!/bin/sh
delta=0.5
level_opt="B3LYP/def2SVP em=GD3BJ"
level_SP="B3LYP/def2TZVP em=GD3BJ"
Gaussian=g16

export inname=$1
filepath=$inname
filename=${filepath##*/}
suffix=${filename##*.}
prefix=${filename%%.*}

if [ $2 ];then
	echo "Net charge = $2"
	chg=$2
else
	echo "Net charge was not defined. Default to 0"
	chg=0
fi

if [ $3 ];then
	echo "Spin multiplicity = $3"
	multi=$3
else
	echo "Spin multiplicity was not defined. Default to 1"
	multi=1
fi

if [ $4 ];then
	echo Solvent is $4
	solvent="scrf(solvent="$4")"
else
	solvent="scrf(solvent=water)"
	echo "Solvent name was not defined. Default to water"
fi

echo delta parameter is $delta

keyword_opt="# "$level_opt" opt=loose "$solvent" nosymm"
keyword_SP_gas="# "$level_SP" pop=MK IOp(6/33=2,6/42=6) nosymm"
keyword_SP_solv="# "$level_SP" "$solvent" pop=MK IOp(6/33=2,6/42=6) nosymm"

#### Convert input file to .xyz file
Multiwfn $1 > /dev/null << EOF
100
2
2
origin.xyz
0
q
EOF

#Optimize geometry
cat << EOF > opt.gjf
%chk=opt.chk
$keyword_opt

structure optimization

$chg $multi
EOF
awk '{if (NR>2) print }' origin.xyz >> opt.gjf
cat << EOF >> opt.gjf


EOF

echo
echo Running optimization task under solvent via Gaussian...
$Gaussian < opt.gjf > opt.out
if grep -Fq "Normal termination" opt.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi

formchk opt.chk > /dev/null

#### Single point in gas
cp opt.chk gas.chk
cat << EOF > gas.gjf
%chk=gas.chk
$keyword_SP_gas geom=allcheck guess=read


EOF

echo
echo Running single point task in gas phase via Gaussian...
$Gaussian < gas.gjf > gas.out

if grep -Fq "Normal termination" gas.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi

echo Running formchk...
formchk gas.chk > /dev/null

echo Running Multiwfn...
Multiwfn gas.fchk > /dev/null << EOF
7
18
8
1
gas.out
y
0
0
q
EOF

echo RESP charge in gas phase has been outputted to gas.chg

#### Single point in solvent
cp gas.chk solv.chk
cat << EOF > solv.gjf
%chk=solv.chk
$keyword_SP_solv geom=allcheck guess=read


EOF

echo
echo Running single point task in solvent phase via Gaussian...
$Gaussian < solv.gjf > solv.out

if grep -Fq "Normal termination" solv.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi

echo Running formchk...
formchk solv.chk > /dev/null

echo Running Multiwfn...
Multiwfn solv.fchk > /dev/null << EOF
7
18
8
1
solv.out
y
0
0
q
EOF

echo RESP charge in solvent phase has been outputted to solv.chg

#### Calculate RESP2
chgname=${prefix}".chg"
outputfile=${prefix}".pqr"

paste gas.chg solv.chg |awk '{printf $1 " " $2 " " $3 " " $4 " " (1-d)*$5+d*$10 "\n"}' d=$delta > $chgname

Multiwfn $chgname > /dev/null << EOF
100
2
1
tmp.pqr
0
q
EOF

mv tmp.pqr $outputfile


echo
echo Finished! The optimized atomic coordinates with RESP2 charges \(the last column\) have been exported to $outputfile
