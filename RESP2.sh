# A script to calculate RESP2 charges by invoking Gaussian and Multiwfn
# Written by Tian Lu (sobereva@sina.com)
# Last update: 2020-Feb-1
# RESP2(0.5) for singlet neutral molecule with water solvent: ./RESP2.sh maki.pdb
# RESP2(0.5) for triplet neutral molecule with water solvent: ./RESP2.sh nozomi.xyz 0 3
# RESP2(0.5) for singlet anion with ethanol solvent: ./RESP2.sh nico.mol -1 1 ethanol


#Edited by YH. W 
#Edit Date: 2021/05/18

#!/bin/bash
delta=0.5
level_opt="B3LYP/def2SVP em=GD3BJ"
level_SP="B3LYP/def2TZVP em=GD3BJ"
Gaussian=g16

export inname=$1
filename=${inname%.*}
suffix=${inname##*.}

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

keyword_opt="# "$level_opt" opt=loose "$solvent
keyword_SP_gas="# "$level_SP" pop=MK IOp(6/33=2,6/42=6)"
keyword_SP_solv="# "$level_SP" "$solvent" pop=MK IOp(6/33=2,6/42=6)"

#### Convert input file to .xyz file
Multiwfn $1 > /dev/null << EOF
100
2
2
tmp.xyz
0
q
EOF

#Optimize geometry
cat << EOF > gau.gjf
%chk=gau.chk
$keyword_opt

test

$chg $multi
EOF
awk '{if (NR>2) print }' tmp.xyz >> gau.gjf
cat << EOF >> gau.gjf


EOF
rm tmp.xyz

echo
echo Running optimization task under solvent via Gaussian...
$Gaussian < gau.gjf > gau.out
if grep -Fq "Normal termination" gau.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi

#### Single point in gas
cat << EOF > gau.gjf
%chk=gau.chk
$keyword_SP_gas geom=allcheck guess=read


EOF

echo
echo Running single point task in gas phase via Gaussian...
$Gaussian < gau.gjf > gau.out

if grep -Fq "Normal termination" gau.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi

echo Running formchk...
formchk gau.chk > /dev/null

echo Running Multiwfn...
Multiwfn gau.fchk > /dev/null << EOF
7
18
8
1
gau.out
y
0
0
q
EOF

mv gau.chg gas.chg
echo RESP charge in gas phase has been outputted to gas.chg

#### Single point in solvent
cat << EOF > gau.gjf
%chk=gau.chk
$keyword_SP_solv geom=allcheck guess=read


EOF

echo
echo Running single point task in solvent phase via Gaussian...
$Gaussian < gau.gjf > gau.out

if grep -Fq "Normal termination" gau.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi

echo Running formchk...
formchk gau.chk > /dev/null

echo Running Multiwfn...
Multiwfn gau.fchk > /dev/null << EOF
7
18
8
1
gau.out
y
0
0
q
EOF

mv gau.chg solv.chg
echo RESP charge in solvent phase has been outputted to solv.chg

#### Calculate RESP2
chgname=${filename}".chg"
pqrname=${filename}".pqr"

paste gas.chg solv.chg |awk '{printf $1 " " $2 " " $3 " " $4 " " (1-d)*$5+d*$10 "\n"}' d=$delta > $chgname

Multiwfn $chgname > /dev/null << EOF
100
2
1
tmp.pqr
0
q
EOF

mv tmp.pqr $pqrname


echo
echo Finished! The optimized atomic coordinates with RESP2 charges \(the last column\) have been exported to $pqrname in current folder
