#!/bin/bash
#OAR -n mixture-nanopore
#OAR -l /nodes=1/gpu=1/cpu=1/core=8,walltime=48:00:00
#OAR -p gpumodel='A100'
#OAR --stdout log.out
#OAR --stderr log.err
#OAR --project tamtam

set -e

export GMX_MAXBACKUP=-1

gmx=/home/gravells/softwares/gromacs-2023/build-gpu/bin/gmx

# create conf and topol
${gmx} genconf -f conf.gro -o conf-nemd.gro -nbox 1 1 2

cp topol.top topol-nemd.top
tail -n 3 "topol-nemd.top" >> "topol-nemd.top"
sed -i '/silica./c\#include "ff/silica.itp" \n#include "ff/posres_silica_nemd.itp"' topol-posres.top

# detect silica resid
sil1=1
id1=1
while IFS= read -r line
do
    resid=${line:0:5}
    res="${line:5:5}"
    id="${line:15:5}"
    if [ ${res} = "Sil" ];
    then
        sil2=${resid}
        id2=${id}
    fi
done < conf-nemd.gro
echo ${sil1}-${sil2}
echo ${id1}-${id2}

# update mdp files
newline='pull-group1-pbcatom = '$id1
oldline=$(cat input/nemd.mdp | grep 'pull-group1-pbcatom =')
sed -i '/'"$oldline"'/c\'"$newline" input/nemd.mdp
newline='pull-group2-pbcatom = '$id2
oldline=$(cat input/nemd.mdp | grep 'pull-group2-pbcatom =')
sed -i '/'"$oldline"'/c\'"$newline" input/nemd.mdp

newline='pull-group1-pbcatom = '$id1
oldline=$(cat input/nemd-eq.mdp | grep 'pull-group1-pbcatom =')
sed -i '/'"$oldline"'/c\'"$newline" input/nemd-eq.mdp
newline='pull-group2-pbcatom = '$id2
oldline=$(cat input/nemd-eq.mdp | grep 'pull-group2-pbcatom =')
sed -i '/'"$oldline"'/c\'"$newline" input/nemd-eq.mdp

# create index
${gmx} make_ndx -f conf-nemd.gro<<EOF
    r ${sil1} 
    name 5 Sil1
    r ${sil2}
    name 6 Sil2
    q
EOF

${gmx} grompp -f input/nemd-eq.mdp -c conf-nemd.gro -p topol-nemd.top -o nemd-eq -pp nemd-eq -po nemd-eq -n index.ndx -r conf-nemd.gro -maxwarn 2
${gmx} mdrun -deffnm nemd-eq -v -rdd 1 -nt 8 -pin on
cp nemd-eq.gro conf-nemd.gro

${gmx} grompp -f input/nemd.mdp -c conf-nemd.gro -p topol-nemd.top -o nemd -pp nemd -po nemd -n index.ndx -r conf-nemd.gro -maxwarn 2
${gmx} mdrun -deffnm nemd -v -rdd 1 -nt 8 -pin on
cp nemd.gro conf-nemd.gro
