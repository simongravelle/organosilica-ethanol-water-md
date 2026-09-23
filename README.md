# Molecular Dynamics Simulations Scripts for the Study of Water-Ethanol Mixture Adsorption Selectivity at Organosilica Surfaces

This repository accompanies the manuscript:

**Adsorption Selectivity of Water-Ethanol Mixtures on Organosilica Surfaces: Role of Hydrophilicity**
S. Gravelle and B. Coasne

## Overview

Water-ethanol mixtures exhibit preferential ethanol enrichment at non-polar
interfaces, a phenomenon relevant to microporous membrane
separation. This repository contains the GROMACS and LAMMPS input files
used to study this enrichment at organosilica surfaces
with tunable hydrophobicity (varying the ratio of methyl to hydroxyl
surface groups, $r_{C/O}$).

## GROMACS input files

Input files for each surface type (Q2, Q3, Q4, Q2/Q3, Q3/Q4, amorphous, and
liquid-vapor) are provided under `simulations/`. For each folder, force field
parameters are provided in the `ff/` subfolder, and inputs are provided in the
`inputs/` subfolder. The topology corresponds to the `conf.gro` file. Inputs
are compatible with GROMACS version 2025.3.

```
gmx grompp -f input/prod.mdp -p topol-posres.top -o prod -pp prod -po prod -r conf.gro -maxwarn 2
gmx mdrun -deffnm prod -v -rdd 1 -nt 8 -pin on
```

## LAMMPS input files

Input files for the [graphite](graphite) slit pores are provided and should
be run as follows. These inputs are compatible with LAMMPS version 2Aug2023.

```
lmp -in input.lmp
```