#!/usr/bin/env python
# coding: utf-8
"""
Build 2D in-plane density maps of surface species at each liquid-solid
interface (top and bottom), using pytim to identify the outermost layers.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import MDAnalysis as mda
import pytim

# ---- inputs ----------------------------------------------------------
tpr = "path-to-tpr"
xtc = "path-to-xtc"

matrix_size = (100, 100)

species_of_interest = [
    ["oxwat", "name OW"],
    ["careth", "name Co Ch"],
    ["eth", "resname Eth"],
    ["oxcareth", "name Co Ch Oh"],
    ["oxeth", "name Oh"],
    ["wat", "resname Sol"],
]


def detect_groups(u,
                   liquid_names=["Sol", "Eth", "Sol Eth Na Cl", "Na", "Cl"],
                   liquid_types=["4 5", "7 8 9 10 11 12", "2 3 4 5 7 8 9 10 11 12", "2", "3"],
                   solid_names=["Sil", "Gra"],
                   solid_types=["0", "1"],
                   lammps=False):
    labels = ["Sol", "Eth", "Liq", "Na", "Cl"]

    if not lammps:
        liquids = [u.select_atoms("resname " + names) for names in liquid_names]
        solids = [u.select_atoms("resname " + names) for names in solid_names]
    else:
        liquids = [u.select_atoms("type " + types) for types in liquid_types]
        solids = [u.select_atoms("type " + types) for types in solid_types]

    return liquids, labels, solids, solid_names


# ---- load trajectory ---------------------------------------------------
u = mda.Universe(tpr, xtc)
u.transfer_to_memory()

liquids, liquid_names, solids, solid_names = detect_groups(u)
liquid_ref = liquids[liquid_names.index("Liq")]

gCH3_1 = u.select_atoms("name Cm and prop z<25")
gCH3_2 = u.select_atoms("name Cm and prop z>25")

Lx, Ly = u.dimensions[0:2]
n_frames = u.trajectory.n_frames

# ---- initialize accumulators --------------------------------------------
dict_results_1 = {name: np.zeros(matrix_size) for name, _ in species_of_interest}
dict_results_2 = {name: np.zeros(matrix_size) for name, _ in species_of_interest}
dict_results_1["CH3"] = np.zeros(matrix_size)
dict_results_2["CH3"] = np.zeros(matrix_size)

# ---- accumulate 2D histograms per frame ---------------------------------
for ts in u.trajectory:
    interface = pytim.ITIM(u, group=liquid_ref, max_layers=4, molecular=True)
    surface_1 = interface.layers[1][0]
    surface_2 = interface.layers[0][0]

    for name, criteria in species_of_interest:
        group_1 = u.select_atoms(criteria).intersection(surface_1)
        counts_1, _, _ = np.histogram2d(
            group_1.positions[:, 0], group_1.positions[:, 1],
            bins=matrix_size, range=[[0, Lx], [0, Ly]]
        )
        dict_results_1[name] += counts_1

        group_2 = u.select_atoms(criteria).intersection(surface_2)
        counts_2, _, _ = np.histogram2d(
            group_2.positions[:, 0], group_2.positions[:, 1],
            bins=matrix_size, range=[[0, Lx], [0, Ly]]
        )
        dict_results_2[name] += counts_2

    counts_1, _, _ = np.histogram2d(
        gCH3_1.positions[:, 0], gCH3_1.positions[:, 1],
        bins=matrix_size, range=[[0, Lx], [0, Ly]]
    )
    dict_results_1["CH3"] += counts_1

    counts_2, _, _ = np.histogram2d(
        gCH3_2.positions[:, 0], gCH3_2.positions[:, 1],
        bins=matrix_size, range=[[0, Lx], [0, Ly]]
    )
    dict_results_2["CH3"] += counts_2

# ---- normalize and save --------------------------------------------------
for name in dict_results_1:
    dict_results_1[name] /= n_frames
    dict_results_2[name] /= n_frames

np.save("surface_density_map_bottom.npy", dict_results_1)
np.save("surface_density_map_top.npy", dict_results_2)
