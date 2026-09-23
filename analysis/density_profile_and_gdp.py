#!/usr/bin/env python
# coding: utf-8
"""
Compute mass-density profiles (MAICoS) for a water-ethanol/salt mixture
confined between solid walls, and locate the Gibbs dividing planes (GDP)
from the liquid density profile.
"""

import numpy as np
import MDAnalysis as mda
import maicos

# ---- inputs ----------------------------------------------------------
tpr = "path-to-tpr"
xtc = "path-to-xtc"

bin_width = 0.1
zmin, zmax = -55, 55


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


def GDP(z, rho_liq):
    """Locate the two Gibbs dividing planes from a liquid density profile."""
    z_l, z_r = z[: len(z) // 2], z[len(z) // 2 :]
    rho_l, rho_r = rho_liq[: len(z) // 2], rho_liq[len(z) // 2 :]

    # assume the middle of the density profile is bulk density
    rho_bulk = np.mean(rho_liq[int(len(z) * 0.4) : int(len(z) * 0.6)])

    # right boundary
    ref = 1e5
    zlim_r = None
    for cpt_z, z_i in enumerate(z_r):
        A = np.trapz(rho_bulk - rho_r[:cpt_z], z_r[:cpt_z])
        B = np.trapz(rho_r[cpt_z:], z_r[cpt_z:])
        diff = np.abs(A - B)
        if diff < ref:
            ref = diff
            zlim_r = z_i

    # left boundary
    ref = 1e5
    zlim_l = None
    for cpt_z, z_i in enumerate(z_l):
        A = np.trapz(rho_l[:cpt_z], z_l[:cpt_z])
        B = np.trapz(rho_bulk - rho_l[cpt_z:], z_l[cpt_z:])
        diff = np.abs(A - B)
        if diff < ref:
            ref = diff
            zlim_l = z_i

    return zlim_l, zlim_r


# ---- load trajectory ---------------------------------------------------
u = mda.Universe(tpr, xtc)
n_frames = u.trajectory.n_frames

liquids, liquid_names, solids, solid_names = detect_groups(u)
liquid_ref = liquids[liquid_names.index("Liq")]

# ---- density profile per species ---------------------------------------
for group, name in zip(liquids, liquid_names):
    if group.n_atoms == 0:
        continue

    planar = maicos.DensityPlanar(
        group, refgroup=liquid_ref, dens="mass",
        dim=2, bin_width=bin_width, unwrap=False,
        zmin=zmin, zmax=zmax, jitter=0.01,
    ).run()

    z = planar.results["bin_pos"]
    profile = planar.results["profile"]

    np.savetxt(f"density_{name}.dat", np.column_stack([z, profile]),
               header="z density")

    if name == "Liq":
        zlim_l, zlim_r = GDP(z, profile.T[0])
        np.savetxt("gdp.dat", [[zlim_l, zlim_r]], header="zlim_l zlim_r")
    