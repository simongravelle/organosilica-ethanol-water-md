#!/usr/bin/env python
# coding: utf-8
"""
For ethanol's hydrophobic carbons (Ch, Co), count neighboring carbons
(ethanol carbons + hydrophobic surface carbons) per interfacial layer
(ITIM), and build the corresponding radial distribution functions.
"""

import numpy as np
import MDAnalysis as mda
import pytim

# ---- inputs ----------------------------------------------------------
tpr = "path-to-tpr"
xtc = "path-to-xtc"
step = 10

r_min, r_max = 2, 15       # A, excludes intramolecular distances
neighbor_cutoff = 4.54     # A
n_bins = 150
n_regions = 5               # 4 ITIM layers + bulk


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


def carbons_by_region(interface, carbon_ids, u):
    """Split carbon ids into the 4 ITIM layers + bulk, as atom groups."""
    groups = []
    ids_by_layer = []
    for i in range(4):
        ids_i = (interface.layers[0][i][(interface.layers[0][i].names == "Ch")
                                         | (interface.layers[0][i].names == "Co")].ids.tolist()
                 + interface.layers[1][i][(interface.layers[1][i].names == "Ch")
                                           | (interface.layers[1][i].names == "Co")].ids.tolist())
        ids_by_layer.append(ids_i)
        groups.append(u.select_atoms("id " + " ".join(map(str, ids_i))) if ids_i else u.select_atoms("id -1"))

    flat_ids = [i for layer in ids_by_layer for i in layer]
    bulk_ids = [int(i) for i in carbon_ids if i not in flat_ids]
    groups.append(u.select_atoms("id " + " ".join(map(str, bulk_ids))) if bulk_ids else u.select_atoms("id -1"))
    return groups


def count_neighbors(mygroup, other_positions, box, counts_bin):
    """Count neighbors of mygroup within [r_min, r_max) in other_positions,
    accumulating a pairwise-distance histogram into counts_bin, and return
    the number of neighbors within neighbor_cutoff."""
    n_neighbors = 0
    for p in mygroup.positions:
        d = (np.remainder(p - other_positions + box[:3] / 2., box[:3]) - box[:3] / 2.).T
        r = np.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
        n_neighbors += np.sum((r < neighbor_cutoff) & (r > r_min))
        r_in_range = r[(r > r_min) & (r < r_max)]
        for ri in r_in_range:
            counts_bin[int(ri * 10)] += 1
    return n_neighbors


# ---- load trajectory ---------------------------------------------------
u = mda.Universe(tpr, xtc)
u.transfer_to_memory(step=step)
n_frames = u.trajectory.n_frames

liquids, liquid_names, solids, solid_names = detect_groups(u)
liquid_ref = liquids[liquid_names.index("Liq")]
eth = liquids[liquid_names.index("Eth")]

surface = u.select_atoms("resname Sil and type C*")
carbon_ids = eth.ids[(eth.names == "Ch") | (eth.names == "Co")]

neighbor_count = np.zeros(n_regions)
carbon_count = np.zeros(n_regions)
counts = [np.zeros(n_bins) for _ in range(n_regions)]

# ---- accumulate over the trajectory --------------------------------------
for ts in u.trajectory:
    box = ts.dimensions
    interface = pytim.ITIM(u, group=liquid_ref, max_layers=4, molecular=True)
    groups = carbons_by_region(interface, carbon_ids, u)

    for i, group in enumerate(groups):
        if group.n_atoms == 0:
            carbon_count[i] += 0
            continue
        n_neighbors = count_neighbors(group, group.positions, box, counts[i])
        if surface.n_atoms > 0:
            n_neighbors += count_neighbors(group, surface.positions, box, counts[i])
        neighbor_count[i] += n_neighbors
        carbon_count[i] += group.n_atoms

# ---- normalize ------------------------------------------------------------
neighbor_count = neighbor_count / n_frames / 2  # each pair counted twice
carbon_count = carbon_count / n_frames
molecule_count = carbon_count / 2  # 2 hydrophobic carbons per ethanol

gofrs = []
r_vec = np.linspace(0, r_max, n_bins)
shell_volumes = (4 / 3) * np.pi * (r_vec[1:] ** 3 - r_vec[:-1] ** 3)
for count in counts:
    count = count / n_frames
    gofrs.append(count[1:] / shell_volumes)

# ---- save ----
np.savetxt("carbon_neighbors.dat", np.vstack([molecule_count, neighbor_count, carbon_count]).T,
           header="molecule_count neighbor_count carbon_count")
np.savetxt("carbon_rdf.dat", np.column_stack([r_vec[1:]] + gofrs),
           header="r " + " ".join(f"region{i}" for i in range(n_regions)))
np.savetxt("carbon_rdf_counts.dat", np.column_stack([r_vec] + counts),
           header="r " + " ".join(f"region{i}" for i in range(n_regions)))
