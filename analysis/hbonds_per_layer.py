#!/usr/bin/env python
# coding: utf-8
"""
Count hydrogen bonds per interfacial layer (ITIM) between a liquid species
(water or ethanol) and its environment, and normalize by the number of
molecules found in each layer.
"""

import numpy as np
import MDAnalysis as mda
import pytim
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

# ---- inputs ----------------------------------------------------------
tpr = "path-to-tpr"
xtc = "path-to-xtc"
step = 500
n_layers = 4  # + 1 bulk layer


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


def define_donor_acceptor(name, short=False):
    """Return candidate donor atoms for `name` and all possible acceptors."""
    donors = {"Sol": ["OW"], "Eth": ["Oh", "Co"]}
    if name not in donors:
        raise ValueError(f"Unknown group '{name}'")

    atoms1 = donors[name]
    atoms2 = ["OW", "Oh", "Co", "Cm", "Oy", "Cg"]

    if not short:
        atoms1 = [f"name {a}" for a in atoms1]
        atoms2 = [f"name {a}" for a in atoms2]

    return atoms1, atoms2


def find_combination(u, atoms1, atoms2):
    """Find valid donor/acceptor selection pairs and their O-H/O-C weight."""
    def pair_weight(name1, name2):
        if name1[-2] == "O" and name2[-2] == "O":
            return 1
        if {name1[-2], name2[-2]} == {"O", "C"}:
            return 0.4
        return 0

    combinations = []
    for name1 in atoms1:
        for name2 in atoms2:
            if do_analyse(u, name1, name2):
                combinations.append([name1, name2, pair_weight(name1, name2)])
    for name1 in atoms2:
        for name2 in atoms1:
            if name1 != name2 and do_analyse(u, name1, name2):
                combinations.append([name1, name2, pair_weight(name1, name2)])
    return combinations


def do_analyse(u, donor_sel, acceptor_sel):
    """Return True if both selections are non-empty for this system."""
    return (u.select_atoms(donor_sel).n_atoms > 0
            and u.select_atoms(acceptor_sel).n_atoms > 0)


def find_peak(donor, acceptor):
    """Distance cutoff (A) for the donor-acceptor pair. Adjust as needed."""
    return 3.5


def eval_hbond(u, donor, acceptor):
    hbonds = HBA(
        universe=u, donors_sel=donor, acceptors_sel=acceptor,
        d_h_cutoff=1.2, d_a_cutoff=find_peak(donor, acceptor), d_h_a_angle_cutoff=150,
    ).run()
    return hbonds.results["hbonds"]


def find_location(atom_id, interface):
    """Layer index (0-3) of the atom, or n_layers if it's in the bulk."""
    for n in range(n_layers):
        if (atom_id in interface.layers[0][n].atoms.ids
                or atom_id in interface.layers[1][n].atoms.ids):
            return n
    return n_layers


def assign_hb(hb, loc_donor, loc_acceptor, resname_donor, resname_acceptor):
    loc_donor += 1
    loc_acceptor += 1
    if resname_donor in ["Sil", "Gra"]:
        loc_donor = 0
    if resname_acceptor in ["Sil", "Gra"]:
        loc_acceptor = 0
    hb[loc_donor, loc_acceptor] += 1
    return hb


def hbonds_per_layer(u, liquid_ref, name, n_mol_group, n_frames):
    """Compute total and surface-only HB counts per layer for species `name`."""
    all_results = {}
    n_molecules = np.zeros(n_layers + 1)

    atoms1, atoms2 = define_donor_acceptor(name)
    combinations = find_combination(u, atoms1, atoms2)

    for name1, name2, pref in combinations:
        hb_name1_name2 = np.zeros((n_layers + 2, n_layers + 2), dtype="int32")
        frame_ids, donor_ids, _, acceptor_ids, _, _ = eval_hbond(u, name1, name2).T

        for ts in u.trajectory:
            interface = pytim.ITIM(u, group=liquid_ref, max_layers=n_layers, molecular=True)

            # molecule count per layer, this frame
            n_molecules_frame = np.zeros(n_layers + 1)
            for n, group1, group2 in zip(range(n_layers), interface.layers[0], interface.layers[1]):
                group_sum = group1 + group2  # assumes symmetric interfaces
                n_molecules_frame[n] = len(np.unique(group_sum[group_sum.resnames == name].resids))
                n_molecules[n] += n_molecules_frame[n]
            n_bulk = n_mol_group - np.sum(n_molecules_frame)
            n_molecules_frame[-1] = n_bulk
            n_molecules[-1] += n_bulk

            # hydrogen-bond count per layer, this frame
            cond = frame_ids == ts.frame
            for id_donor, id_acceptor in zip(donor_ids[cond].astype(int), acceptor_ids[cond].astype(int)):
                resname_donor = u.select_atoms(f"id {id_donor}").resnames[0]
                resname_acceptor = u.select_atoms(f"id {id_acceptor}").resnames[0]
                loc_donor = find_location(id_donor, interface)
                loc_acceptor = find_location(id_acceptor, interface)
                hb_name1_name2 = assign_hb(hb_name1_name2, loc_donor, loc_acceptor,
                                            resname_donor, resname_acceptor)

        if pref == 1:  # keep only O-H...O bonds
            all_results[name1[-2:] + "-" + name2[-2:]] = np.round(hb_name1_name2 * pref / n_frames, 2)

    n_molecules /= n_frames

    atoms_of_interest, _ = define_donor_acceptor(name, short=True)
    total_hb, surface_hb = [], []
    for layer_id in range(1, n_layers + 2):
        hb_layer_id = hb_surface_id = 0
        for couple, hb in all_results.items():
            donor_name, acceptor_name = couple[:2], couple[-2:]
            for donor_i in range(n_layers + 2):
                for acceptor_j in range(n_layers + 2):
                    n_ij = hb[donor_i, acceptor_j]
                    if donor_i == layer_id and donor_name in atoms_of_interest:
                        hb_layer_id += n_ij
                        if acceptor_j == 0:
                            hb_surface_id += n_ij
                    if acceptor_j == layer_id and acceptor_name in atoms_of_interest:
                        hb_layer_id += n_ij
                        if donor_i == 0:
                            hb_surface_id += n_ij
        total_hb.append(hb_layer_id)
        surface_hb.append(hb_surface_id)

    total_hb = np.array(total_hb) / n_molecules
    surface_hb = np.array(surface_hb) / n_molecules
    return total_hb, surface_hb, n_molecules


# ---- load trajectory ---------------------------------------------------
u = mda.Universe(tpr, xtc)
u.transfer_to_memory(step=step)
n_frames = u.trajectory.n_frames

liquids, liquid_names, solids, solid_names = detect_groups(u)
liquid_ref = liquids[liquid_names.index("Liq")]

# ---- HB-per-layer analysis for each liquid species (skip Liq/Na/Cl) -----
for group, name in zip(liquids, liquid_names):
    if name in ("Liq", "Na", "Cl") or group.n_atoms == 0:
        continue

    total_hb, surface_hb, n_molecules = hbonds_per_layer(
        u, liquid_ref, name, group.n_residues, n_frames
    )

    np.savetxt(f"hb_layers_{name}.dat",
               np.vstack([total_hb, surface_hb, n_molecules]).T,
               header="total_hb surface_hb n_molecules")