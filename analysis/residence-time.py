#!/usr/bin/env python
# coding: utf-8
"""
For water and ethanol in the first ITIM interfacial layer:
  1) compute the continuous survival probability of the "in-layer" (A) and
     "out-of-layer" (B) states, their characteristic desorption/adsorption
     times (tauA, tauB), and the equilibrium in-layer occupancy (etaA);
  2) compute the occupancy autocorrelation function C(t);
  3) compute the lateral (x,y) MSD of molecules during each continuous
     residence period spent in the first layer.
"""

import numpy as np
import MDAnalysis as mda
import pytim

# ---- inputs ----------------------------------------------------------
tpr = "path-to-tpr"
xtc = "path-to-xtc"
min_segment_length = 5  # frames, minimum residence-segment length kept for MSD


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


def rle(inarray):
    """ 
    Run length encoding. 
    Multi datatype arrays catered for including non Numpy
    returns: tuple (runlengths, startpositions, values) 
    """
    ia = np.asarray(inarray)  # force numpy
    n = len(ia)
    if n == 0:
        return (None, None, None)
    else:
        y = np.array(ia[1:] != ia[:-1]) # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1) # must include last element posi
        z = np.diff(np.append(-1, i)) # run lengths
        p = np.cumsum(np.append(0, z))[:-1] # positions
        return (z, p, ia[i])


def calculate_first_passage_proba(It, dt):
    """Calculate survival probability"""
    # initialisation 
    PsiA = np.zeros((len(It)))
    PsiB = np.zeros((len(It)))
    tauA = 0
    tauB = 0
    # loop over all water molecules
    for i, moltraj in enumerate(It.T):
        # individual first passage times to be averaged over
        thisPsiA = np.zeros(len(It))
        thisPsiB = np.zeros(len(It))
        thisTauA = 0
        thisTauB = 0
        # T will contain the lifetimes in domain 'which'
        [T, _, which] = rle(moltraj)
        # PhiA
        if len(T) > 1: # exclude molecule that did not exchange with the interface (should be irrellvant for long trajectory)
            # integrate lifetimes to get first passage times
            for Tmax in T[which == True]:
                thisPsiA[:Tmax] += np.ones(Tmax)
                thisTauA += Tmax
            # divide by the number of events considered in trajectory
            thisPsiA /= len(T[which == True])
            thisTauA /= len(T[which == True])
            # add to average
            PsiA += thisPsiA
            tauA += thisTauA
            # PhiB (same procedure as PsiA, with False <-> True)
            for Tmax in T[which == False]:
                thisPsiB[:Tmax] += np.ones(Tmax)
                thisTauB += Tmax
            thisPsiB /= len(T[which == False])
            thisTauB /= len(T[which == False])
            PsiB += thisPsiB
            tauB += thisTauB
    for Psi in [PsiA, PsiB]:
        # divide by the number of molecules
        Psi /= len(It.T)
        # finally, normalize \int_0^\infty \mathrm{d}t \Psi(t) = 1
        Psi /= np.sum(Psi) * dt
    tauA /= len(It.T)
    tauB /= len(It.T)
    tauA *= dt
    tauB *= dt
    etaA = tauA/(tauA + tauB)
    return PsiA, PsiB, tauA, tauB, etaA


def calculate_correlation(a,b=None,subtract_mean=False, normed=False):
    meana = int(subtract_mean)*np.mean(a,axis=0)
    if len(a.shape) > 1:
        a2 = np.append(a-meana, np.zeros((2**int(np.ceil((np.log(len(a))/np.log(2))))-len(a),a.shape[1])),axis=0)
    else:
        a2 = np.append(a-meana, np.zeros(2**int(np.ceil((np.log(len(a))/np.log(2))))-len(a)),axis=0)
    data_a = np.append(a2, np.zeros(a2.shape),axis=0)
    fra = np.fft.fft(data_a,axis=0)
    
    if b is None:
        sf = np.conj(fra)*fra
    else:
        meanb = int(subtract_mean)*np.mean(b,axis=0)
        if len(b.shape) > 1:
            b2 = np.append(b-meanb, np.zeros((2**int(np.ceil((np.log(len(b))/np.log(2))))-len(b),b.shape[1])),axis=0)
        else:
            b2 = np.append(b-meanb, np.zeros(2**int(np.ceil((np.log(len(b))/np.log(2))))-len(b)),axis=0)
        data_b = np.append(b2, np.zeros(b2.shape),axis=0)
        frb = np.fft.fft(data_b,axis=0)
        sf = np.conj(fra)*frb
    res = np.fft.ifft(sf,axis=0)
    
    if len(a.shape) > 1:
        # average over all particles/molecules
        cor = (np.real(res[:len(a)])/np.array(range(len(a),0,-1))[:,np.newaxis]).mean(axis=1)
    else:
        cor = np.real(res[:len(a)])/np.array(range(len(a),0,-1))
    if normed:
        return cor/cor[0]
    else:
        return cor


def extract_interface_segments(inout, positions, min_length=2):
    """Continuous periods during which each molecule stays in the layer."""
    n_frames, n_molecules = inout.shape
    segments = []

    for mol in range(n_molecules):
        state = inout[:, mol]
        changes = np.diff(np.r_[False, state, False].astype(int))
        starts = np.where(changes == 1)[0]
        stops = np.where(changes == -1)[0]

        for start, stop in zip(starts, stops):
            if stop - start >= min_length:
                segments.append({
                    "molecule": mol,
                    "start": start,
                    "stop": stop,
                    "positions": positions[start:stop, mol].copy(),
                })

    return segments


def calculate_msd(position, box):
    """Lateral (x,y) MSD with periodic boundary correction."""
    Lx, Ly = box[0], box[1]
    n = len(position)
    msd = np.zeros(n)

    for lag in range(n):
        displacement = position[lag:] - position[: n - lag]
        displacement[:, 0] -= Lx * np.rint(displacement[:, 0] / Lx)
        displacement[:, 1] -= Ly * np.rint(displacement[:, 1] / Ly)
        msd[lag] = np.mean(displacement[:, 0] ** 2 + displacement[:, 1] ** 2)

    return msd


def average_msds(all_msds):
    max_length = max(len(msd) for msd in all_msds)
    msd_sum = np.zeros(max_length)
    msd_count = np.zeros(max_length)

    for msd in all_msds:
        n = len(msd)
        msd_sum[:n] += msd
        msd_count[:n] += 1

    msd_average = np.divide(msd_sum, msd_count,
                             out=np.full(max_length, np.nan), where=msd_count > 0)
    return msd_average, msd_count


# ---- load trajectory ---------------------------------------------------
u = mda.Universe(tpr, xtc)
dt = np.round(u.trajectory.dt, 2)
n_frames = u.trajectory.n_frames

liquids, liquid_names, solids, solid_names = detect_groups(u)
liquid_ref = liquids[liquid_names.index("Liq")]

# only water and ethanol are analyzed; name_atom selects the O used to
# track in/out-of-layer membership for that species (OW for water, Oh for
# the ethanol hydroxyl oxygen)
species = [("Sol", "name OW"), ("Eth", "name Oh")]

# ---- desorption time + lateral MSD, per species --------------------------
for name, name_atom in species:
    group = liquids[liquid_names.index(name)]
    if group.n_atoms == 0:
        continue

    n_res = group.n_residues
    inout = np.zeros((n_frames, n_res))
    positions = np.zeros((n_frames, n_res, 3))

    oxygen = u.select_atoms(f"resname {name} and name O*")
    id_subject = group.select_atoms(name_atom).atoms.ids

    for ts in u.trajectory:
        interface = pytim.ITIM(u, group=liquid_ref, max_layers=4, molecular=True)
        first_layer = (interface.layers[0].tolist()[0]
                       + interface.layers[1].tolist()[0]).select_atoms(name_atom)

        condition = [aid in first_layer.ids for aid in id_subject]
        assert np.sum(condition) == len(first_layer)
        inout[ts.frame] = condition
        positions[ts.frame] = oxygen.positions

    # --- desorption time / survival probability ---
    PsiA, PsiB, tauA, tauB, etaA = calculate_first_passage_proba(inout, dt)
    C = calculate_correlation(inout) / etaA
    t = np.arange(len(C)) * dt

    np.savetxt(f"survival_{name}.dat", np.column_stack([t, C, PsiA, PsiB]),
               header="t C PsiA PsiB")
    np.savetxt(f"survival_{name}_tau.dat", [[tauA, tauB, etaA]],
               header="tauA tauB etaA")

    # --- lateral MSD within first-layer residence segments ---
    segments = extract_interface_segments(inout, positions, min_length=min_segment_length)
    all_msds = [calculate_msd(seg["positions"], u.dimensions) for seg in segments]
    msd, counts = average_msds(all_msds)
    time = np.arange(len(msd)) * dt

    np.savetxt(f"lateral_msd_{name}.dat", np.column_stack([time, msd]),
               header="t msd")