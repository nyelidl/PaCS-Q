import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
from MDAnalysis.analysis import align, rms
import numpy as np


def min_rmsd_single(top, nc, ref_path="./orca/F.pdb", atom_selection="resname CHA"):
    reference = mda.Universe(ref_path)
    u = mda.Universe(top, nc)
    mobile_atoms = u.select_atoms(atom_selection)
    reference_atoms = reference.select_atoms(atom_selection)
    aligner = align.AlignTraj(u, reference, select=atom_selection, in_memory=True)
    aligner.run()
    R = rms.RMSD(mobile_atoms, reference_atoms, select=atom_selection)
    R.run()
    rmsd_data = R.results.rmsd  # R.rmsd is a numpy array with columns: (frame, time, RMSD)
    rmsd_values = rmsd_data[:, 2]  # Extracts just the RMSD column
    d_d = rmsd_values.tolist()
    min_value = min(d_d)
    index = d_d.index(min_value)
    return min_value, index

def min_dis_single(top, nc, atom_selection1, atom_selection2, choose = 1):
    u = mda.Universe(top, nc)
    atom1 = u.select_atoms(atom_selection1)
    atom2 = u.select_atoms(atom_selection2)
    d_d = []
    for ts in u.trajectory:
        center1 = atom1.center_of_mass()
        #print(center1)
        center2 = atom2.center_of_mass()
        #print(center2)
        distance = np.linalg.norm(center1 - center2)
        #print(distance)
        d_d.append(distance)
    if choose == 1:
        min_value = min(d_d)
    if choose == 0:
        min_value = max(d_d)
    index = d_d.index(min_value)
    return min_value, index

"""
def last_dis(top, nc, atom1, atom2):
    u = mda.Universe(top, nc)
    atom1 = u.select_atoms("resid 349 and name C3")  # 第一个残基的 Cα 原子
    atom2 = u.select_atoms("resid 349 and name O2")  # 第十个残基的 Cα 原子
    d_d = []
    for ts in u.trajectory:
        distance = dist(atom1, atom2)
        #d_d.append(distance[2][0].item())
        d_d.append(distance[2][0])
    min_value = min(d_d)
    index = d_d.index(min_value)
    return min_value, index
"""
