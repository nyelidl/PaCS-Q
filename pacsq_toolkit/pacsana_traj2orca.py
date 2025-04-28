#!/usr/bin/python
import warnings
from Bio import BiopythonDeprecationWarning
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
import MDAnalysis as mda
import argparse
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Atom
import warnings


# Suppress specific MDAnalysis PDB warnings
warnings.filterwarnings("ignore", message=r"Found no information for attr: '.*'")
warnings.filterwarnings("ignore", message=r"Found missing chainIDs.*")
warnings.filterwarnings("ignore", message=r"Atom with index >=100000 cannot write bonds to PDB CONECT records.")
warnings.filterwarnings("ignore", message=r"Unit cell dimensions not found. CRYST1 record set to unitary values.")



parser = argparse.ArgumentParser(description="""PaCS-Q-Analysis v1.0.0 by L.Duan 2025.2.19""", epilog="""
example: pacsana_traj2g16.py -n sum.nc -t complex.top -f 500 -s "resid 200 331 444 540" 

Please cite paper: 
    XXX. et al. JXXX XXX:XXX 2024 DOI:XXXXXXX""", formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-n','--nc', type=str, help='traj file name')
parser.add_argument('-t','--top', type=str, help='top file name')
parser.add_argument('-f','--frame', type=int, help='frame')
parser.add_argument('-s','--selection', type=str, help='selection')


args = parser.parse_args()

location = os.getcwd()
os.chdir(location)


u = mda.Universe(args.top, args.nc)
frame_index = args.frame


if frame_index < 0 or frame_index >= len(u.trajectory):
    raise ValueError(
        f"Frame {args.frame} is out of range. "
        f"Trajectory contains {len(u.trajectory)} frames."
    )

u.trajectory[frame_index]
output_filename = f"frame_{args.frame}.pdb"
with mda.Writer(output_filename, u.atoms.n_atoms) as W:
    W.write(u.atoms)

print(f"Frame {args.frame} saved to {output_filename}")

user_selected = args.selection
selected_residue = u.select_atoms(f"{user_selected} and not (name N or name C or name O or name H)")
remaining_atoms = selected_residue
new_u = mda.Merge(remaining_atoms)
new_u.atoms.write("edited.pdb")

def compute_rotation_matrix(ideal_cb, ideal_hstar, actual_cb, actual_hstar):
    """
    Compute a rotation matrix that aligns the two ideal vectors (ideal_cb and ideal_hstar)
    with the corresponding actual directions (actual_cb and actual_hstar).
    """
    # Build ideal local frame:
    u_i = ideal_cb / np.linalg.norm(ideal_cb)
    proj_i = np.dot(ideal_hstar, u_i)
    v_i = ideal_hstar - proj_i * u_i
    if np.linalg.norm(v_i) < 1e-6:
        v_i = np.cross(u_i, np.array([1, 0, 0]))
        if np.linalg.norm(v_i) < 1e-6:
            v_i = np.cross(u_i, np.array([0, 1, 0]))
    v_i = v_i / np.linalg.norm(v_i)
    w_i = np.cross(u_i, v_i)

    # Build actual local frame:
    u_a = actual_cb / np.linalg.norm(actual_cb)
    proj_a = np.dot(actual_hstar, u_a)
    v_a = actual_hstar - proj_a * u_a
    if np.linalg.norm(v_a) < 1e-6:
        v_a = np.cross(u_a, np.array([1, 0, 0]))
        if np.linalg.norm(v_a) < 1e-6:
            v_a = np.cross(u_a, np.array([0, 1, 0]))
    v_a = v_a / np.linalg.norm(v_a)
    w_a = np.cross(u_a, v_a)

    R_ideal = np.column_stack((u_i, v_i, w_i))
    R_actual = np.column_stack((u_a, v_a, w_a))
    R = np.dot(R_actual, R_ideal.T)
    return R

# ---------------- Parameters ----------------
bond_length_CH = 1.09  # typical C–H bond length in Å

# Define four ideal tetrahedral vertices.
# These four unit vectors have the property that any two distinct vertices have a dot product of -1/3.
ideal_vertices = np.array([
    [ 1,  1,  1],
    [ 1, -1, -1],
    [-1,  1, -1],
    [-1, -1,  1]
], dtype=float)
ideal_vertices = ideal_vertices / np.linalg.norm(ideal_vertices[0])

# ---------------- Parse the Input PDB ----------------
parser = PDBParser(QUIET=True)
structure = parser.get_structure("my_structure", "edited.pdb")

# Loop over all residues in all models/chains.
for model in structure:
    for chain in model:
        for residue in chain:
            # Require that the residue has CA and CB.
            if "CA" not in residue or "CB" not in residue:
                continue

            CA = residue["CA"]
            CB = residue["CB"]
            CA_coord = CA.get_coord()

            # -------- Select the shortest H (H*) --------
            hstar = None
            hstar_dist = None
            for atom in residue:
                if atom.element == "H" and atom.get_id() not in ["H1", "H2"]:
                    d = np.linalg.norm(atom.get_coord() - CA_coord)
                    if hstar is None or d < hstar_dist:
                        hstar = atom
                        hstar_dist = d
            if hstar is None:
                # Skip this residue if no hydrogen is found.
                continue

            # Compute the actual (normalized) vectors from CA toward CB and H*
            v_CB = CB.get_coord() - CA_coord
            v_CB /= np.linalg.norm(v_CB)
            v_Hstar = hstar.get_coord() - CA_coord
            v_Hstar /= np.linalg.norm(v_Hstar)

            # -------- Ideal Tetrahedron Matching --------
            # Make a working copy of the ideal vertices.
            verts = ideal_vertices.copy()
            # For the CB substituent: select the ideal vertex that best aligns with v_CB.
            dots_cb = np.dot(verts, v_CB)
            idx_cb = np.argmax(dots_cb)
            ideal_cb = verts[idx_cb]
            verts = np.delete(verts, idx_cb, axis=0)

            # For the H* substituent: select the ideal vertex (from the remaining ones)
            # that best aligns with v_Hstar.
            dots_hstar = np.dot(verts, v_Hstar)
            idx_hstar = np.argmax(dots_hstar)
            ideal_hstar = verts[idx_hstar]
            verts = np.delete(verts, idx_hstar, axis=0)

            # The two remaining vertices will be assigned to the new hydrogens H1 and H2.
            ideal_h1 = verts[0]
            ideal_h2 = verts[1]

            # -------- Compute Rotation Matrix --------
            # Compute a rotation that maps the ideal directions for CB and H* to the actual ones.
            R = compute_rotation_matrix(ideal_cb, ideal_hstar, v_CB, v_Hstar)

            # -------- Compute New Positions for H1 and H2 --------
            new_H1_coord = CA_coord + bond_length_CH * np.dot(R, ideal_h1)
            new_H2_coord = CA_coord + bond_length_CH * np.dot(R, ideal_h2)

            # -------- Create and Add New Hydrogen Atoms --------
            h1 = Atom.Atom(name="H1", coord=new_H1_coord,
                           bfactor=CA.get_bfactor(), occupancy=1.0,
                           altloc=" ", fullname=" H1 ", serial_number=0,
                           element="H")
            h2 = Atom.Atom(name="H2", coord=new_H2_coord,
                           bfactor=CA.get_bfactor(), occupancy=1.0,
                           altloc=" ", fullname=" H2 ", serial_number=0,
                           element="H")
            residue.add(h1)
            print(f"add hydrogen at {new_H1_coord}")
            residue.add(h2)
            print(f"add hydrogen at {new_H2_coord}")

# ---------------- Write Out the Modified Structure ----------------
io = PDBIO()
io.set_structure(structure)
io.save(f"{frame_index}_split.pdb")

# Load the final PDB file
u = mda.Universe(f"{frame_index}_split.pdb")

# Extract atomic coordinates
atoms = []
ca_indices = []  # Store indices of all CA atoms

for i, atom in enumerate(u.atoms):
    symbol = atom.element  # MDAnalysis stores the element
    x, y, z = atom.position
    atoms.append((symbol, x, y, z, atom.name))

    if atom.name == "CA":  # Locate **all** CA atoms
        ca_indices.append(i)  # Store index of each CA atom

# Write the Gaussian input file
with open(f"{frame_index}.inp", "w") as f:
    f.write(f"! B3LYP D3 def2-TZVP def2/J RIJCOSX noautostart miniprint nopop\n")
    f.write(f"%maxcore     6000\n")
    f.write(f"%pal nprocs   24 end\n")
    f.write(f"* xyz   0   1\n")  # Number of processors

    # Write atomic coordinates with proper alignment
    for i, (symbol, x, y, z, name) in enumerate(atoms):
        f.write(f"{symbol:<2}    {x:>11.8f}   {y:>11.8f}   {z:>11.8f}\n")

    f.write(" *\n")  # Blank line at the end for Gaussian

print(f"orca input file saved as '{frame_index}.inp'")









