#!/usr/bin/python
import os
import re
import argparse
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm
from MDAnalysis.analysis import distances

try:
    from MDAnalysis.analysis.sasa import ShrakeRupley
except ImportError:
    ShrakeRupley = None


parser = argparse.ArgumentParser(
    description="LB-PaCS-MD selection analysis: ligand SASA + contact number + VFD2 opening angle"
)

parser.add_argument("-n", "--nc", type=str, required=True, help="trajectory directory")
parser.add_argument("-t", "--top", type=str, required=True, help="topology directory")
parser.add_argument("-l", "--ligand", type=str, required=True, help="ligand selection")
parser.add_argument("-p", "--protein", type=str, default="protein", help="protein selection")

parser.add_argument("--lobe1", type=str, required=True, help="VFD2 lobe 1 selection")
parser.add_argument("--hinge", type=str, required=True, help="VFD2 hinge selection")
parser.add_argument("--lobe2", type=str, required=True, help="VFD2 lobe 2 selection")

parser.add_argument("--cutoff", type=float, default=4.0, help="contact cutoff in Angstrom")
parser.add_argument("-s", "--save", type=str, default="lbpacs_score.dat", help="output file")

args = parser.parse_args()


def find_files(directory=".", suffix=".nc"):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(suffix):
                file_list.append(os.path.join(root, file))
    return file_list


def extract_sort_keys(filename):
    """
    Sort files like:
    /1/1/qmmm1_1.nc
    /1/2/qmmm1_2.nc
    /2/1/qmmm2_1.nc
    """
    match = re.search(r"/(\d+)/(\d+)/qmmm(\d+)_(\d+)\.nc$", filename)
    if match:
        n1 = int(match.group(1))
        n2 = int(match.group(2))
        n3 = int(match.group(3))
        n4 = int(match.group(4))
        if n1 == n3 and n2 == n4:
            return (n1, n2)
    return (float("inf"), float("inf"))


def normalize(values):
    values = np.array(values, dtype=float)
    vmin = np.min(values)
    vmax = np.max(values)
    if np.isclose(vmax, vmin):
        return np.zeros_like(values)
    return (values - vmin) / (vmax - vmin)


def calc_angle(p1, p2, p3):
    """
    angle p1-p2-p3 in degrees
    p2 is hinge point
    """
    v1 = p1 - p2
    v2 = p3 - p2

    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)

    if norm1 == 0 or norm2 == 0:
        return np.nan

    cosang = np.dot(v1, v2) / (norm1 * norm2)
    cosang = np.clip(cosang, -1.0, 1.0)

    return np.degrees(np.arccos(cosang))


def calc_ligand_sasa(u, ligand_sel):
    """
    Calculate ligand SASA for current frame.
    Requires MDAnalysis ShrakeRupley.
    """
    if ShrakeRupley is None:
        raise ImportError(
            "MDAnalysis.analysis.sasa.ShrakeRupley is not available. "
            "Please update MDAnalysis or use FreeSASA-based calculation."
        )

    ligand = u.select_atoms(ligand_sel)

    sr = ShrakeRupley(ligand, probe_radius=1.4, n_sphere_points=960)
    sr.run(start=u.trajectory.frame, stop=u.trajectory.frame + 1)

    return np.sum(ligand.atoms.tempfactors)


def calc_contact_number(ligand, protein, cutoff=4.0):
    """
    Residue-based ligand-protein contact number.
    Counts protein residues having at least one atom within cutoff of ligand.
    """
    dist_mat = distances.distance_array(ligand.positions, protein.positions)
    contact_atom_indices = np.where(dist_mat < cutoff)[1]

    if len(contact_atom_indices) == 0:
        return 0

    contacted_atoms = protein[contact_atom_indices]
    contacted_resids = set(contacted_atoms.resids)

    return len(contacted_resids)


def analyze_all(nc_location, top_location):
    nc_files = find_files(nc_location, ".nc")
    nc_files = sorted(nc_files, key=extract_sort_keys)

    top_files = find_files(top_location, ".top")

    if len(nc_files) == 0:
        raise FileNotFoundError("No .nc trajectory files found.")

    if len(top_files) == 0:
        raise FileNotFoundError("No .top topology files found.")

    top_file = top_files[0]

    print(f"Topology: {top_file}")
    print(f"Number of trajectories: {len(nc_files)}")

    all_sasa = []
    all_contacts = []
    all_angles = []
    all_traj = []
    all_frame = []

    for nc in tqdm(nc_files, desc="Analyzing trajectories"):
        u = mda.Universe(top_file, nc)

        ligand = u.select_atoms(args.ligand)
        protein = u.select_atoms(args.protein)
        lobe1 = u.select_atoms(args.lobe1)
        hinge = u.select_atoms(args.hinge)
        lobe2 = u.select_atoms(args.lobe2)

        if len(ligand) == 0:
            raise ValueError(f"Ligand selection returned 0 atoms: {args.ligand}")
        if len(protein) == 0:
            raise ValueError(f"Protein selection returned 0 atoms: {args.protein}")
        if len(lobe1) == 0:
            raise ValueError(f"Lobe1 selection returned 0 atoms: {args.lobe1}")
        if len(hinge) == 0:
            raise ValueError(f"Hinge selection returned 0 atoms: {args.hinge}")
        if len(lobe2) == 0:
            raise ValueError(f"Lobe2 selection returned 0 atoms: {args.lobe2}")

        for ts in u.trajectory:
            sasa = calc_ligand_sasa(u, args.ligand)

            contacts = calc_contact_number(
                ligand=ligand,
                protein=protein,
                cutoff=args.cutoff
            )

            p1 = lobe1.center_of_mass()
            p2 = hinge.center_of_mass()
            p3 = lobe2.center_of_mass()

            angle = calc_angle(p1, p2, p3)

            all_sasa.append(sasa)
            all_contacts.append(contacts)
            all_angles.append(angle)
            all_traj.append(nc)
            all_frame.append(ts.frame)

    sasa_norm = normalize(all_sasa)
    contact_norm = normalize(all_contacts)
    angle_norm = normalize(all_angles)

    score = sasa_norm - contact_norm + angle_norm

    return all_traj, all_frame, all_sasa, all_contacts, all_angles, sasa_norm, contact_norm, angle_norm, score


results = analyze_all(args.nc, args.top)

(
    trajs,
    frames,
    sasa,
    contacts,
    angles,
    sasa_norm,
    contact_norm,
    angle_norm,
    scores,
) = results


with open(args.save, "w") as f:
    f.write(
        "traj\tframe\tligand_SASA\tcontact_number\tVFD2_opening_angle\t"
        "SASA_norm\tContact_norm\tAngle_norm\tScore\n"
    )

    for i in range(len(scores)):
        f.write(
            f"{trajs[i]}\t{frames[i]}\t"
            f"{sasa[i]:.6f}\t{contacts[i]}\t{angles[i]:.6f}\t"
            f"{sasa_norm[i]:.6f}\t{contact_norm[i]:.6f}\t"
            f"{angle_norm[i]:.6f}\t{scores[i]:.6f}\n"
        )

print(f"Done. Output saved to {args.save}")
