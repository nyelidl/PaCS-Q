#!/usr/bin/python
import warnings
from Bio import BiopythonDeprecationWarning
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
import numpy as np
import os
import MDAnalysis as mda
from MDAnalysis.analysis import pca, align
import argparse


parser = argparse.ArgumentParser(description="""PaCS-Q-Analysis v1.0.0 by L.Duan 2024.11.27""", epilog="""
example: pacsana_procupine.py -n sum.nc -t complex.top -x "resname CHA" -s save -b 0 -l 100

Please cite paper: 
    XXX. et al. JXXX XXX:XXX 2024 DOI:XXXXXXX""", formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-n','--nc', type=str, help='traj file name')
parser.add_argument('-t','--top', type=str, help='top file name')
parser.add_argument('-x','--sel', type=str, help='atoms selection')
parser.add_argument('-s','--save', type=str, help='name of saved data file')
parser.add_argument('-b','--str', type=int, help='begin frame')
parser.add_argument('-l','--stp', type=int, help='last frame')

args = parser.parse_args()

location = os.getcwd()
os.chdir(location)


u = mda.Universe(args.top, args.nc)
start_point = args.str
stop_point = args.stp
select = args.sel
selection = u.select_atoms(select)
file_name = args.save

pc = pca.PCA(u, select=select,
             align=False, mean=None,
             n_components=None).run(start=start_point, stop=stop_point)

modes = pc.results.p_components
variances = pc.results.variance

frame_index = args.str
u.trajectory[frame_index]
atoms = u.atoms.select_atoms(args.sel)
coordinates = atoms.positions.flatten()  # 原始坐标数据
atomnames = [atom.name for atom in atoms]
resnames = [atom.resname for atom in atoms]
resids = [atom.resid for atom in atoms]
chainids = [atom.segid for atom in atoms]
bfactors = np.zeros(len(atoms))


with open(f'./{file_name}_result.nmd', 'w') as file:
    file.write(f'pca_result\n')

    file.write('names ')
    file.write(' '.join(atomnames))
    file.write('\n')

    file.write('resnames ')
    file.write(' '.join(resnames))
    file.write('\n')

    file.write('resids ')
    file.write(' '.join(map(str, resids)))
    file.write('\n')

    file.write('coordinates ')
    file.write(' '.join(f'{c:.3f}' for c in coordinates))
    file.write('\n')

    for i, (mode, variance) in enumerate(zip(modes.T, variances)):
        file.write(f'mode {i + 1} {variance:.2f} ')
        file.write(' '.join(f'{m:.3f}' for m in mode))
        file.write('\n')




