# PaCS-Q
PaCS-Q is a Python toolkit designed to assist with Parallel Cascade Selection simulations (PaCS) for studying protein structural transitions in MD and QM/MM MD levels.

Welcome to PaCS-Q v1.0.5 by L.Duan 2025.4.15

**PaCS-Q** is a Python toolkit designed to assist with
Parallel Cascade Selection simulations (PaCS) for studying
protein structural transitions in MD and QM/MM MD level.

Please cite paper:
-Lian Duan, Kowit Hengphasatporn, Ryuhei Harada, Yasuteru Shigeta. JCTC https://doi.org/10.1021/acs.jctc.5c00169
-Lian Duan, Kowit Hengphasatporn, Yasuteru Shigeta. JCTC https://doi.org/xx.xxx/acs.jcim.xxxxx

                    ██████╗░░█████╗░░█████╗░░██████╗░░░░░░░░░██████╗░
                    ██╔══██╗██╔══██╗██╔══██╗██╔════╝░░░░░░░░██╔═══██╗
                    ██████╔╝███████║██║░░╚═╝╚█████╗░░█████╗║██╗██░██║
                    ██╔═══╝░██╔══██║██║░░██╗░╚═══██╗░╚════╝░╚██████╔╝
                    ██║░░░░░██║░░██║╚█████╔╝██████╔╝░░░░░░░░░╚═██╔═╝░
                    ╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░░░░░░░░░░░░╚═╝░░░

## Prerequisites
#Before installing and running PaCS-Q, ensure you have the following:

- **AMBER** (version 22 or later) compiled with **MPI** support
- **Miniconda** (for environment management)
- Python packages:
  - [MDAnalysis](https://www.mdanalysis.org/)
  - [Biopython](https://biopython.org/)

## Installation
# 1. Clone or Download the Toolkit
# Using the following command:

```bash
wget https://github.com/DLCCS/PaCS-Q/archive/refs/heads/main.zip
unzip main.zip
mv PaCS-Q-main PaCS-Q
```
Or clone it using Git:
```bash
git clone https://github.com/DLCCS/PaCS-Q.git
```

# 2. Export the Path in ~/.bashrc
export PATH="/path/to/PaCS-Q:$PATH"
export PATH="/path/to/PaCS-Q/pacsq_toolkit:$PATH"
```bash
source ~/.bashrc
```

# 3. Create and Activate Conda Environment (recommend)
```bash
conda create -n pacs-q
conda init
conda activate pacs-q
```

# 4. Install Required Python Packages
```bash
pip install --upgrade MDAnalysis
pip install biopython
```

# 5. Input Files
- For LB-PaCS-MD (Distance-based PaCS-Q)
        Topology file from (.top) tLEaP and the coordinate file (.rst or .crd) obtained from after heating up step.
        You can adjust any MD parameter in "md.in".
  
- For QM/MM MD (RMSD-based PaCS-Q)
         In this simulation we need to prepare initial and reference structures, as follows.
         Topology file from tLEaP (.top), the coordinate file (.rst or .crd) obtained from after heating up step, and reference structure in PDB file format (.pdb).
         You can adjust any MD parameter in "qmmm.in".

# 6. Extend and Generate Trajectory
- To extend the simulation, please use "--rerun" in the command line.
- To generate the trajectory and the lastframe in pdb file, please use "cpp.sh" and "pdb_last.sh".
- To clean all of MD directory, please use "clean.sh".


# 7. MD Analysis
- Generate CV and 2D-PES by "pacsana_dis_collection.py"
- PCA analysis by "pacsana_procupine.py"
- Create QM input from PaCS-Q Trajectory "pacsana_QM_input.py"


## For More Detail of PaCS-Q Usage:
PaCS-Q -h
                      

  
