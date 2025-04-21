# PaCS-Q (dev)
PaCS-Q is a Python toolkit designed to assist with Parallel Cascade Selection simulations (PaCS) for studying protein structural transitions in MD and QM/MM MD levels.

Welcome to PaCS-Q v1.0.5 by L.Duan 2025.4.15

**PaCS-Q** is a Python toolkit designed to assist with
Parallel Cascade Selection simulations (PaCS) for studying
protein structural transitions in MD and QM/MM MD level.

## Please cite paper:
1. Lian Duan, Kowit Hengphasatporn, Ryuhei Harada, Yasuteru Shigeta. JCTC https://doi.org/10.1021/acs.jctc.5c00169
2. Lian Duan, Kowit Hengphasatporn, Yasuteru Shigeta. JCIM https://doi.org/xx.xxx/acs.jcim.xxxxx

                    ██████╗░░█████╗░░█████╗░░██████╗░░░░░░░░░██████╗░
                    ██╔══██╗██╔══██╗██╔══██╗██╔════╝░░░░░░░░██╔═══██╗
                    ██████╔╝███████║██║░░╚═╝╚█████╗░░█████╗║██╗██░██║
                    ██╔═══╝░██╔══██║██║░░██╗░╚═══██╗░╚════╝░╚██████╔╝
                    ██║░░░░░██║░░██║╚█████╔╝██████╔╝░░░░░░░░░╚═██╔═╝░
                    ╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░░░░░░░░░░░░╚═╝░░░

## Prerequisites
**Before installing and running PaCS-Q, ensure you have the following:**

- **AMBER** (version 22 or later) compiled with **MPI(sander.MPI)** support and **CUDA(pmemd.cuda)** support
- **Miniconda** (for environment management)
- Python packages:
  - [MDAnalysis](https://www.mdanalysis.org/)
  - [Biopython](https://biopython.org/)

## Installation
### 1. Clone or Download the Toolkit
**Using the following command:**

Clone it using Git:
```bash
git clone https://github.com/nyelidl/PaCS-Q.git
```

### 2. Export the Path in ~/.bashrc
export PATH="/path/to/PaCS-Q:$PATH"
export PATH="/path/to/PaCS-Q/pacsq_toolkit:$PATH"
```bash
source ~/.bashrc
```

### 3. Create and Activate Conda Environment (recommend)
```bash
conda create -n pacs-q
conda init
conda activate pacs-q
```

### 4. Install Required Python Packages
```bash
pip install --upgrade MDAnalysis
pip install biopython
```

### 5. Input Files
- For LB-PaCS-MD (Distance-based PaCS-Q)
        Topology file from (.top) tLEaP and the coordinate file (.rst or .crd) obtained from after heating up step.
        You can adjust any MD parameter in "md.in".
  
- For QM/MM MD (RMSD-based PaCS-Q)
         In this simulation we need to prepare initial and reference structures, as follows.
         Topology file from tLEaP (.top), the coordinate file (.rst or .crd) obtained from after heating up step, and reference structure in PDB file format (.pdb).
         You can adjust any MD parameter in "qmmm.in".

### 6. Extend and Generate Trajectory
- To extend the simulation, please use "--rerun" in the command line.
- To generate the trajectory and the lastframe in pdb file, please use "cpp.sh" and "pdb_last.sh".
- To clean all of MD directory, please use "clean.sh".


### 7. MD Analysis
- Generate CV and 2D-PES by "pacsana_dis_collection.py"
- PCA analysis by "pacsana_procupine.py"
- Create QM input from PaCS-Q Trajectory "pacsana_QM_input.py"


# Example
## Case I
### Introduction
The **Grotthuss mechanism** reveals the unique process by which protons move in water: water molecules first absorb an extra proton to form a hydronium ion, and this ion rapidly transfers the proton to a neighboring molecule through hydrogen bond interactions. This results in a series of consecutive proton hops that greatly enhance the efficiency of proton conduction. Today, this model has become one of the key theories for explaining fast proton conduction in water. In this work, I will reconstruct this mechanism using the PaCS-Q method.

### Preparation
As with many sampling techniques, having a good initial configuration is the key to success. Therefore, I strongly recommend that you perform initial sampling through molecular dynamics simulations, selecting a relatively appropriate hydrogen-bond network from the simulation for your study. There are already plenty of tutorials available on molecular dynamics, so I will not take up the space here to describe how molecular dynamics is carried out. Through molecular dynamics, we can obtain a topology file (.top).

As shown in the figure, we selected a short hydrogen-bond network from the molecular dynamics trajectory file, and saved the coordinates of that frame (in a .rst file).

Since we need to study the process of proton transfer, we optimize the structure after the reaction using any QM software (in this example, ORCA 6.0 is used), examine the optimized structure, and save it as a PDB file. **It is important to emphasize that the number of atoms and the order of atoms in each residue in the PDB file must be identical to those in the topology file—even if the proton has transferred and no longer belongs to that residue, the atom must still be included within the original residue range.**

![fig1](case_1/fig1.png)


## For More Detail of PaCS-Q Usage:
pacs_q.py -h
                      

  
