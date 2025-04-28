# qmmm templete
from os import write

qmmm_tem = f"""Short QM/MM simulation
&cntrl
imin=0, nstlim=5,
ntx=1,
irest=0,
dt=0.001,
ntt=3, tempi=300.0, temp0=300.0,
gamma_ln=2,
ntb=2, ntp=1,
cut=8.0,
ntpr=1,
ntwx=1,
ntwr=100,
iwrap=1,
pres0=1.0,
barostat=1,
ig=-1,
ifqnt=1
/
&qmmm
qmmask=':73,157',
qmcharge=0,
qm_theory='DFTB',
qmcut=8.0,
qmshake=0,
writepdb=1,
/
"""


qmmm_tem_s = f"""Short QM/MM simulation
&cntrl
imin=0, nstlim=10,
ntx=1,
irest=0,
dt=0.001,
ntt=3, tempi=300.0, temp0=300.0,
gamma_ln=2,
ntf=2, ntc=2,
ntb=2, ntp=1,
cut=8.0,
ntpr=10,
ntwx=1,
ntwr=100,
iwrap=1,
pres0=1.0,
barostat=1,
ig=-1,
ifqnt=1
/
&qmmm
qmmask=':349',
qmcharge=0,
qm_theory='PM3',
qmcut=8.0,
qmshake=1,
writepdb=1, 
/
"""

md_tem_s = f"""MD simulation
#MD simulation for 100 ps
 &cntrl
  imin=0,
  ntx=1,
  irest=0,
  nstlim=50000,
  dt=0.002,
  ntf=2,
  ntc=2,
  temp0=300.0,
  ntpr=500,
  ntwx=500,
  cut=12.0,
  ntb=2,
  ntp=1,
  ntt=3,
  barostat=1,
  gamma_ln=2.0,
  ig=-1,
  iwrap=1
 /
"""

bo = """
center !:WAT mass origin
image origin center familiar

trajout sum.nc
go
"""

cpp = """#!/usr/bin/bash

cat sum-all.sh bo.dat > sum2.sh

cpptraj sum2.sh
"""

pdb_last = """#!/usr/bin/bash
head -n1 sum-all.sh > hd.dat
cat <<EOF > last.dat
trajin ./sum.nc        lastframe
strip :WAT
trajout last.pdb
go
EOF

cat hd.dat last.dat > pdb.sh

cpptraj pdb.sh
rm hd.dat last.dat
"""

clean = """#!/usr/bin/bash

rm -r MDrun
rm dis*
rm sum-all*
rm bo.dat
rm cpp.sh

"""
qmset = """# QMMM
%pal nprocs 24 end
! BLYP SVP
"""



def write_to_sh(file_name, name="qmmm_0_1", inx=10, top="complex.top"):
    New_input = f"""parm ../../../{top}
trajin {name}.nc {inx + 1} {inx + 1} 1
trajout min.rst
go
"""
    with open(f"{file_name}", "w") as file:
        file.write(New_input)


def write_to_dat(file_name, text):
    with open(f"{file_name}", "w") as file:
        file.write(text)

# Open the file in read mode
def read_dat(file_name):
    with open(f'{file_name}', 'r') as file:
        # Read the contents of the file
        contents = file.read()
    return contents





