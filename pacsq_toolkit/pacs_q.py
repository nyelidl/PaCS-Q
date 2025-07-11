#!/usr/bin/python
import argparse
import os
import warnings



def main():
    from pacsq_toolkit.pacsq_run import pacsq_run
    from pacsq_toolkit.pacsq_exq_run import pacsq_exq_run
    from pacsq_toolkit.pacsq_rerun import get_latest_folder_name, pacsq_rerun
    from pacsq_toolkit.file_find import find_top_files, find_nc_files, find_crd_files
    
    crd_file = find_crd_files()
    default_crd = crd_file[0] if crd_file else None
    top_file = find_top_files()
    default_top = top_file[0] if top_file else None


    parser = argparse.ArgumentParser(description="""Welcome to PaCS-Q v1.0.10 by L.Duan 2025.7.11
    
    
    
                    ██████╗░░█████╗░░█████╗░░██████╗░░░░░░░░░██████╗░
                    ██╔══██╗██╔══██╗██╔══██╗██╔════╝░░░░░░░░██╔═══██╗
                    ██████╔╝███████║██║░░╚═╝╚█████╗░░█████╗║██╗██░██║
                    ██╔═══╝░██╔══██║██║░░██╗░╚═══██╗░╚════╝░╚██████╔╝
                    ██║░░░░░██║░░██║╚█████╔╝██████╔╝░░░░░░░░░╚═██╔═╝░
                    ╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░░░░░░░░░░░░╚═╝░░░
                    ///////////////QM/MM MD-SIMULATION///////////////
    """, epilog="""
example: pacs_q -cy 4000 -cd 5 -r ./F.pdb -s "resname CHA" -qm qmmm.in
         pacs_q --rerun -cy 2000 -cd 5 -r ./nonc.pdb -s "resid 73 97 157" -q qmmm.in
         
!!! Warning !!!
    Don't name your files starting with 'dis' or 'sum-all', they will be deleted by clean code!
         
Please cite paper: 
    1. Lian Duan, Kowit Hengphasatporn, Ryuhei Harada, and Yasuteru Shigeta Journal of Chemical Theory and Computation 2025 21 (8), 4309-4318 DOI: 10.1021/acs.jctc.5c00169
    2. Lian Duan, Kowit Hengphasatporn, and Yasuteru Shigeta. Journal of Chemical Information and Modeling 2025 DOI: 10.1021/acs.jcim.5c00936



    """, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-cy', '--cyc', type=int, help='How many cycles to run?')

    parser.add_argument('-cd', '--candi', type=int, help='How many candidates to run?')

    parser.add_argument('-qm', '--qms', type=str, default=None, help=f"""QM/MM MD input file""")

    parser.add_argument('-p','--top', type=str, default=default_top, help=f"""PaCS-Q will be automatically specify your topology file. If you need to specify topology file by yourself, please use this keyword.
    \033[1mdefault: {default_top}\033[0m 
    ***Warning: If you want let program detect your coordinate file automatically, you should name like XXX.top""")

    parser.add_argument('-c','--crd', type=str, default=default_crd, help=f"""PaCS-Q will be automatically specify your rst or crd file. If you need to specify rst or crd file by yourself, please use this keyword.
    \033[1mdefault: {default_crd}\033[0m 
    ***Warning: If you want let program detect your coordinate file automatically, you should name like XXX.rst or XXX.crd""")

    parser.add_argument('-r','--ref', type=str, default=None, help=f"""For RMSD based selection PaCS-Q: specify your reference structure file name in PDB, example: ./ref.pdb""")

    parser.add_argument('-s','--sel', type=str, help=f"""Specify atom or residue for PaCS-Q selection, example: resid 5-7""")

    parser.add_argument('-e', '--exq', type=str, default=None, help=f"""Name of your extend qm input (dev)""")

    parser.add_argument('--rerun', action='store_true',
                        help="This section can rerun your calculation from the died point")

    parser.add_argument('--exqm', action='store_true',
                        help="Run with extend qm software (dev)")

    parser.add_argument('-d','--dir', type=str, default="MDrun", help="""Specify your run directory
    default: MDrun""")

    parser.add_argument('-l','--loc', type=str, default=os.getcwd(), help=f"""Path to PaCS-Q work directory
    default: {os.getcwd()}""")

    args = parser.parse_args()

    # print
    parser.print_help()

    if args.rerun:
        pacsq_rerun(args.cyc, args.candi, args.dir, args.loc, args.crd, args.top, args.ref, args.sel, args.qms)
    elif args.exqm:
        pacsq_exq_run(args.cyc, args.candi, args.dir, args.loc, args.crd, args.top, args.ref, args.sel, args.qms, args.exq)
    else:
        pacsq_run(args.cyc, args.candi, args.dir, args.loc, args.crd, args.top, args.ref, args.sel, args.qms)

if __name__ == "__main__":
    main()



