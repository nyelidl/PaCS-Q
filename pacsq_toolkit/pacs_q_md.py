#!/usr/bin/python
import argparse
import os
import warnings



def main():
    from pacsq_toolkit.pacsq_run import pacsq_run
    from pacsq_toolkit.pacsq_exq_run import pacsq_exq_run
    from pacsq_toolkit.pacsq_pmemd_rerun import pacsq_pmemd_rerun
    from pacsq_toolkit.pacsq_pmemd_run import pacsq_pmemd_run_rmsd, pacsq_pmemd_run_dis
    from pacsq_toolkit.pacsq_rerun import get_latest_folder_name, pacsq_rerun
    from pacsq_toolkit.file_find import find_top_files, find_nc_files, find_crd_files

    crd_file = find_crd_files()
    default_crd = crd_file[0] if crd_file else None
    top_file = find_top_files()
    default_top = top_file[0] if top_file else None

    parser = argparse.ArgumentParser(description="""Welcome to PaCS-Q v1.0.9 by L.Duan 2025.4.28
    
    
    
                    ██████╗░░█████╗░░█████╗░░██████╗░░░░░░░░░██████╗░
                    ██╔══██╗██╔══██╗██╔══██╗██╔════╝░░░░░░░░██╔═══██╗
                    ██████╔╝███████║██║░░╚═╝╚█████╗░░█████╗║██╗██░██║
                    ██╔═══╝░██╔══██║██║░░██╗░╚═══██╗░╚════╝░╚██████╔╝
                    ██║░░░░░██║░░██║╚█████╔╝██████╔╝░░░░░░░░░╚═██╔═╝░
                    ╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░░░░░░░░░░░░╚═╝░░░
                    //////////////////MD-SIMULATION//////////////////
    """, epilog="""

example: 
RMSD based PaCS-Q:
    Mandatory files: Reference structure (ref.pdb), MD input file (md.in), topology (.top) and coordinate (.rst or .crd) files
         pacs_q_md_test.py -cy 100 -cd 5 -r ./ref.pdb -s "resname MOL" -md md.in
         pacs_q_md_test.py --rerun -cy 100 -cd 5 -s "resname MOL" -md md.in
         
Distance based PaCS-Q:
    Mandatory files: MD input file (md.in), topology (.top) and coordinate (.rst or .crd) files 
         pacs_q_md_test.py -cy 100 -cd 5 -s "resid 73" -s2 "resid 150" -md md.in -m b
         pacs_q_md_test.py --rerun -cy 100 -cd 5 -s "resid 73" -s2 "resid 150" -md md.in -m b
         
!!! Warning !!!
    Don't name your files starting with 'dis' or 'sum-all', they will be deleted by clean code!
         
Please cite paper: 
    Lian Duan, Kowit Hengphasatpron, Ryuhei Harada, Yasuteru Shigeta. JCTC https://doi.org/10.1021/acs.jctc.5c00169 

    """, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-cy', '--cyc', type=int, help='How many cycles to run?')

    #parser.add_argument('-r','--rep', type=int, help='How many candidates to run')
    parser.add_argument('-cd', '--candi', type=int, help='How many candidates to run?')

    parser.add_argument('-md', '--mds', type=str, default=None, help=f"""MD input file""")

    parser.add_argument('-p','--top', type=str, default=default_top, help=f"""PaCS-Q will be automatically specify your topology file. If you need to specify topology file by yourself, please use this keyword.
    \033[1mdefault: {default_top}\033[0m 
    ***Warning: If you want let program detect your coordinate file automatically, you should name like XXX.top""")

    parser.add_argument('-c','--crd', type=str, default=default_crd, help=f"""PaCS-Q will be automatically specify your rst or crd file. If you need to specify rst or crd file by yourself, please use this keyword.
    \033[1mdefault: {default_crd}\033[0m 
    ***Warning: If you want let program detect your coordinate file automatically, you should name like XXX.rst or XXX.crd""")

    parser.add_argument('-r','--ref', type=str, default=None, help=f"""For RMSD based selection PaCS-Q: specify your reference structure file name in PDB, example: ./ref.pdb""")

    parser.add_argument('-s','--sel', type=str, help=f"""Specify atom or residue for PaCS-Q selection, example: resid 5-7;
    Specify only this selection for RMSD based selection; Specify -s as the first selection and -s2 as the second selection for Distance based selection""")

    parser.add_argument('-s2','--sel2', type=str, default=None, help=f"""For distance based selection PaCS-Q: Specify atom or residue for the second selection, example: resid 8""")

    parser.add_argument('-m', '--set', type=str, default="b", help=f"""For distance based selection PaCS-Q: type b for binding simulation or u for unbinding simulation""")

    parser.add_argument('--rerun', action='store_true',
                        help="This section can rerun your calculation from the died point")

    parser.add_argument('-d','--dir', type=str, default="MDrun", help="""Specify your run directory
    default: MDrun""")

    parser.add_argument('-l','--loc', type=str, default=os.getcwd(), help=f"""Path to PaCS-Q work directory
    default: {os.getcwd()}""")

    args = parser.parse_args()

    # print
    parser.print_help()

    if args.sel2 is None:
        print("run pmemd by RMSD")
        pacsq_pmemd_run_rmsd(args.cyc, args.candi, args.dir, args.loc, args.crd, args.top, args.ref, args.sel, args.mds)
    elif args.ref is None:
        if args.set == "b":
            if args.rerun:
                print("rerun pmemd by Distance (binding)")
                pacsq_pmemd_rerun(args.cyc, args.candi, args.dir, args.loc, args.crd, args.top, args.sel, args.sel2,
                                    args.mds, 1)
            else:
                print("run pmemd by Distance (binding)")
                pacsq_pmemd_run_dis(args.cyc, args.candi, args.dir, args.loc, args.crd, args.top, args.sel, args.sel2,
                                    args.mds, 1)

        if args.set == "u":
            if args.rerun:
                print("rerun pmemd by Distance (unbinding)")
                pacsq_pmemd_rerun(args.cyc, args.candi, args.dir, args.loc, args.crd, args.top, args.sel, args.sel2,
                                  args.mds, 0)
            else:
                print("run pmemd by Distance (unbinding)")
                pacsq_pmemd_run_dis(args.cyc, args.candi, args.dir, args.loc, args.crd, args.top, args.sel, args.sel2,
                                        args.mds, 0)


    elif args.sel2 or args.ref is None:
        print("should select rmsd or distace")
    #if args.rerun:
    #    pacsq_rerun(args.cyc, args.rep, args.fol, args.loc, args.crd, args.top, args.ref, args.sel, args.qms)
    #elif args.exqm:
    #    pacsq_exq_run(args.cyc, args.rep, args.fol, args.loc, args.crd, args.top, args.ref, args.sel, args.qms, args.exq)
    #else:
    #    pacsq_run(args.cyc, args.rep, args.fol, args.loc, args.crd, args.top, args.ref, args.sel, args.qms)

if __name__ == "__main__":
    main()



