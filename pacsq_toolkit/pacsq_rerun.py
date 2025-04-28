from pacsq_toolkit.pacsrmsd import min_rmsd_single
from pacsq_toolkit.qmmm_setting import (qmmm_tem_s, qmmm_tem, write_to_sh, write_to_dat, read_dat)
from pacsq_toolkit.sander_run_mpi import sander_run_mpi, sander_run_mpi_cyc
import os
from tqdm import tqdm



def get_latest_folder_name(directory):
    folders = [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]

    numeric_folders = [int(f) for f in folders if f.isdigit()]

    if numeric_folders:
        latest_folder = str(max(numeric_folders))
        return latest_folder
    else:
        return None



def pacsq_rerun(cyc_plus, rep, foldername="MDrun", location=os.getcwd(), crd="qmmm.crd",
                top="qmmm.top", ref_location="F.pdb", selection="./orca/F.pdb", qmmm_int=None):
    if qmmm_int is None:
        write_to_dat("qmmm.in", qmmm_tem)
    else:
        qm_int = read_dat(qmmm_int)
        write_to_dat("qmmm.in", qm_int)
    input = f"parm {location}/{top}"
    print("### Pacs-Q rerunning... ###")
    died_cyc = int(get_latest_folder_name(f"./{foldername}"))
    for i in tqdm(range(died_cyc, cyc_plus + 2)):
        os.chdir(location)
        os.system(f"mkdir ./{foldername}/{i}")
        if i == 0:
            for j in range(1, rep + 1):
                #print(f"rep {j}")
                os.system(f"mkdir ./{foldername}/{i}/{j}")
                os.chdir(f"{location}/{foldername}/{i}/{j}")
                sander_run_mpi(crd, top, i, j)
                os.chdir(location)
        else:
            c_top = f"{location}/{top}"
            dis_list = []
            index_list = []
            #print(f"Calculating Cycle: {i-1}")
            for j in range(1, rep + 1):
                c_nc = f"{location}/{foldername}/{i-1}/{j}/qmmm{i-1}_{j}.nc"
                c_dis, c_index = min_rmsd_single(c_top, c_nc, ref_location, selection)
                dis_list.append(c_dis)
                index_list.append(c_index)
            min_dis = min(dis_list)
            min_index = dis_list.index(min_dis)
            min_index_label = index_list[min_index]
            ref = min_index + 1
            cpp_text = f"\ntrajin {location}/{foldername}/{i-1}/{ref}/qmmm{i-1}_{ref}.nc 1 {min_index_label + 1} 1"
            with open(f"sum-all.sh", "a") as file:
                file.write(cpp_text)
            dis = f"\nUsing Rep: {ref}\nDis: {min_dis}\n"
            with open(f"dis.dat", "a") as file:
                file.write(f"\ncyc {i}\n")
                file.write(dis)
            with open(f"dis_plot.dat", "a") as file:
                file.write(f"{min_dis}\n")
            os.chdir(f"{location}/{foldername}/{i-1}/{ref}")
            write_to_sh("new.sh", f"qmmm{i-1}_{ref}", min_index_label, top)
            os.system(f"cpptraj new.sh > cpp.log")
            os.chdir(location)

            for j in range(1, rep + 1):
                #print(f"running rep {j}")
                os.system(f"mkdir ./{foldername}/{i}/{j}")
                os.chdir(f"{location}/{foldername}/{i}/{j}")
                sander_run_mpi_cyc(crd, top, i, j, location, foldername, ref)
                os.chdir(location)



    with open("sum-all.sh", "r", encoding="utf-8") as f:
        lines = f.readlines()

    seen = set()
    unique_lines = []

    for line in lines:
        if line not in seen:
            unique_lines.append(line)
            seen.add(line)

    with open("sum-all.sh", "w", encoding="utf-8") as f:
        f.writelines(unique_lines)


    with open("dis.dat", "r", encoding="utf-8") as f:
        lines = f.readlines()

    seen = set()
    unique_lines = []

    for line in lines:
        if line not in seen:
            unique_lines.append(line)
            seen.add(line)

    with open("dis.dat", "w", encoding="utf-8") as f:
        f.writelines(unique_lines)



    with open("dis_plot.dat", "r", encoding="utf-8") as f:
        lines = f.readlines()

    seen = set()
    unique_lines = []

    for line in lines:
        if line not in seen:
            unique_lines.append(line)
            seen.add(line)

    with open("dis_plot.dat", "w", encoding="utf-8") as f:
        f.writelines(unique_lines)
