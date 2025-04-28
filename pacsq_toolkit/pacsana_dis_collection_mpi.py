#!/usr/bin/python
import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
import os
from tqdm import tqdm
import argparse
from multiprocessing import Pool
import re


parser = argparse.ArgumentParser(description="""PaCS-Q-Analysis v1.0.0 by L.Duan 2024.11.27""", epilog="""
example: pacsana_dis_collection_mpi.py -n ./MDrun/1/ -t . -x "resid 349 and name C" -y "resid 349 and name C5" -s dis1.dat -c 10
Please cite paper: 
    XXX. et al. JXXX XXX:XXX 2024 DOI:XXXXXXX""", formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-n','--nc', type=str, help='traj file location')
parser.add_argument('-t','--top', type=str, help='top file location')
parser.add_argument('-x','--sel1', type=str, help='traj file location')
parser.add_argument('-y','--sel2', type=str, help='traj file location')
parser.add_argument('-s','--save', type=str, help='name of saved data file')
parser.add_argument('-c','--core', type=int, help='how many core to calculation')

args = parser.parse_args()

location = os.getcwd()
os.chdir(location)

def find_nc_files(directory='.'):
    print("Searching for nc file...")
    file_list = []
    for root, dirs, files in tqdm(os.walk(directory)):
        for file in files:
            if file.endswith('.nc'):
                file_list.append(os.path.join(root, file))
    return file_list
    print("Done!")

def find_top_files(directory='.'):
    print("Searching for top file...")
    file_list = []
    for root, dirs, files in tqdm(os.walk(directory)):
        for file in files:
            if file.endswith('.top'):
                file_list.append(os.path.join(root, file))
    return file_list
    print("Done!")


def cal_dis_all(sel_1="para 1", sel_2="para 2", nc_location=".", top_location="."):
    nc_files = find_nc_files(nc_location)
    sorted_nc_files = sorted(nc_files, key=extract_sort_keys)
    print(sorted_nc_files)
    top_files = find_top_files(top_location)
    sel_1 = sel_1
    sel_2 = sel_2
    d_list = []
    print(f"calculation the reaction {sel_1} & {sel_2} ...")
    for i in tqdm(range(len(nc_files))):
        u = mda.Universe(top_files[0], nc_files[i])
        atom1 = u.select_atoms(sel_1)
        atom2 = u.select_atoms(sel_2)
        dis_traj = []
        for ts in u.trajectory:
            dis_ts = dist(atom1, atom2)
            dis_traj.append(dis_ts[2][0])
        d_list = d_list + dis_traj
    return d_list


def extract_sort_keys(filename):
    """
    从文件路径中提取用于排序的数字。
    返回一个 (num1, num2) 的元组，供排序使用。
    """
    # 使用正则表达式提取数字
    match = re.search(r'/(\d+)/(\d+)/qmmm(\d+)_(\d+)\.nc$', filename)
    if match:
        num1 = int(match.group(1))
        num2 = int(match.group(2))
        # 验证文件名中的数字是否与路径中的数字一致
        if num1 == int(match.group(3)) and num2 == int(match.group(4)):
            return (num1, num2)
            print("Done!")
        else:
            # 如果不一致，可能是格式不同，返回一个大数，使其排在后面
            return (float('inf'), float('inf'))
    else:
        # 如果匹配失败，返回一个大数，使其排在后面
        return (float('inf'), float('inf'))


def process_file(args):
    nc_file, topology, sel_1, sel_2 = args
    u = mda.Universe(topology, nc_file)
    atom1 = u.select_atoms(sel_1)
    atom2 = u.select_atoms(sel_2)
    dis = []
    for ts in u.trajectory:
        dis_ts = dist(atom1, atom2)
        dis.append(dis_ts[2][0])
    return dis


def cal_dis_mpi(sel_1="para 1", sel_2="para 2", nc_location=".", top_location=".", num_cores=1):
    nc_files = find_nc_files(nc_location)
    # 对 nc_files 进行排序
    sorted_nc_files = sorted(nc_files, key=extract_sort_keys)
    top_files = find_top_files(top_location)
    topology = top_files[0]  # 假设使用第一个拓扑文件

    # 准备多处理的参数列表
    args_list = [(nc_file, topology, sel_1, sel_2) for nc_file in sorted_nc_files]

    dis_list = []

    # 使用多处理池，并指定进程数量
    with Pool(processes=num_cores) as pool:
        # 使用 tqdm 显示进度条
        for dis in tqdm(pool.imap_unordered(process_file, args_list), total=len(args_list)):
            dis_list = dis_list + dis

    return dis_list


list_test = cal_dis_mpi(args.sel1, args.sel2, args.nc, args.top, args.core)

with open(args.save, 'w') as file:
    for i in list_test:
        file.write(str(i))
        file.write('\n')


