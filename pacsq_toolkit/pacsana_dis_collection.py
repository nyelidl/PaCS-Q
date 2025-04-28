#!/usr/bin/python
import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
import os
from tqdm import tqdm
import re
import argparse



parser = argparse.ArgumentParser(description="""PaCS-Q-Analysis v1.0.0 by L.Duan 2024.11.27""", epilog="""
example: pacsana_dis_collection.py -n ./MDrun/1/ -t . -x "resid 349 and name C" -y "resid 349 and name C5" -s dis1.dat
Please cite paper: 
    XXX. et al. JXXX XXX:XXX 2024 DOI:XXXXXXX""", formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-n','--nc', type=str, help='traj file location')
parser.add_argument('-t','--top', type=str, help='top file location')
parser.add_argument('-x','--sel1', type=str, help='traj file location')
parser.add_argument('-y','--sel2', type=str, help='traj file location')
parser.add_argument('-s','--save', type=str, help='name of saved data file')

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

list_test = cal_dis_all(args.sel1, args.sel2, args.nc, args.top)

with open(args.save, 'w') as file:
    for i in list_test:
        file.write(str(i))
        file.write('\n')


