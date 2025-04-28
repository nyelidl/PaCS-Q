import subprocess




def pmemd_run(crd, top, i, j):
    command = [
        'pmemd.cuda',
        '-O',
        '-i', '../../../md.in',
        '-o', 'md.out',
        '-ref', f'../../../{crd}',
        '-c', f'../../../{crd}',
        '-p', f'../../../{top}',
        '-r', f'md{i}_{j}.rst',
        '-x', f'md{i}_{j}.nc',
        '-v', 'mdvel'
    ]
    #print(command)


    with open('run.log', 'w') as stderr_file:
        try:
            result = subprocess.run(
                command,
                stderr=stderr_file,
                text=True
            )
            result.check_returncode()
        except subprocess.CalledProcessError as e:
            print(f"Error，code：{e.returncode}")
            print("Please check 'run.log' for detail")
        except Exception as e:
            print(f"Error：{e}")
            print("Please check 'run.log' for detail")



def pmemd_run_cyc(crd, top, i, j, location, foldername, ref):
    command = [
        'pmemd.cuda',
        '-O',
        '-i', '../../../md.in',
        '-o', f'md{i}_{j}.out',
        '-ref', f'{location}/{foldername}/{i-1}/{ref}/min.rst',
        '-c', f'{location}/{foldername}/{i-1}/{ref}/min.rst',
        '-p', f'../../../{top}',
        '-r', f'md{i}_{j}.rst',
        '-x', f'md{i}_{j}.nc',
        '-v', 'mdvel'
    ]
    #print(command)


    with open('run.log', 'w') as stderr_file:
        try:
            result = subprocess.run(
                command,
                stderr=stderr_file,
                text=True
            )
            result.check_returncode()
        except subprocess.CalledProcessError as e:
            print(f"Error，code：{e.returncode}")
            print("Please check 'run.log' for detail")
        except Exception as e:
            print(f"Error：{e}")
            print("Please check 'run.log' for detail")

