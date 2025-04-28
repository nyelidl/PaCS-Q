import subprocess




def sander_run(crd, top, i, j):
    command = [
        'sander',
        '-O',
        '-i', '../../../qmmm.in',
        '-o', 'qmmm.out',
        '-ref', f'../../../{crd}',
        '-c', f'../../../{crd}',
        '-p', f'../../../{top}',
        '-r', f'qmmm{i}_{j}.rst',
        '-x', f'qmmm{i}_{j}.nc'
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



def sander_run_cyc(crd, top, i, j, location, foldername, ref):
    command = [
        'sander',
        '-O',
        '-i', '../../../qmmm.in',
        '-o', f'qmmm{i}_{j}.out',
        '-ref', f'{location}/{foldername}/{i-1}/{ref}/min.rst',
        '-c', f'{location}/{foldername}/{i-1}/{ref}/min.rst',
        '-p', f'../../../{top}',
        '-r', f'qmmm{i}_{j}.rst',
        '-x', f'qmmm{i}_{j}.nc'
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

#sander_run_mpi_cyc("test", "test", 1, 3, "/usr/location", "MDrun", 3)
