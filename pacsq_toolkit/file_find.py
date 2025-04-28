import os

def find_nc_files(directory='.'):
    file_list = []
    for file in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, file)) and file.endswith('.nc'):
            file_list.append(os.path.join(directory, file))
    return file_list

def find_top_files(directory='.'):
    file_list = []
    for file in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, file)) and file.endswith('.top'):
            file_list.append(os.path.join(directory, file))
    return file_list

def find_crd_files(directory='.'):
    file_list = []
    for file in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, file)) and (file.endswith('.crd') or file.endswith('.rst')):
            file_list.append(os.path.join(directory, file))
    return file_list

"""
def find_nc_files(directory='.'):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.nc'):
                file_list.append(os.path.join(root, file))
    return file_list

def find_top_files(directory='.'):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.top'):
                file_list.append(os.path.join(root, file))
    return file_list

def find_crd_files(directory='.'):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.crd') or file.endswith('.rst'):
                file_list.append(os.path.join(root, file))
    return file_list
"""