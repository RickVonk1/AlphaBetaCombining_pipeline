#------------- Explanation -------------#
''''
Author : Rick Vonk

This code is made for the Alpha/beta combining pipeline. 
The intent of this code is to standardize the method of tarring large data file locaitons created within the pipeline. 
The code required only 1 input: the directory wanting to tar.

Input:
    - The directory wanting to tar

Output:
    - A .tar.gz file with the name of the direcotry wanting to tar, in the directory above it.
    - The direcotry tarred, removed.

- Usage:
    - python tar_and_remove.py /path/to/dir

'''
# ----------------- Imports ------------------ #
import os
import subprocess
import sys
from pathlib import Path
# ----------------- code ------------------ #

def tar_remove(dir_to_tar):
    script_location = Path(__file__).resolve()
    template = os.path.join(script_location.parent.parent, '4_haddock', 'required_for_scripts', 'zipping_template.sh')

    dir_to_tar = os.path.abspath(dir_to_tar)
    parent_dir = os.path.dirname(dir_to_tar)
    folder_name = os.path.basename(dir_to_tar)
    
    files_dir = os.path.join(parent_dir, 'Intermediate_files/zipping')
    os.makedirs(files_dir, exist_ok=True)

    with open(template, 'r') as f:
        contents = f.read()
    
    contents = contents.replace('$Result_location', files_dir)
    contents = contents.replace('$main_loc', parent_dir)
    contents = contents.replace('$folder_to_zip', folder_name)
    
    if os.path.isfile(os.path.join(dir_to_tar, f'{folder_name}.tar.gz')):
        contents = contents.replace('$name', f'{folder_name}_2')
    else:
        contents = contents.replace('$name', f'{folder_name}')
        
    contents = contents.replace('$haddock_output_loc', dir_to_tar)

    script_path = os.path.join(files_dir, f'{folder_name}_AF_zipping')
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(contents)

    result = subprocess.run(['sbatch', script_path], capture_output=True, text=True)
    return result

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Only 1 input is required')
        sys.exit()

    tar_remove(sys.argv[1])