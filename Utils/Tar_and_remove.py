#------------- Explanation -------------#
''''
Author : Rick Vonk


This code is made for the Alpha/beta combining pipeline. The intent of this code s to standardize the method of tarring large data file locaitons created within the pipeline. The code required only 2 inputs: a template and the directory wanting to tar.

Input:
    - Template, standardized
    - The direcotry wanting to tar

Output:
    - A .tar.gz file with the anme fo the direcotry wantign to tar
    - The direcotry tarred, removed.

'''
# ----------------- code ------------------ #

def tar_remove(template, dir_to_tar):
    """File to zip SwiftTCR files to reduce sproject space"""
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

    subprocess.run(['sbatch', script_path])
    return
