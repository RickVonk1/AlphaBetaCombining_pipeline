#------------- Explanation -------------#
"""
Author : Rick Vonk

This is the second script in the haddock3 pipeline and is able to make cfg files for each pdb entry and activate them in batches of 50.

Input:
    - The Experiment directory
    - An experiment name

Output:
    - Haddock-configs for all models
    - A haddock run job

"""
#------------ Import -------------#
from pathlib import Path
import os
import subprocess
import sys

#--------- input -----------#
working_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'

experiment = 'tcr_alpha_3'

#--------- Code -----------#


def dir_maker(pdb_dir):
    # Makes more directories and saves their location
    pdb_working_dir = Path(pdb_dir)
    exp_name = pdb_working_dir.parent.name

    config_dir = pdb_working_dir.parent/'Haddock_output'/'config_files'
    config_dir.mkdir(parents=True, exist_ok=True)

    runs_dir = pdb_working_dir.parent/'Haddock_output'/'runs'
    runs_dir.mkdir(parents=True, exist_ok=True)
    
    return config_dir, runs_dir, exp_name

def config_maker(template, pdb_dir, config_dir):
    # Makes config files for each structure

    pdb_working_dir = Path(pdb_dir)
    for path in pdb_working_dir.glob('*') :
        new_config = config_dir/f'Haddock3_workflow_{path.stem}.cfg'

        experiment_name = path.stem
        working_dir = pdb_working_dir.parent/'Haddock_output'/experiment_name
        working_dir.mkdir(parents=True, exist_ok=True)

        with open (template,'r', encoding='utf-8') as f:
            contents = f.read()
            
        contents = contents.replace('$dir', str(working_dir))
        contents = contents.replace('$Case', str(path))
        with open (new_config, 'w', encoding = 'utf-8') as f:
            f.write(contents)

    return

def runs_maker(experiment_name,runs_template, runs_dir, config_dir, working_dir, group_size=50 ):
    # This combined haddock run activations to have a set number activated

    files = sorted([str(config_dir / f) for f in os.listdir(config_dir)])
    iteration = 1
    for i in range(0,len(files), group_size):
        group = files[i:i + group_size]
        run_iteration = runs_dir/f'Haddock3_run_{iteration}.sh'
        command_block = "\n".join([f"haddock3 {cfg} &" for cfg in group])

        with open (runs_template,'r', encoding='utf-8') as f:
            contents = f.read()

        contents = contents.replace('$working_dir', str(working_dir))
        contents = contents.replace('$iteration', str(iteration))
        contents = contents.replace('$exp_name', str(experiment_name))
        contents = contents.replace('$cases_to_run', command_block)

        with open(run_iteration, 'w', encoding='utf-8') as f:
            f.write(contents)
        
        iteration += 1
    print(f'\nMade a runs file with iteration: {iteration-1} \n')
    return 

def Activation(runs_dir):
    # Activates the run files
    run = Path(runs_dir)
    for script in sorted(run.glob('*.sh')) :
        subprocess.run(['sbatch', str(script)])


#--------- Activation -----------#

script_location = Path(__file__).resolve()
script_dir = script_location.parent

template_dir = os.path.join(script_dir,"required_for_scripts")
runs_template_loc = os.path.join(template_dir,'runs_template.sh')
cfg_template_loc = os.path.join(template_dir, 'Template.cfg')
pdb_loc = os.path.join(working_dir,experiment,'input_pdb')

if __name__ == "__main__":

    config, runs, experiment_name = dir_maker(pdb_loc)
    config_maker(cfg_template_loc, pdb_loc, config)
    runs_maker(experiment_name, runs_template_loc, runs, config, working_dir)
    Activation(runs)




