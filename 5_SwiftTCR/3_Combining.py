#------------- Explanation -------------#
"""
Author : Rick Vonk

Extract chain information for SwiftTCR from existing structure data

Input:
    -

Output:
    -

"""
#------------ Import -------------#
import os
import glob
import subprocess
import pandas as pd
from pathlib import Path


#----------------- inputs ------------------#
exp_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'

Experiment_name = 'tcr_alpha_11'

# ----------------- code ------------------ #
def combine_results(experiment_dir, experiment_name, template):
    swift_dir = os.path.join(experiment_dir, experiment_name,'SwiftTCR')
    for case in os.listdir(swift_dir):
        if 'runs_dir' in case:
            continue 
        bash_template = os.path.join(template,'combine_models_run.sh')
        python_script = os.path.join(template,'combine_models_ensemble.py')
        output_dir = os.path.join(swift_dir, case,'Swift_output')
        combined_dir = os.path.join(swift_dir, case,'Swift_combined')
        os.makedirs(combined_dir, exist_ok=True)


        with open(bash_template, 'r') as f:
            contents = f.read()
        
        contents = contents.replace('[output_dir]', output_dir)
        contents = contents.replace('[case]', case)
        contents = contents.replace('[combined_dir]', combined_dir)
        contents = contents.replace('[python_script]', python_script)

        combine_run = os.path.join (combined_dir, f'run_{case}.sh')
        with open (combine_run, 'w') as f:
            f.write(contents)

    return

def combine_activation(experiment_dir, experiment_name):
    combined_run = os.path.join(experiment_dir, experiment_name,'SwiftTCR')
    for case in os.listdir(combined_run):
        if 'runs_dir' in case :
            continue

        case_path = os.path.join(combined_run, case,'Swift_combined', f'run_{case}.sh')
        command = ['sbatch', case_path]

        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit {case}. Error: {e.stderr.strip()}")
        
    return


#----------------- Activation ------------------#
script_location = Path(__file__).resolve()
script_dir = script_location.parent
template = os.path.join(script_dir,'Template')


if __name__ == "__main__":
    combine_results(exp_dir, Experiment_name, template)
    combine_activation(exp_dir, Experiment_name)
