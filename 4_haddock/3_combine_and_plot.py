#------------- Explanation -------------#
"""
Author : Rick Vonk

This script is the last step in the haddock3 pipeline, it extracts relavant information, centralises it for the whole experiment, makes a plot and then then tars the Haddock_ouput file to reduce inode usage

Input:
    - The Experiment directory
    - An experiment name
    - the metrics you wan to use
    - if you want the files ot be zipped

Output:
    - if zipping than a tarred haddock_output file
    - a plot with output ranked on total
    - csv-file with plot data
    - txt-file with analysis

"""
#------------ Import -------------#
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os
import shutil
import subprocess
import tarfile
from pathlib import Path
#------------ locations / Inputs-------------#
exp_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'

Experiment_name = 'tcr_alpha_3'

# total possibilities : ['AF3_confidence', 'total_energy', 'elec_energy','deltaG', 'haddock_score']
metrics = ['AF3_confidence', 'haddock_score', 'deltaG']
zipping_post = True

#------------ Functions -------------#

def extract_haddock_results(haddock_output, prod=False):
    datasets = {'total_energy': {}, 'elec_energy': {}}
    prodigy = {'deltaG': {}} if prod else {} 
    for case in glob.glob(os.path.join(haddock_output, "*")):
        experiment_id = os.path.basename(case)
        if experiment_id in {'config_files', 'log','runs', 'store'}:
            continue
        try:
            pattern = os.path.join(case, "analysis", "*_caprieval_analysis", "capri_clt.tsv")
            files = glob.glob(pattern)
            
            df = pd.read_csv(files[0], sep="\t", comment="#")            
            datasets['total_energy'][experiment_id] = float(df["total"].iloc[0])
            datasets['elec_energy'][experiment_id] = float(df["elec"].iloc[0])
            if prod :
                prod_pattern = os.path.join(case, '*_prodigyprotein','prodigyprotein.tsv')
                file = glob.glob(prod_pattern)
                df = pd.read_csv(file[0], sep="\t", comment="#")
                prodigy['deltaG'][experiment_id] = float(df['score'].iloc[0])
            
        except Exception:
            if experiment_id not in global_errors:
                global_errors.append(experiment_id)
    return datasets, prodigy

def extract_AF3_confidence(haddock_output, csv_file):
    datasets = {'AF3_confidence': {}}
    confidence = {}
    df = pd.read_csv(csv_file).set_index('experiment_id')

    for path in glob.glob(os.path.join(haddock_output, "*")):
        experiment_id = os.path.basename(path)
        if experiment_id in {'config_files', 'log','runs', 'store'}:
            continue

        if experiment_id in df.index:
            if (df.at[experiment_id, 'docking_orientation'] == True) and (df.at[experiment_id, 'Reverse_true_docking'] == True):
                confidence[experiment_id] = df.at[experiment_id, 'AF3_confidence_score']
            else:
                if experiment_id not in global_errors:
                    global_errors.append(experiment_id)
        else:
            if experiment_id not in global_errors:
                global_errors.append(experiment_id)

    datasets['AF3_confidence'] = confidence
    return datasets

def Haddock_score(haddock_output):
    data = {'haddock_score':{}}
    for case in glob.glob(os.path.join(haddock_output, "*")):
        experiment_id = os.path.basename(case)
        if experiment_id in {'config_files', 'log','runs', 'store'}:
            continue
        pattern = glob.glob(os.path.join(case,'analysis', '3_caprieval_analysis', 'summary.tgz'))[0]
        with tarfile.open(pattern, 'r:gz') as tar:
            with tar.extractfile('model_1.pdb') as f:
                for line in f:
                    line_str = line.decode('utf-8')
                    if 'REMARK HADDOCK score:' in line_str:
                        H_score = float(line_str.split()[-1])
                        data['haddock_score'][experiment_id] = H_score
                        continue
    return data

def normalise_data(data_dicts, prod_data):
    norm_datasets = {}
    data = data_dicts | prod_data
    for name, dataset in data.items():
        filtered_items = {k: v for k, v in dataset.items() if k not in global_errors}
        if not filtered_items:
            continue

        values = list(filtered_items.values())
        maximum = max(values)
        minimum = min(values)
        range_val = (maximum - minimum) if maximum != minimum else 1
        
        norm_datasets[name] = {}
        for k, v in filtered_items.items():
            norm = (v - minimum) / range_val
            if name in ['elec_energy', 'total_energy','deltaG','haddock_score']:
                norm = 1 - norm
            norm_datasets[name][k] = round(norm, 3)
            
    return norm_datasets

def csv_writer(data_dict, normalised_data_dict,prod_data, output_dir, name):
    output_path = os.path.join(output_dir, f'{name}.csv')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'experiment', 'total_energy', 'elec_energy', 'confidence', 'deltaG', 'haddock_score',
            'norm_total_energy', 'norm_elec_energy', 'norm_AF3_confidence', 'norm_deltaG', 'norm_haddock_score'
        ])
        
        ids_x = sorted(set().union(*[d.keys() for d in normalised_data_dict.values()]))
        for exp in ids_x:
            writer.writerow([
                exp,
                data_dict.get('total_energy', {}).get(exp, ''),
                data_dict.get('elec_energy', {}).get(exp, ''),
                data_dict.get('AF3_confidence', {}).get(exp, ''),
                prod_data.get('deltaG', {}).get(exp, ''),
                data_dict.get('haddock_score', {}).get(exp, ''),
                normalised_data_dict.get('total_energy', {}).get(exp, ''),
                normalised_data_dict.get('elec_energy', {}).get(exp, ''),
                normalised_data_dict.get('AF3_confidence', {}).get(exp, ''),
                normalised_data_dict.get('deltaG', {}).get(exp, ''),
                normalised_data_dict.get('haddock_score', {}).get(exp, ''),
            ])

def get_plotting_data(normalised_data_dict, selection, post_sel=False, avg_sel=False):
    plot_data = {k: v for k, v in normalised_data_dict.items() if k in selection}
    all_ids = set().union(*[d.keys() for d in plot_data.values()])
    
    temp_groups = {}
    for i in all_ids:
        base_id = i.split('_rs')[0]
        if base_id not in temp_groups: temp_groups[base_id] = []
        temp_groups[base_id].append(i)

    final_processed_data = {}

    for base_id, sub_ids in temp_groups.items():
        if avg_sel:
            avg_components = {}
            for feat in selection:
                vals = [normalised_data_dict[feat].get(sid, 0) for sid in sub_ids]
                avg_components[feat] = sum(vals) / len(vals) if vals else 0
            final_processed_data[base_id] = avg_components

        elif post_sel:
            best_sub_id = max(sub_ids, key=lambda sid: sum(normalised_data_dict[f].get(sid, 0) for f in selection))
            final_processed_data[best_sub_id] = {f: normalised_data_dict[f].get(best_sub_id, 0) for f in selection}

    return final_processed_data

def generate_rank_plot(processed_data, output_dir, name, selection):
    sort_totals = {idx: sum(components.values()) for idx, components in processed_data.items()}
    ids = sorted(sort_totals, key=sort_totals.get, reverse=True)

    x = np.arange(len(ids))
    bottom = np.zeros(len(ids))
    # True combination from test set
    highlight = {'tcra2b3', 'tcra5b2', 'tcra1b1'}
    
    edge_colors = ['black' if any(h in i for h in highlight) else 'none' for i in ids]
    color_map = {
        'AF3_confidence': '#1f77b4',
        'total_energy': '#ff7f0e',
        'elec_energy': '#2ca02c',
        'deltaG' : '#e8e337',
        'haddock_score' : '#FC9483'
    }

    plt.figure(figsize=(12, 6))
    
    for dataset_name in selection:
        values = [processed_data[i].get(dataset_name, 0) for i in ids]        
        plt.bar(x, values, bottom=bottom, edgecolor=edge_colors, 
                linewidth=1.5, color=color_map.get(dataset_name, 'gray'), label=dataset_name)
        
        bottom += np.array(values)

    is_averaged = all('_rs' not in i for i in ids)
    title = f"Ranking for {name}, sorted on {'average' if is_averaged else 'best Model'} per case"
    
    plt.xticks(x, ids, rotation=45, ha='right', fontsize=6)
    plt.xlabel('Cases')
    plt.ylabel('Normalized Combined Score')
    plt.title(title)
    plt.legend()
    plt.tight_layout()

    suffix = "average" if is_averaged else "selected"
    plot_out = os.path.join(output_dir, f"{name}_{suffix}.png")
    os.makedirs(os.path.dirname(plot_out), exist_ok=True)
    plt.savefig(plot_out, dpi=300)

def run_analytics(haddock_output, error_list, output_dir, name, normalised_data, prod_data, prod=False):
    output_path = os.path.join(output_dir, f'{name}.txt')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        print('\n#--------------- Analytics ---------------#', file=f)
        # entry data
        pdb_loc = os.listdir(os.path.join(os.path.dirname(haddock_output), 'input_pdb'))
        print(f'Total number of pdb entries supplied is : {len(pdb_loc)}\n', file=f)

        if error_list:
            print(f'Total failed cases: {len(set(error_list))}\n', file=f)

        # Module counting
        counts = {}
        rates = {}
        base_modules = ['topoaa', 'mdref', 'prodigyligand', 'clustfcc', 'seletopclusts', 'caprieval']
        highlight = {'tcra2b3', 'tcra5b2', 'tcra1b1'}
        
        successful_ids = set(normalised_data.get('total_energy', {}).keys())

        for path in glob.glob(os.path.join(haddock_output, "*")):
            full_id = os.path.basename(path)
            if full_id in {'config_files', 'log','runs'}: 
                continue
            subdirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
            
            for m_name in base_modules:
                for s_dir in subdirs:
                    if s_dir.endswith(m_name):
                        module_key = f"{s_dir}"
                        counts[module_key] = counts.get(module_key, 0) + 1
            
            base_id = full_id.split('_rs')[0]
            if base_id not in rates:
                rates[base_id] = [0, 0]
            
            rates[base_id][0] += 1
            if full_id in successful_ids:
                rates[base_id][1] += 1
        
        print(f'Module presence counts: {counts}', file=f)
        print('Detailed result sorted by success rate, based on total amount of models:', file=f)
        sorted_rates = sorted(
            rates.items(), 
            key=lambda item: (item[1][1] / item[1][0] if item[1][0] > 0 else 0), 
            reverse=True
        )
        highest = max(total for exp, (total, success) in sorted_rates)
        for exp, (total, success) in sorted_rates:
            rate = (success / highest * 100) if total > 0 else 0
            if exp in highlight:
                print(f'{exp} success rate : {rate:.1f}% ({success}/{highest}) <-- "real"', file=f)
            else:
                print(f'{exp} success rate : {rate:.1f}% ({success}/{highest})', file=f)
    
        # Copy cfg file to results page
        try :
            shutil.copy2(f'/home/rvonk1/3_Jolanda_data_Pipeline/Script_pipeline/3_haddock/Configs/{name}.cfg',os.path.join(output_dir, f'{name}.cfg'))
        except Exception :
            print('\ncfg file could not be found, require manual copying', file=f)

        # prodigy results analysis
        if prod :
            scores = prod_data.get('deltaG', {})
            sorted_pre = sorted(scores.items(), key=lambda item: item[1])
            sorted_prod = [(k, v) for k, v in sorted_pre if k not in error_list]
            print(f'\nTop 10 results by deltaG:', file=f)
            for exp_id, score in sorted_prod[:15]:
                is_real = any(target in exp_id for target in highlight)
                if is_real:
                    print(f'{exp_id}: {score} kcal/mol <-- "real"', file=f)
                else:
                    print(f'{exp_id}: {score} kcal/mol', file=f)
            print('\nA DeltaG table: \n-8/-10 is moderate \n-11/-13 is strong \n-14+ is very strong',file=f)
        
    return rates

def zip_and_remove(haddock_output, zip_template, name):
    haddock_output = os.path.abspath(haddock_output)
    base_dir = os.path.dirname(haddock_output)
    folder_name = os.path.basename(haddock_output)
    
    Result_dir = os.path.join(os.path.dirname(haddock_output), 'Result')
    zip_loc = os.path.join(Result_dir, 'zipping')
    zipping_file = f'{zip_loc}/Zipping_file_{name}.sh'
    os.makedirs(zip_loc, exist_ok=True)
    
    with open(zip_template, 'r') as f:
        contents = f.read()
    
    contents = contents.replace('$Result_location', str(Result_dir))
    contents = contents.replace('$main_loc', str(base_dir))
    contents = contents.replace('$folder_to_zip', str(folder_name))
    
    if os.path.isfile(os.path.join(base_dir, 'Haddock_output.tar.gz')):
        contents = contents.replace('$name', 'Haddock_output_2')
    else:
        contents = contents.replace('$name', 'Haddock_output')
        
    contents = contents.replace('$haddock_output_loc', str(haddock_output))
    
    with open(zipping_file, 'w', encoding='utf-8') as f:
        f.write(contents)
    
    subprocess.run(['sbatch', zipping_file])

    print(f'Zipping and removing Haddock_output of {name}\nDo not touch the Haddock_output file until job is done')

#------------ Main Execution -------------#
confidence_csv = os.path.join(exp_dir,Experiment_name,'Intermediate_files', 'Confidence_docking_csv.csv')

script_location = Path(__file__).resolve()
script_dir = script_location.parent

template_dir = os.path.join(script_dir,"required_for_scripts")
zip_template_loc = os.path.join(template_dir,'zipping_template.sh')

Result_output = os.path.join(exp_dir, Experiment_name, 'Result')
working_dir = os.path.join(exp_dir, Experiment_name, 'Haddock_output')

global_errors = []

Post_selection = False
Average_slection = True
Prodigy = True

if __name__ == "__main__":
    
    # 1. Extraction
    haddock_data, prod_data = extract_haddock_results(working_dir, Prodigy)
    af3_data = extract_AF3_confidence(working_dir, confidence_csv)
    H_score = Haddock_score(working_dir)
    
    combined_raw_data = af3_data | haddock_data | H_score 
    
    # 2. Normalization & CSV
    norm_datasets = normalise_data(combined_raw_data, prod_data)
    csv_writer(combined_raw_data, norm_datasets, prod_data, Result_output, Experiment_name)   

    # 3. Plotting
    totals = get_plotting_data(norm_datasets, metrics, Post_selection, Average_slection)
    generate_rank_plot(totals, Result_output, Experiment_name, metrics)

    # 4. Analytics
    run_analytics(working_dir, global_errors, Result_output, Experiment_name, norm_datasets, prod_data, Prodigy)

    #5. Cleaning
    if zipping_post:
        zip_and_remove(working_dir, zip_template_loc, Experiment_name)