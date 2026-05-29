import os
import json
import pandas as pd
import numpy as np
from numpy import random
import argparse


def create_parser():

    parser = argparse.ArgumentParser(
        description="Generate <-num-seeds> AlphaFold3 JSON input files for each case (row) in the input CSV file."
    )

    parser.add_argument(
        "--mode", "-m",
        type=str,
        required=True,
        choices=["make_json", "run_MSA", "run_inference"],
        help="Mode of operation. 'make_json' to generate AlphaFold3 JSON input files; 'run_MSA' to submit MSA generation jobs; 'run_inference' to submit AF3 inference jobs.",
    )

    parser.add_argument(
        "--tcrpmhc-specific", "-t",
        action='store_true',
        help="If set, uses template scripts that call TCRmodel2 Database for speed. Only applicable for 'run_MSA' and 'run_inference' modes."
    )
    
    parser.add_argument(
        "--input-csv", "-i",
        type=str,
        required=True,
        help="CSV file containing the cases and their sequences.",
    )

    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        required=True,
        help="Directory to save the generated AlphaFold3 JSON input files.",
    )

    parser.add_argument(
        "--num-seeds", "-n",
        type=int,
        default=5,
        help="Number of different runs (each using a different random seed) for each input case (i.e. each CSV row)."
    )

    parser.add_argument(
        "--ID-column", "-d",
        type=str,
        default="PDBID",
        help="Column name for case ID in the CSV file."
    )

    parser.add_argument(
        "--chainID-columns", "-c",
        type=str,
        nargs='+',
        default=["A", "B"],
        help="Column name(s) for chain IDs in the CSV file. Can provide multiple column names.")

    parser.add_argument(
        "--model-weights", "-w",
        type=str,
        help="Path to the directory containing the AF3 model weights. Only applicable for 'run_inference' mode."
    )

    return parser.parse_args()

def add_protein_chain(chain_id, seq):
    """
    Inputs :
        chain_id : str
            Single-character chain identifier used by AlphaFold3.
        seq : str
            One letter amino acid sequence of the protein chain.

    Output: Dictionary formatted for AlphaFold3 JSON input containing the
        protein chain ID and amino acid sequence.
    """
    return {
        "protein": {
            "id": chain_id,
            "sequence": seq,
        }
    }

def make_AF3_json(chains_dict, ID, output_dir, seednumber=None):
    """
    Inputs:
    chains_dict : dict
        Dictionary containing chain IDs as keys and one letter amino acid sequences as values.
    ID : str
        Unique identifier used as the AlphaFold3 job name and JSON filename.
    output_dir : str
        Directory path where the AlphaFold3 input JSON file is written.
    seednumber : int or None, optional
        Random seed for AlphaFold3 inference. If a numeric value is provided,
        it is used directly. If None, a random seed is generated (range = 1 - 10000).
        Defaults to None.

    Output:
        Writes a single AlphaFold3-compatible JSON file to disk containing
        the protein chains specified in the input dictionary.

    """
    #If no seednumber is provided, generate a random one
    if not isinstance(seednumber, (int, float)):
        rng = np.random.default_rng()
        seednumber = int(rng.integers(10000))

    #Template for Alphafold3 submission json file designed for AF3.
    af3_setup = {
        "name": ID,
        "modelSeeds": [seednumber],
        "sequences": [],
        "dialect": "alphafold3",  # Or "alphafold3" based on your need
        "version": 1
    }

    #Add chain IDs and sequences to the json structure
    for chain_ID, seq in chains_dict.items():
        af3_setup["sequences"].append(add_protein_chain(chain_ID, seq))

    #generate unique json file
    filename = f"{ID}.json"
    output_path = os.path.join(output_dir,filename)

    with open(output_path, 'w') as f:
        f.write(json.dumps(af3_setup, indent=2))

def submit_MSA(input_folder, case_ID, tcrpmhc_specific=False):
    '''
    Submits a job to generate MSAs for the given case using the template_data_process.sh script.
    Inputs:
        input_folder : str
            Directory containing the input JSON file(s) for the case.
        case_ID : str
            Unique identifier for the case, used in job naming and logging.
        tcrpmhc_specific : bool, optional
            If True, uses the template script calling TCR-pMHC optimized databases for speed. 
            Defaults to False.
    '''
    # Make log directory
    if not os.path.exists(f"{input_folder}/log"):
        os.makedirs(f"{input_folder}/log")

    if not os.path.exists(f"{input_folder}/output_AF3_MSA"):
        os.makedirs(f"{input_folder}/output_AF3_MSA")

    if tcrpmhc_specific:
        template_script = "./template_scripts/template_data_process_spedup.sh"
    else:
        template_script = "./template_scripts/template_data_process.sh"

    # Definse submission command
    command = (' ').join(["sbatch",  "--job-name", f"AF3_MSA_{case_ID}",
                            "--output", f"{input_folder}/log/%x_%j.out",
                            "--error", f"{input_folder}/log/%x_%j.err",
                            template_script,
                            f"{input_folder}",
                ])

    # Log submission command       
    with open(f"{input_folder}/log/{case_ID}_MSA_command.txt", 'w') as f:
        f.write(command)

    # Submit the job
    print(f"Submitting MSA job for case {case_ID} with command:\n{command}")
    os.popen(command).read()

def submit_AF3_inference(input_folder, case_ID, model_weights, tcrpmhc_specific=False):
    '''
    Submits a job to run AF3 inference for the given case using the template_inference_a100.sh script.
    Inputs:
        input_folder : str
            Directory containing the input JSON file(s) for the case.
        case_ID : str
            Unique identifier for the case, used in job naming and logging.
        tcrpmhc_specific : bool, optional
            If True, uses the template script calling TCR-pMHC optimized databases for speed.
            Defaults to False.
        model_weights : str, optional
            Path to the directory containing the AF3 model weights. Defaults to '/home/ddiepenbroek'.
    '''
    # Make log directory
    if not os.path.exists(f"{input_folder}/log"):
        os.makedirs(f"{input_folder}/log")

    # Definse submission command
    if tcrpmhc_specific:
        template_script = "./template_scripts/template_inference_a100_spedup.sh"
    else:
        template_script = "./template_scripts/template_inference_a100.sh"

    command = (' ').join(["sbatch",  "--job-name", f"AF3_inference_{case_ID}",
                            "--output", f"{input_folder}/log/%x_%j.out",
                            "--error", f"{input_folder}/log/%x_%j.err",
                            template_script,
                            f"{model_weights}",
                            f"{input_folder}",
                            
                ])

    # Log submission command       
    with open(f"{input_folder}/log/{case_ID}_AF3_inference_command.txt", 'w') as f:
        f.write(command)

    # Submit the job
    print(f"Submitting AF3 inference job for case {case_ID} with command:\n{command}")
    os.popen(command).read()


if __name__ == "__main__":

    args = create_parser()

    #Make output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    #Load data
    df = pd.read_csv(args.input_csv)

    #Collect all the cases 
    cases = {}
    for case in df.iloc():
        cases[case[args.ID_column]] = {x : case[x] for x in args.chainID_columns}

    #%%
    ## Run according to mode
    if args.mode == "make_json":

        #generate the AF3 submission files
        for case_ID, chains_dict in cases.items():
            if not os.path.exists(f"{args.output_dir}/{case_ID}"):
                os.makedirs(f"{args.output_dir}/{case_ID}")

            for i in range(args.num_seeds):
                make_AF3_json(chains_dict, f"{case_ID}_rs{i}", f"{args.output_dir}/{case_ID}",
                             seednumber=None)

        print(f"All done! JSON files generated in {args.output_dir}")

    elif args.mode == "run_MSA":
        # Submit MSA generation jobs
        for case_ID, chains_dict in cases.items():
            submit_MSA(f"{args.output_dir}/{case_ID}", case_ID, tcrpmhc_specific=args.tcrpmhc_specific)

    elif args.mode == "run_inference":
        # Submit AF3 inference jobs
        for case_ID, chains_dict in cases.items():
            submit_AF3_inference(f"{args.output_dir}/{case_ID}", case_ID, 
                                    model_weights=args.model_weights, tcrpmhc_specific=args.tcrpmhc_specific)