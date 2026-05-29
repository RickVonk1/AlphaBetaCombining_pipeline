import glob
import os
import shutil
from argparse import ArgumentParser

"""
This copies the models in the alphafold3 output directory into a new output directory and renames them.

Reason: the model names of the local alphafold3 output is identical for every case. 

input:  model_dir = path to Alphafold3's output directory
        output_dir = the location where the models are copied to.

output: copied models with the directory name and sample number 
        example:    dir= /output_af3/8shi_rs2
                    name = 8shi_rs2_m* (* = sample number)

"""

parser = ArgumentParser(description="Copies and renames Alphafold3 models.")

parser.add_argument(
    "--model-dir", "-m",
    type=str,
    required=True,
    help="Path to the Alphafold3 output directory",
)
parser.add_argument(
    "--output-dir", "-o",
    type=str,
    required=True,
    help="Directory to save the copied and renamed models",
    )

args = parser.parse_args()
model_dir = args.model_dir
output_dir = args.output_dir

for dir1 in glob.glob(f"{model_dir}/*"):
    basename = os.path.basename(dir1)[-8:]
    called = os.path.basename(model_dir)
    for dir2 in glob.glob(f"{dir1}/*"):
        model_id = os.path.basename(dir2).split("-")[-1]
        af_id = os.path.basename(dir2)
        for path in glob.glob(f"{dir2}"):
            if path.endswith(".cif"):
                name = f"{called}_{basename}_{model_id}"
                destination_path = os.path.join(output_dir, name)
                shutil.copy(path,destination_path)

