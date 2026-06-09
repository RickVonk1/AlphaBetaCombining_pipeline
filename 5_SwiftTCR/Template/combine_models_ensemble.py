import glob
import os
import pandas as pd
from collections import defaultdict
import shutil
import sys
import argparse
import tarfile
from tqdm import tqdm

# Local pipeline scripts in the same folder
import clustering
import pairwise_rmsd
# From lthijs

"""
SwiftTCR FT-run combiner + optional clustering pipeline

This script merges multiple SwiftTCR docking runs per PDB case into a single ranked dataset.
For each run directory inside the input folder, it reads the ft.000.00 file, selects the
top 1000 models, and combines them across runs belonging to the same PDB ID (first 4
characters of the directory name). The combined list is re-ranked by "total weighted energy"
and optionally limited to a maximum of 5000 models for memory stability.

For each PDB case, the script then extracts the corresponding PDB structures from each
run's merged.tar archive into a temporary working directory. The extracted PDBs are renamed
according to their new global ranking index.

Optionally, the script performs pairwise RMSD calculation (using SwiftTCR's pairwise_rmsd)
and clustering (using SwiftTCR's clustering module) on the combined set of models, producing
irmsd.csv and clustering.txt.

Finally, it writes a new output folder per PDB containing:
  - a combined merged.tar archive with the re-ranked PDB structures
  - a combined ft.000.00 file with the updated ranking
  - optional analysis outputs (irmsd.csv and clustering.txt)

"""

def main(input_dir, output_dir, cores, cluster=True, pdb_filter=None, run_filter=None):
    merged_dfs = combine_ft_runs(input_dir, pdb_filter=pdb_filter, run_filter=run_filter)

    for pdb_id, merged_df in merged_dfs.items():
        final_dest = os.path.join(output_dir, pdb_id)
        if os.path.exists(final_dest):
            print(f"\nSkipping {pdb_id}: Output map bestaat al in {output_dir}")
            continue

        print(f"\n--- Verwerken van PDB: {pdb_id} ---")

        # Work map TMPDIR
        tmp_pdb_work = os.path.join(os.environ.get("TMPDIR", "/tmp"), f"work_{pdb_id}")
        if os.path.exists(tmp_pdb_work):
            shutil.rmtree(tmp_pdb_work)
        os.makedirs(tmp_pdb_work, exist_ok=True)

        # 1. PDB's extraheren
        copy_merged_to_temp(merged_df, input_dir, tmp_pdb_work)
        print(f"PDB's extracted")
        
        # NIEUW: Schrijf de ft.000.00 alvast in de tmp_dir zodat de clustering hem kan lezen
        temp_ft = os.path.join(tmp_pdb_work, "ft.000.00")
        merged_df.to_csv(temp_ft, sep="\t", header=False, index=False)
        
        # 2. Clustering (kan nu de ft.000.00 file vinden)
        if cluster:
            swifttcr_clustering_combined(pdb_id, tmp_pdb_work, 3, cores)
            
        # 3. Verplaats alles naar de definitieve output map
        finalize_output(pdb_id, tmp_pdb_work, output_dir, merged_df)
        
        # 4. Tidying up
        shutil.rmtree(tmp_pdb_work)
        print(f"--- Klaar met {pdb_id} ---")

def combine_ft_runs(input_dir, pdb_filter=None, run_filter=None):
    headers = ["Rotation Index", "Translation (x)", "Translation (y)", "Translation (z)",
               "total weighted energy", "repulsive vdW energy (unweighted)",
               "attractive vdW energy (unweighted)", "coulombic electrostatic energy (unweighted)",
               "generalized Born approximation electrostatics energy (unweighted)",
               "pairwise potential energy (unweighted)"]
    dfs = defaultdict(list)
    print("Running combine_ft_runs")

    for dir1 in glob.glob(f"{input_dir}/*"):
        print(f"{dir1} exists")
        if not os.path.isdir(dir1): continue
        
        run_name = os.path.basename(dir1).lower()
        pdb_id = run_name[:4]
        
        if run_filter and run_name not in run_filter: continue
        if pdb_filter and pdb_id not in [p.lower() for p in pdb_filter]: continue

        ft_path = os.path.join(dir1, "ft.000.00")
        if os.path.exists(ft_path) and os.path.getsize(ft_path) > 0:
            print(f"{ft_path} exists and size is {os.path.getsize(ft_path)}")
            df = pd.read_csv(ft_path, sep="\t", header=None)
            df.columns = headers
            df = df.iloc[:1000] # Takes top 1000 models per run
            df["source_dir"] = os.path.basename(dir1)
            df["original_index"] = df.index
            dfs[pdb_id].append(df)

    merged_ranked = {}
    
    for pdb_id, df_list in dfs.items():
        print(f"{dfs.keys()} are the keys of created dfs")
        merged_df = pd.concat(df_list, ignore_index=True)
        merged_df = merged_df.sort_values("total weighted energy").reset_index(drop=True)
        
        # MEMORY PROTECTION: Limit to a total of x for clustering to avoid OOM
        if len(merged_df) > 5000:
            print(f"Note: {pdb_id} has {len(merged_df)} models. Limiting to 5000 for stability.")
            merged_df = merged_df.iloc[:5000]
            

        merged_df["new_id"] = merged_df.index
        merged_ranked[pdb_id] = merged_df
    print(f"Merged ranked is created with keys: {merged_ranked.keys()}")    
    return merged_ranked

def get_models(merged_df, inDIR):
    # inDIR contains the path to the directory of untarred merged files in subdirectories
    # merged_df contains path of 5000 lowest energie models

    # Path to your tar file
    tar_path = "archive.tar"
    
    # List of files you want to extract (paths inside the tar archive)
    files_to_extract = ["file1.txt", "folder/file2.csv"]

    if os.path.exists(tar_path):
        if tar_path not in opened_tars:
            opened_tars[tar_path] = tarfile.open(tar_path)
            

    # Path to extract files to
    extract_path = "./extracted_files"

    # Open the tar file in read mode
    with tarfile.open(tar_path, "r") as tar:
        # Iterate over each member
        for member in tar.getmembers():
            if member.name in files_to_extract:
                tar.extract(member, path=extract_path)
                print(f"Extracted: {member.name}")


    for _, case in merged_df.iterrows():
        source_run_dir = os.path.join(input_dir, case.source_dir)
        tarFl = os.path.join(source_run_dir, "merged.tar")
        
        # De structuur binnen de tar is 'merged/merged_X.pdb'
        in_file_name = f"merged_{case.original_index}.pdb"
        out_file = os.path.join(out_merg, f"merged_{case.new_id}.pdb")

        with tarfile.open(tarFl, 'r') as tar:
        # Get the TarInfo object for 'docs/report.pdf'
            members = tar.getmembers()
            for m in members:
                file_name = os.path.basename(m.name)
                print(file_name)
                if file_name == in_file_name:
                # if file_name != "merged":
                    out_file = f"/scratch-local/75131/{file_name}"
                    f = tar.extractfile(m)
                    with open(out_file, 'wb') as target:
                        print(f"{out_file}")
                        target.write(f.read())


def copy_merged_to_temp(merged_df, input_dir, tmp_pdb_work):
    
    out_merg = os.path.join(tmp_pdb_work, "merged")
    os.makedirs(out_merg, exist_ok=True)
    
    opened_tars = {}
    print(f"Extracting {len(merged_df)} PDBs from tar files into {out_merg}")
    
    for _, case in tqdm(merged_df.iterrows()):
        source_run_dir = os.path.join(input_dir, case.source_dir)
        tar_path = os.path.join(source_run_dir, "merged.tar")
        
        # De structuur binnen de tar is 'merged/merged_X.pdb'
        in_file_name = f"merged/merged_{case.original_index}.pdb"
        out_file = os.path.join(out_merg, f"merged_{case.new_id}.pdb")
        print(f"inputfile: {in_file_name}\nOutputfile: {out_file}")


        if os.path.exists(tar_path):
            if tar_path not in opened_tars:
                opened_tars[tar_path] = tarfile.open(tar_path)
            
            try:
                member = opened_tars[tar_path].getmember(in_file_name)
                f = opened_tars[tar_path].extractfile(member)
                with open(out_file, 'wb') as target:
                    target.write(f.read())

            except KeyError:
                # Soms zit het bestand zonder 'merged/' prefix in de tar
                try:
                    alt_name = f"merged_{case.original_index}.pdb"
                    member = opened_tars[tar_path].getmember(alt_name)
                    f = opened_tars[tar_path].extractfile(member)
                    with open(out_file, 'wb') as target:
                        target.write(f.read())
                except KeyError:
                    print(f"Error: {in_file_name} niet gevonden in {tar_path}")

    for t in opened_tars.values():
        t.close()

def finalize_output(pdb_id, tmp_pdb_work, output_dir, merged_df):
    final_dest = os.path.join(output_dir, pdb_id)
    os.makedirs(final_dest, exist_ok=True)
    
    # 1. Create new combined merged.tar
    merged_dir = os.path.join(tmp_pdb_work, "merged")
    if os.path.exists(merged_dir):
        with tarfile.open(os.path.join(final_dest, "merged.tar"), "w") as tar:
            tar.add(merged_dir, arcname="merged")
        
    # 2. Saves new ft.000.00 op
    merged_df.to_csv(os.path.join(final_dest, "ft.000.00"), sep="\t", header=False, index=False)
    
    # 3. Copy analysis results
    for f_name in ["irmsd.csv", "clustering.txt"]:
        src = os.path.join(tmp_pdb_work, f_name)
        if os.path.exists(src):
            shutil.copy(src, final_dest)

# from combine_models_ensemble_old.py
def swifttcr_clustering_combined(pdb_id, work_dir, threshold, cores):
    print("Clustering started")
    os.chdir(work_dir)
    if not os.path.exists("merged") or len(os.listdir("merged")) == 0:
        print("Geen PDBs gevonden om te clusteren.")
        return

    print(f"Start pairwise rmsd met {cores} cores")
    pairwise_rmsd.calc_rmsd("merged", "irmsd.csv", "A", "D", 10, n_cores=cores)
    clustering.clustering_main("irmsd.csv", threshold, "clustering.txt")

# def swifttcr_clustering_combined(pdb_id, work_dir, threshold, cores):
#     os.chdir(work_dir)
#     if not os.path.exists("merged") or len(os.listdir("merged")) == 0:
#         print("Geen PDBs gevonden om te clusteren.")
#         return
        
#     print(f"Start pairwise rmsd met {cores} cores")

#     # 1. Bereken de afstanden (iRMSD) tussen alle gecombineerde modellen
#     pairwise_rmsd.calc_rmsd("merged", "irmsd.csv", "A", "D", 10, n_cores=cores)
    
#     # 2. De gecombineerde energie-data staat al klaar (door de main functie gemaakt)
#     # In finalize_output wordt merged_df opgeslagen als ft.000.00, 
#     # maar we hebben het hier al nodig voor de ranking.
#     energy_file = "ft.000.00" 
    
#     # Check of de energie-file al is aangemaakt, anders even tijdelijk aanmaken
#     # (Dit is nodig omdat in je huidige main 'finalize_output' pas NA clustering komt)
    
#     print(f"Start Ranked Clustering (Density Scoring) voor {pdb_id}...")
    
#     # Belangrijk: Gebruik de nieuwe functie naam en geef de energy_file mee
#     clustering.clustering_main_ranked(
#         input_file="irmsd.csv", 
#         energy_file=energy_file, 
#         threshold=threshold, 
#         output_file="clustering.txt"
#     )
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--cores", type=int, default=24)
    parser.add_argument("--pdb", nargs="+")
    args = parser.parse_args()
    print("Running combine_models_ensemble.py")

    main(input_dir=args.input_dir, output_dir=args.output_dir, cores=args.cores, pdb_filter=args.pdb)