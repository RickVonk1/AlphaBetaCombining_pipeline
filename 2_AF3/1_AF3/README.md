
# Run modified database for AlphaFold3
> **note:** the AF3 model weights need to be requested at [AlphaFold3 repository](https://github.com/google-deepmind/alphafold3)

This folder contains all requirements to run AlphaFold3 (AF3) with a restricted database. constructed to reduce AF3 most time consuming step, the dataprocessing step.
The new database consists of reduced fasta files required for the Multiple Sequence Alignment (MSA) and a reduced mmcif directory. All based on [TCRmodel2 repository](https://github.com/piercelab/tcrmodel2). 

**Time difference**
- original database = ~2 hours,
- TCRmodel2 database = ~6 min

## Run Alphafold3 with the reduced database
AlphaFold3 must be runned using a submission script. Generally the data processing and inference step require different resources, as the dataprocessing step cannot speed up by allocating a GPU source. Each step will be explained in general, as each system requires alternative paths.

### Steps
1. generate an AlphaFold3 input `.json` file
    - example: `./test/input_8shi.json`
2. Modify `1_run_AF3_dataproc_smalldb.slurm` paths
    - log and output inside folder `./test`
3. Modify `2_run_AF3_inf_smalldb.slurm` output paths
    - log and output inside folder `./test`

This section shows the paths that need to be modified to run AlphaFold3 with a modified database.

**Modify paths  `1_run_AF3_dataproc_smalldb.slurm`:**

This is a simplefied version of `1_run_AF3_dataproc_smalldb.slurm`. Modify the paths to the right directory. The data_root and python_file are present in this folder. While the model root and json path is case dependent. 
```bash
#path to reduced database
data_root=/modified_database

# current modified script
python_file=run_alphafold.py

# path to directory of saved AF3 model weights
model_root=/home/ddiepenbroek

# path to submission script AF3
JSON_PATH=/test/input_8shi.json

cmd_args="--json_path ${JSON_PATH}
 --output_dir ${project_root}/test/
 --db_dir ${data_root}
 --pdb_database_path ${data_root}/mmcif_files_AF3_selection
 --model_dir ${model_root}"

# run commands with container
 cmd="singularity exec --bind ${bind_root} ${container_root} python ${python_file} ${cmd_args}" 
```

**Modify paths  `2_run_AF3_inf_smalldb.slurm`:**

Most modified paths can be copied from the `1_run_AF3_dataproc_smalldb.slurm`, except `JSON_path`. 

```bash 
# new path is a new file inside the output directory of AF3
folder_name=8shi_tcrmodel2_data
JSON_PATH="${project_root}/test/${folder_name}/${foder_name}_data.json"

```
The name provided in the input `json` will be the new folder and file name provided

## Overview setup database
As mentioned before the modified database is constucted based on [TCRmodel2's repository](https://github.com/piercelab/tcrmodel2). The overview will be split into a sequence database and mmcif database.

### Sequence Database
The database files were directly copied and the new file names are appended in the `run_alphafold.py` script. Each modification made in the script is listed below. 

**Modified Code:**
```py
_SMALL_BFD_DATABASE_PATH = flags.DEFINE_string(
    'small_bfd_database_path',
    '${DB_DIR}/small_bfd.tcrmhc.fasta',
    'Small BFD database path (tcrmodel2), used for protein MSA search of TCR-pMHC.',
)
_MGNIFY_DATABASE_PATH = flags.DEFINE_string(
    'mgnify_database_path',
    '${DB_DIR}/mgnify.fasta',
    'Mgnify database path (tcrmodel2), used for protein MSA search.',
)
_UNIPROT_CLUSTER_ANNOT_DATABASE_PATH = flags.DEFINE_string(
    'uniprot_cluster_annot_database_path',
    '${DB_DIR}/uniprot.tcrmhc.fasta',
    'UniProt database path (tcrmodel2), used for protein paired MSA search of TCR-pMHC.',
)
_UNIREF90_DATABASE_PATH = flags.DEFINE_string(
    'uniref90_database_path',
    '${DB_DIR}/uniref90.tcrmhc.fasta',
    'UniRef90 database path (tcrmodel2), used for MSA search. The MSA obtained by '
    'searching it is used to construct the profile for template search.',
)
_SEQRES_DATABASE_PATH = flags.DEFINE_string(
    'seqres_database_path',
    '${DB_DIR}/pdb_seqres.fasta',
    'PDB sequence database path (tcrmodel2), used for template search.',
)
```

**Original Code:**
```py
_SMALL_BFD_DATABASE_PATH = flags.DEFINE_string(
    'small_bfd_database_path',
    '${DB_DIR}/bfd-first_non_consensus_sequences.fasta',
    'Small BFD database path, used for protein MSA search.',
)
_MGNIFY_DATABASE_PATH = flags.DEFINE_string(
    'mgnify_database_path',
    '${DB_DIR}/mgy_clusters_2022_05.fa',
    'Mgnify database path, used for protein MSA search.',
)
_UNIPROT_CLUSTER_ANNOT_DATABASE_PATH = flags.DEFINE_string(
    'uniprot_cluster_annot_database_path',
    '${DB_DIR}/uniprot_all_2021_04.fa',
    'UniProt database path, used for protein paired MSA search.',
)
_UNIREF90_DATABASE_PATH = flags.DEFINE_string(
    'uniref90_database_path',
    '${DB_DIR}/uniref90_2022_05.fa',
    'UniRef90 database path, used for MSA search. The MSA obtained by '
    'searching it is used to construct the profile for template search.',
)
_SEQRES_DATABASE_PATH = flags.DEFINE_string(
    'seqres_database_path',
    '${DB_DIR}/pdb_seqres_2022_09_28.fasta',
    'PDB sequence database path, used for template search.',
)
```

