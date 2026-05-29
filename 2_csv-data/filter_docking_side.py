"""Postfilter

Author: Jan Aarts, Dario Marzella
"""
from pdb2sql import pdb2sql
import pandas as pd
import numpy as np
import math
import os
from glob import glob
from joblib import Parallel, delayed


AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def aa3to1(res_name):
    """Convert three-letter amino acid code to one-letter code."""
    return AA3_TO_1.get(res_name.upper(), 'X')  # 'X' for unknown residues
            
def postfilter(path_to_pdb, pept_seq, cdr3a_seq, cdr3b_seq, 
               min_num_residues = 3,  threshold= 8, verbose=True,
               criteria='and', tcr_dist_threshold=25):
    """Boolean function to check if TCR is on top of pMHC

    Checks if TCR is docked on top of pMHC.
    Returns True if min_num_residues residue in cdr3 loops is within <threshold> distance of any residue of peptide.
    Ignores multiple matches for now and takes the first one. #TODO fix multiple matches
    """
    pdb = pdb2sql(path_to_pdb)
    chain_dict = lookup_sequence(pdb)
    #find coordinates peptide
    peptide = find_sequence_in_pdb(chain_dict, pept_seq)
    if peptide == {}:
        if verbose:
            print(f"Peptide sequence not found: {pept_seq}, for {os.path.basename(path_to_pdb)}")
        return False
    try:
        first_chain_peptide, first_numbering_peptide = list(peptide.items())[0]
        start_pep_i, end_pep_i = first_numbering_peptide[0]
    except Exception as e:
        print("Unexpected error ", e)
        return False
    if len(peptide.items()) > 1:
        raise Exception(f'More than one sequence instance for peptide found.')
    coordinates_peptide = pdb.get('x,y,z', chainID=[first_chain_peptide], resSeq=[i for i in range(start_pep_i, end_pep_i+1)], name= "CA")
    
    #find coordinates cdr3a
    cdr3a = find_sequence_in_pdb(chain_dict, cdr3a_seq)
    if cdr3a == {}:
        if verbose:
            print(f"Cdr3a sequence not found: {cdr3a_seq} for {os.path.basename(path_to_pdb)}")
        return False
    try:
        first_chain_cdr3a, first_numbering_cdr3a = list(cdr3a.items())[0]
        start_cdr3a_i, end_cdr3a_i = first_numbering_cdr3a[0]
    except Exception as e:
        print("Unexpected error ", e)
        return False
    if len(cdr3a.items()) > 1:
        raise Exception(f'More than one sequence instance for CDR3a found in TCR')
    coordinates_cdr3a = pdb.get('x,y,z', chainID=[first_chain_cdr3a], resSeq=[i for i in range(start_cdr3a_i, end_cdr3a_i+1)], name= "CA")

    #find coordinates cdr3b
    cdr3b = find_sequence_in_pdb(chain_dict, cdr3b_seq)
    if cdr3b == {}:
        if verbose:
            print(f"Cdr3b sequence not found: {cdr3b_seq} for {os.path.basename(path_to_pdb)}")
        return False
    try:
        first_chain_cdr3b, first_numbering_cdr3b = list(cdr3b.items())[0]
        start_cdr3b_i, end_cdr3b_i = first_numbering_cdr3b[0]
    except Exception as e:
        print("Unexpected error ", e)
        return False
    if len(cdr3b.items()) > 1:
        raise Exception(f'More than one sequence instance for CDR3b found in TCR')
    coordinates_cdr3b = pdb.get('x,y,z', chainID=[first_chain_cdr3b], resSeq=[i for i in range(start_cdr3b_i, end_cdr3b_i+1)], name = "CA")

    #calc distance between any cdr3a and peptide
    distances_a_pep = calc_min_distance(coordinates_peptide, coordinates_cdr3a)
    
    #calc distance between any cdr3b and peptide
    distances_b_pep = calc_min_distance(coordinates_peptide, coordinates_cdr3b)

    #if verbose add print, match cdr3a or cdr3b
    passed_a = has_n_res_in_thr(distances_a_pep, threshold, min_num_residues)
    passed_b = has_n_res_in_thr(distances_b_pep, threshold, min_num_residues)
    
    tcr_centroids_distance = get_centroid_distance(pdb, first_chain_cdr3a, first_chain_cdr3b)
    passed_tcr_centroid = tcr_centroids_distance <= tcr_dist_threshold

    if criteria=='and':
        return passed_a and passed_b 
    elif criteria=='or':
        return passed_a or passed_b 
    else:
        raise Exception(f'Criteria not recognized {criteria}')

def has_n_res_in_thr(distances, threshold, min_num_residues):
    dist_sort= sorted(distances)
    subset = dist_sort[:min_num_residues]
    if subset[-1]> threshold:
        return False
    else:
        return True

def postfilter_anarci(path_to_pdb, peptide_chain, peptide_resids, 
                      cdr3a_chain, cdr3a_resids, cdr3b_chain, cdr3b_resids,
                      min_num_residues=3, dist_threshold=6, verbose=True,
                      criteria='and'):
    
    pdb = pdb2sql(path_to_pdb)
    coordinates_peptide = pdb.get('x,y,z', chainID=[peptide_chain], resSeq=peptide_resids, name="CA")
    coordinates_cdr3a = pdb.get('x,y,z', chainID=[cdr3a_chain], resSeq=cdr3a_resids, name="CA")
    coordinates_cdr3b = pdb.get('x,y,z', chainID=[cdr3b_chain], resSeq=cdr3b_resids, name="CA")
    #calc distance between any cdr3a and peptide
    distances_a_pep = calc_min_distance(coordinates_peptide, coordinates_cdr3a)
    
    #calc distance between any cdr3b and peptide
    distances_b_pep = calc_min_distance(coordinates_peptide, coordinates_cdr3b)
    #if verbose add print, match cdr3a or cdr3b
    passed_a = has_n_res_in_thr(distances_a_pep, dist_threshold, min_num_residues)
    passed_b = has_n_res_in_thr(distances_b_pep, dist_threshold, min_num_residues)
    if verbose:
        if not passed_a:
            print("Distance cdr3a and peptide is larger than: ", dist_threshold)
        if not passed_b:
            print("Distance cdr3b and peptide is larger than: ", dist_threshold)
    if criteria=='and':
        return passed_a and passed_b 
    elif criteria=='or':
        return passed_a or passed_b 
    else:
        raise Exception(f'Criteria not recognized {criteria}')

def get_centroid_distance(pdb, first_chain_cdr3a='D', first_chain_cdr3b='E'):
    tcra_coords = pdb.get('x,y,z', chainID=[first_chain_cdr3a], name="CA" )
    tcra_centroid = np.mean(np.asarray(tcra_coords), axis=0)
    tcrb_coords = pdb.get('x,y,z', chainID=[first_chain_cdr3b], name="CA" )
    tcrb_centroid = np.mean(np.asarray(tcrb_coords), axis=0)
    tcr_centroids_distance = float(np.linalg.norm(tcra_centroid - tcrb_centroid))
    return tcr_centroids_distance


def calc_min_distance(coordinates1, coordinates2):
    """Calculates the minimum distance between two coordinate clouds.

    Args: 
    coordinates1 (list), [[x11, x12, x13],[x21, x22, x23]]
    coordinates2 (list), [[x11, x12, x13],[x21, x22, x23]]
    Returns: min_dist, (float) minimum distance between two clouds.
    """
    coords1 = np.array(coordinates1)
    coords2 = np.array(coordinates2)
    distances = np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=2)
    min_distances = np.min(distances, axis=1)
    return min_distances

def calc_min_dist_pure_python(coordinates1, coordinates2):
    """Calculate distances without the dependence on numpy.

    """
    min_dist = float('inf')  # Initialize with a large value

    for x1, y1, z1 in coordinates1:
        for x2, y2, z2 in coordinates2:
            # Calculate Euclidean distance
            dist = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
            min_dist = min(min_dist, dist)  # Update minimum distance

    return min_dist


def find_sequence_in_pdb(chain_dict, sequence):
    """Finds residue numbering(s) matching a sequence in the pdb dictionary.

    Input: chain_dict (dict): (key: chainID, value: tup(sequence, numbering))
    sequence (str): amino acid sequence
    Returns: matches (dict): (key:chainID, value: list((start_index, end_index)))
    """
    matches = {}
    for chain, (pdb_sequence, numbering) in chain_dict.items():
        found_positions = []
        start_idx = pdb_sequence.find(sequence)
        while start_idx != -1:
            start_resnum = numbering[start_idx]
            end_resnum = numbering[start_idx +len(sequence)-1]
            found_positions.append((start_resnum, end_resnum))

            start_idx = pdb_sequence.find(sequence, start_idx +1)
        if found_positions:
            matches[chain] = found_positions
    return matches


def get_all_chains(pdb_file):
    """
    Get all chain identifiers in a PDB file using pdb2sql.

    Parameters:
    - pdb_file (str): Path to the PDB file.

    Returns:
    - List of chain identifiers (e.g., ['A', 'B', 'C']).
    """
    pdb = pdb2sql(pdb_file)
    chains = pdb.get_chains()
    return chains

def lookup_sequence(pdb):
    #pdb = pdb2sql(pdb_file)
    results = pdb.get_residues()
    chain_dict = {}

    for chain, residue_name, number in results:
        # Convert 3-letter AA code to 1-letter
        aa = aa3to1(residue_name)
        if chain not in chain_dict:
            chain_dict[chain] = ([],[])
        chain_dict[chain][0].append(aa)
        chain_dict[chain][1].append(number)
    for chain in chain_dict:
        chain_dict[chain] = ("".join(chain_dict[chain][0]), chain_dict[chain][1])
    return chain_dict


def test_cases(path_to_dir):
    #for this case CDR3a = CAVSGADKLIF, CDR3b = CASSNQAALRRLNTEAFF
    #peptide RLSSCVPVA
    for f in path_to_dir.iterdir():
        if f.suffix == ".pdb":
            passed = postfilter(f, "RLSSCVPVA", "CAVSGADKLIF", "CASSNQAALRRLNTEAFF", verbose=False)

            #passed1 = postfilter_anarci(f, "C", [1,2,3,4,5,6,7,8,9], "D", [i for i in range(106, 140)], "E", [j for j in range(106, 140)])

def get_first_docked_models(job_folder, min_num_residues, threshold, criteria):
    top_models = []
    for folder in glob(f'{job_folder}/*'):
        #print(folder)
        model = folder #ensures to check the models in ranking order
        CaseID = os.path.splitext(os.path.basename(folder))[0]
        #print(f'ranked_{i}')
        case_data = df[df['experiment_id']==CaseID]
        if postfilter(path_to_pdb=model, pept_seq=case_data['Peptide'].to_string(index=False), 
                    cdr3a_seq=case_data['CDR3a'].to_string(index=False), 
                    cdr3b_seq=case_data['CDR3b'].to_string(index=False), 
                    min_num_residues=min_num_residues,  threshold=threshold, verbose=False,
                    criteria=criteria):
            #return model
            #print(f'Rank {i} correct for case {CaseID}')
            top_models.append([1, model])
            break
        else:
            #print(folder)
            #print(f'ranked_{i}')
            #print('\n')
            pass
        
        if i == len(glob(f'{folder}/ranked_*.pdb'))-1:
            top_models.append(['FAILED', model])
    
    return top_models


if __name__ =="__main__":#remove when importing as library
    n_min_res = 2
    dist_threshold = 10
    criteria = 'and'

    data_csv = '/home/rvonk1/4_Haddock_config_experimentation/8_haddock3_complex_selection/other/Data_file.csv'
    df = pd.read_csv(data_csv)
    
    job  = []
    pdb_folder = '/projects/0/prjs1135/report_Rick/4_Haddock_config_experimentation/experiments/New_AF3_struc/AF3_process/AF3_step_2/TCR_post'
    for case in glob(pdb_folder) :
        jobs.append(case)


    #runs_folder = '/projects/0/prjs1135/immrep25/runs'
    #jobs = []
    #for pept_fol in glob(f'{runs_folder}/*'):
    #    for job in glob(f'{pept_fol}/job_?'):
    #        jobs.append(job)
    
    num_cores = 32
    top_models = Parallel(n_jobs = num_cores, verbose = 1)(delayed(get_first_docked_models)(job, 
                            min_num_residues=n_min_res,threshold=dist_threshold, criteria=criteria) for job in jobs)
    top_models = sum(top_models, [])
    
    out_df = pd.DataFrame(top_models, columns=['rank', 'path'])
    out_df.to_csv(f'/projects/0/prjs1135/report_Rick/4_Haddock_config_experimentation/8_Haddock3_complex_selection/csv_files/{n_min_res}t{dist_threshold}{criteria.upper()}.csv', index=False)
    out_df[out_df['rank']!=0].to_csv(f'/projects/0/prjs1135/report_Rick/4_Haddock_config_experimentation/8_Haddock3_complex_selection/csv_files/non_ranked_{n_min_res}t{dist_threshold}{criteria.upper()}.csv', index=False)
    
    
    '''
    n2t8OR = pd.read_csv('/projects/0/prjs1135/immrep25/postfiltering_data/distances/non_rank0_docked_models_n2t8OR.csv')
    n1t8AND = pd.read_csv('/projects/0/prjs1135/immrep25/postfiltering_data/distances/non_rank0_docked_models_n1t8AND.csv')
    merged_df = pd.merge(n2t8OR, n1t8AND, on='path')
    merged_df= merged_df.rename(columns={'rank_x':'n2t8OR', 'rank_y':'n1t8AND'})
    merged_df[merged_df['n2t8OR']!=merged_df['n1t8AND']].to_csv('/projects/0/prjs1135/immrep25/postfiltering_data/distances/n2t8OR_vs_n1t8AND.csv', index=False)
    '''