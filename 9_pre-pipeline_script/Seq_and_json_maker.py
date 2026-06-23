#------------- Explanation -------------#
"""
Author : Rick Vonk

This script is for the extraction and preparation of data for the Jolanda data itself.
This script expects a csv file with atleast :
Number,CDR3 sequence,V chain,D chain,J chain

Input:
    - csv-file with single chain data in the '/1_pre-AF3/required_files' folder

Output:
    - a list of json files for each chain.
    - a csv file with the murine combined sequences
    - an updates input file at '/1_pre-AF3/required_files' with the new sequences
"""
#------------ Import -------------#
import csv
import subprocess
import os
import json
import glob
from pathlib import Path
#----------------- inputs ------------------#

# ----------------- code ------------------ #

def stitchr_json(output_dir, csv_file):
    os.makedirs(output_dir, exist_ok=True)
    error_list = {}

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            output_name = os.path.join(output_dir, row['Number'])
            command = [
                'stitchr',
                '-v', row['V chain'],
                '-j', row['J chain'],
                '-cdr3', row['CDR3 sequence'],
                '-n', output_name,
                '-m', 'json'
            ]

            try:
                print(f"Processing: {row['Number']}...")
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                error_list[row['Number']] = (row['V chain'], row['J chain'], row['CDR3 sequence'], str(e))

        for k,v in error_list.items():
            print(k,v)
    return error_list

def translate_codon(codon):
    codon_map = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    return codon_map.get(codon.upper(), "[Invalid Codon]")

def get_translated_list(sequence):
    num_full_codons = len(sequence) // 3
    return "".join([
        translate_codon(sequence[i:i+3]) 
        for i in range(0, num_full_codons * 3, 3)
    ])

def json_extraction(output_dir):
    case_info = {}
    for case in glob.glob(output_dir + '/*.json'):
        with open(case, 'r') as f:
            data = json.load(f)
            ID = os.path.basename(case).replace('.json', '')
            case_info[ID] = {}
            seq_data = data.get('seqs', {})
            case_info[ID]['l'] = get_translated_list(seq_data.get('l', "").strip().replace(" ", ""))
            case_info[ID]['v'] = get_translated_list(seq_data.get('v', "").strip().replace(" ", ""))
            case_info[ID]['j'] = get_translated_list(seq_data.get('j', "").strip().replace(" ", ""))
            case_info[ID]['c'] = get_translated_list(seq_data.get('c', "").strip().replace(" ", ""))

    return case_info

def Manual_addition(case_info):
    # The following cases failed with the following info:
    # a8 ('TRAV38-2', 'TRAJ31', 'CAYRSALDNARLMF')
    # b12 ('TRBV29-1', 'TRBJ2-7*01', 'CSVEVGMGLTYEQY')
    # Info below was obtained by putting the gene in google's Gemini LLM
    case_info['a8'] = {}
    case_info['a8']['l'] = 'MKKLLAMILWLQLDRLSGE'
    case_info['a8']['v'] = 'LKVEQNPRFLITVKEGKNATLICEVTVPSTTATLQWFRQNRGKGLEFLIYYNNGEKEDGRFTAQVDKSSKYISLFIRDSQPSDSATYLCAM'
    case_info['a8']['j'] = 'NNNARLMFGDGTQLVVKP'
    case_info['a8']['c'] = ''

    case_info['b12'] = {}
    case_info['b12']['l'] = 'MGPQLLGYVVLCLLGAGPLEA'
    case_info['b12']['v'] = 'QVTQNPRYLITVTGKKLTVTCSQNMNHEYMSWYRQDPGLGLRQIYYSMNVEVTDKGDVPEGYKVSRKEKRNFPLILESPSPNQTSLYFCASS'
    case_info['b12']['j'] = 'SYEQYFGPGTRLTVT'
    case_info['b12']['c'] = ''

    return case_info

def csv_maker(output_dir, case_info):
    csv_file_loc = os.path.join(os.path.dirname(os.path.abspath(output_dir)), 'murine_combined_seq_CSV_file.csv')
    with open (csv_file_loc, 'w', newline='') as f:
        csv_file = csv.writer(f)
        murine_c_seq_a = 'DIQNPEPAVYQLKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDKCVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNATYPSSDVPCDATLTEKSFETDMNLNFQNLSVMGLRILLLKVAGFNLLMTLRLWSS'
        murine_l_seq_a = ''
        murine_c_seq_b = 'EDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVCTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADCGITSASYHQGVLSATILYEILLGKATLYAVLVSGLVLMAMVKKKNSGSG'
        murine_l_seq_b = ''
        csv_file.writerow(['ID','sequence'])
        for ID, seq in case_info.items():
            if 'a' in ID:
                sequence = f"{murine_l_seq_a}{seq['v']}{seq['j']}{murine_c_seq_a}"
            if 'b' in ID:
                sequence = f"{murine_l_seq_b}{seq['v']}{seq['j']}{murine_c_seq_b}"
            csv_file.writerow([ID,sequence])
    
    return csv_file_loc

def csv_updater(csv_file, data_csv):
    c_dict = {}
    with open(data_csv, 'r', newline='') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row:
                c_dict[row[0]] = row[1]

    csv_content = []
    with open(csv_file, 'r', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            csv_content.append(row)
    
    for row in csv_content[1:]:
        if row[0] in c_dict:
            row[11] = c_dict[row[0]]
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(csv_content)
    
    return

#----------------- Activation ------------------#

script_location = Path(__file__).resolve()
json_dir = os.path.join(script_location.parent, 'json/')
ab_csv = os.path.join(script_location.parent.parent, '1_pre-AF3','required_files', 'Alpha_Beta_single_chains.csv')

if __name__ == '__main__':
    errors = stitchr_json(json_dir, ab_csv)
    data = json_extraction(json_dir)
    data = Manual_addition(data)
    csv_loc = csv_maker(json_dir, data)
    csv_updater(ab_csv, csv_loc)
