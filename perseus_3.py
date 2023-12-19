import sys
from perseuspy import pd
from Bio import SeqIO
import re
import math
import os
import glob


def get_fasta_file():
    """searches for the FASTA-file in the current working directory"""
    dir_path = os.path.dirname(os.path.abspath(__file__))
    fasta_file = glob.glob(dir_path+r'\*.fasta')
    if len(fasta_file) > 1:
        raise FileExistsError('More than one fasta-file found')
    if len(fasta_file) == 0:
        raise FileNotFoundError('No fasta-file found')
    return fasta_file[0]


def read_FASTA_file(fasta_file):
    """converts FASTA-file to a dictionary termed FASTA_dict"""
    FASTA_dict = dict()
    seq_file = SeqIO.parse(fasta_file, 'fasta')
    for entry in seq_file:
        identifier = entry.id.split('|')[1]
        FASTA_dict[identifier] = str(entry.seq)
    return FASTA_dict


def read_df_and_get_sequnces(FASTA_dict, matrix_file):
    """reads a Perseus matrix file and adds the UniProt sequence based on the Protein.Group column"""
    def get_sequence(item):
        if ';' in item:
            item = item.split(';')[0]
        return FASTA_dict.get(item)
    df = pd.read_perseus(matrix_file)
    df['UniProt_sequence'] = df['Protein.Group'].apply(get_sequence)
    return df


def find_peptide_in_protein(df):
    """finds the position of the peptide in the protein sequence"""
    for index, row in df.iterrows():
        peptide = row['Stripped.Sequence']
        protein = row['UniProt_sequence']
        if protein:
            df.loc[index, 'position_in_protein'] = int(protein.find(peptide))
    return df


def get_UniMod_position_and_AA_one_item(item):
    match_positions = list()
    aa_of_interest = list()
    pattern = r'\(UniMod:(\d+)\)'
    matches = re.finditer(pattern, item)
    total_match_length = 0
    for match in matches:
        match_string = match.group()
        if match_string == '(UniMod:21)':
            match_position = match.span()[0]
            match_positions.append(match_position - total_match_length)
            aa_of_interest.append(item[match_position-1])
        total_match_length += len(match_string)
    return match_positions, aa_of_interest


def get_UniMod_position_and_AA_all_items(df):
    for index, row in df.iterrows():
        position_in_protein = row['position_in_protein']
        if not math.isnan(float(position_in_protein)):
            modified_sequence = row['Modified.Sequence']
            uni_mod_positions, aa_of_interest = get_UniMod_position_and_AA_one_item(modified_sequence)
            if uni_mod_positions:
                uni_mod_positions = [str(int(position_in_protein + i)) for i in uni_mod_positions]
                uni_mod_positions_text = ';'.join(uni_mod_positions)
                aa_of_interest = ';'.join(aa_of_interest)
                df.loc[index, 'mod21'] = uni_mod_positions_text
                df.loc[index, 'mod21_aa'] = aa_of_interest
    return df


def prepare_for_writing(df):
    df.drop(columns=['UniProt_sequence', 'position_in_protein'], inplace=True)
    df['mod21_aa'].replace('n', '', inplace=True)
    df['mod21'].replace('nan', '', inplace=True)
    return df


def perseus_version():
    """works as a perseus plugin"""
    _,  infile, outfile = sys.argv
    fasta_file = get_fasta_file()
    FASTA_dict = read_FASTA_file(fasta_file)
    matrix_file = infile
    df = read_df_and_get_sequnces(FASTA_dict, matrix_file)
    df = find_peptide_in_protein(df)
    df = get_UniMod_position_and_AA_all_items(df)
    df = prepare_for_writing(df)
    df.to_perseus(outfile)


def standalone_version():
    fasta_file = get_fasta_file()
    FASTA_dict = read_FASTA_file(fasta_file)
    matrix_file = 'report.pr_matrix.tsv'
    df = read_df_and_get_sequnces(FASTA_dict, matrix_file)
    df = find_peptide_in_protein(df)
    df = get_UniMod_position_and_AA_all_items(df)
    df = prepare_for_writing(df)
    df.to_csv('report.pr_matrix_mod_21.tsv', sep='\t')


if len(sys.argv) != 3:
    standalone_version()
else:
    perseus_version()
