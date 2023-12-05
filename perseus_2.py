import sys
from perseuspy import pd
from Bio import SeqIO
import re
import math
import os
import glob


dir_path = os.path.dirname(os.path.abspath(__file__))
fasta_file = glob.glob(dir_path+r'\*.fasta')[0]



_,  infile, outfile = sys.argv # read arguments from the command line
df = pd.read_perseus(infile)

print(os.getcwd())

sequence_dict = dict()
seq_file = SeqIO.parse(fasta_file, 'fasta')
for entry in seq_file:
    identifier = entry.id.split('|')[1]
    sequence = str(entry.seq)
    sequence_dict[identifier] = sequence


def get_sequence(item):
    if ';' in item:
        item = item.split(';')[0]
    return sequence_dict.get(item)


df['UniProt_sequence'] = df['Protein.Ids'].apply(get_sequence)

for index, row in df.iterrows():
    peptide = row['Stripped.Sequence']
    protein = row['UniProt_sequence']
    if protein:
        df.loc[index, 'position_in_protein'] = int(protein.find(peptide))


def get_UniMod_position(item):
    match_positions = list()
    aa_of_interest = list()
    pattern = r'\(UniMod:21\)'
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


for index, row in df.iterrows():
    position_in_protein = row['position_in_protein']
    if not math.isnan(float(position_in_protein)):
        modified_sequence = row['Modified.Sequence']
        uni_mod_positions, aa_of_interest = get_UniMod_position(modified_sequence)
        if uni_mod_positions:
            uni_mod_positions = [str(int(position_in_protein + i)) for i in uni_mod_positions]
            uni_mod_positions_text = ';'.join(uni_mod_positions)
            aa_of_interest = ';'.join(aa_of_interest)
            df.loc[index, 'mod21'] = uni_mod_positions_text
            df.loc[index, 'mod21_aa'] = aa_of_interest


df.drop(columns=['UniProt_sequence', 'position_in_protein'], inplace=True)
df['mod21_aa'].replace('n', '', inplace=True)
df['mod21'].replace('nan', '', inplace=True)
df.to_perseus(outfile)
print('Done. Good bye!')


if __name__ == '__main__':
    match_positions, aa_of_interest = get_UniMod_position(item='AAS(UniMod:21)PYSQRPAS(UniMod:21)PTAVR')
