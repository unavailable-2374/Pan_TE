#!/usr/bin/env python3
import sys

# Define the classification hierarchy based on the image information.
classification_hierarchy = {
    'gypsy': 'LTR/Gypsy',
    'Copia': 'LTR/Copia',
    'Bel-Pao': 'LTR/Bel-Pao',
    'Retrovirus': 'LTR/Retrovirus',
    'DIRS': 'non-LTR/DIRS',
    'Ngaro': 'non-LTR/Ngaro',
    'VIPER': 'non-LTR/VIPER',
    'LINE': 'non-LTR/LINE',
    'SINE': 'non-LTR/SINE',
    'LTR': 'LTR',
    'non-LTR': 'non-LTR',
    'TIR': 'DNA/TIR',
    'Helitron': 'DNA/Helitron',
    'Maverick': 'DNA/Maverick',
    'R2': 'non-LTR/R2',
    'RTE': 'non-LTR/RTE',
    'Jockey': 'non-LTR/Jockey',
    'L1': 'non-LTR/L1',
    'Tc1-Mariner': 'DNA/Tc1-Mariner',
    'hAT': 'DNA/hAT',
    'Mutator': 'DNA/Mutator',
    'Merlin': 'DNA/Merlin',
    'PiggyBac': 'DNA/PiggyBac',
    'CACTA': 'DNA/CACTA',
    'SubclassI': 'DNA',
    'SubclassII': 'DNA',
}

def process_te_sequences(file_path, output_file_path):
    processed_lines = []
    with open(file_path, 'r') as file:
        # Skip the header line
        next(file)
        # Iterate over the file two lines at a time
        for line in file:
            id_line = line.strip()
            class_line = next(file, '').strip()

            if id_line.startswith('TE_'):
                # Get the classification keyword (remove comma)
                classification_keyword = class_line[1:]
                # Map the keyword to the full classification hierarchy
                full_classification = classification_hierarchy.get(classification_keyword, classification_keyword)
                
                # Concatenate the ID with the full classification
                processed_line = f"{id_line}#{full_classification}"
                processed_lines.append(processed_line)

    # Write the processed lines to the output file
    with open(output_file_path, 'w') as file:
        for line in processed_lines:
            file.write(line + '\n')

# Command line arguments for input and output file paths
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# Call the processing function
process_te_sequences(input_file_path, output_file_path)

