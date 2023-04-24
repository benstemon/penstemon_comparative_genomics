import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    min_len = int(sys.argv[3])
except:
    min_len = 1
    
print(f'Input File: {input_file}')
print(f'Output File: {output_file}')
if min_len:
    print(f'Minimum Length Cutoff: {min_len}')
    
contigs_dict = {}
total_length = 0
saved_length = 0
removed_length = 0

total_contigs = 0
saved_contigs = 0
removed_contigs = 0

with open(input_file, 'r') as fasta_input:
    current_contig = None
    current_sequence = ''
    current_length = 0
    has_lines = True
    while has_lines:
        fasta_line = fasta_input.readline()
        if len(fasta_line) == 0:
            has_lines = False
        else:
            fasta_line = fasta_line.strip()
            if fasta_line.startswith('>'):
                print(total_contigs)
                if current_contig:
                    if not min_len or current_length >= min_len:
                        contig_data = (current_contig, current_sequence, current_length)
                        while current_length in contigs_dict:
                            current_length -= 1
                        contigs_dict[current_length] = contig_data
                        saved_length += current_length
                        saved_contigs += 1
                    else:
                        removed_length += current_length
                        removed_contigs += 1
                    total_length += current_length
                    total_contigs += 1
                current_contig = fasta_line
                current_sequence = ''
                current_length = 0
            else:
                current_sequence = current_sequence + fasta_line
                current_length += len(fasta_line)
    
    print('\n<End of File>\n')
    if not min_len or current_length >= min_len:
        contig_data = (current_contig, current_sequence, current_length)
        while current_length in contigs_dict:
            current_length -= 1
        contigs_dict[current_length] = contig_data
        saved_length += current_length
        saved_contigs += 1
    else:
        removed_length += current_length
        removed_contigs += 1
        total_length += current_length
        total_contigs += 1

contig_order = sorted(contigs_dict.keys(), reverse=True)

print("Retained scaffolds:")
for j in contig_order:
    print(contigs_dict[j][0])
    print(f"{len(contigs_dict[j][1])} bp")

print(" === Statistics ===")
print(f"Total: {total_length} bp in {total_contigs} contigs\n")
if min_len:
    print(f"Saved: {saved_length} bp in {saved_contigs} contigs")
    print(f"Removed: {removed_length} bp in {removed_contigs} contigs\n")
    saved_str = "saved "
print(f"Largest Contig: {max(contig_order)}bp")
print(f"Smallest {saved_str} Contig: {min(contig_order)} bp\n")


with open(output_file, 'w') as fasta_output:
    for i in contig_order:
        fasta_output.write(contigs_dict[i][0] + '\n')
        fasta_output.write(contigs_dict[i][1] + '\n\n')
