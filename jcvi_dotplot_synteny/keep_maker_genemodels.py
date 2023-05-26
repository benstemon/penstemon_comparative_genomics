import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        line = line.strip()
        if line.startswith('##'):
            # Retain GFF version header
            if line.startswith('##gff-version'):
                print(line, file=f_out)
        else:
            fields = line.split('\t')
            if len(fields) >= 2 and fields[1] == 'maker':
                print(line, file=f_out)
