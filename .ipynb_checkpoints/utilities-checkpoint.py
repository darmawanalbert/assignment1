    
    
def import_16s_sequences():
    with open('COMP90014_16S_sequences.fasta', 'r') as fp:
        lines = fp.readlines()
        lines = [line.rstrip('\n') for line in lines]
    
    headers_sequences = {}
    header = ''
    sequence = ''
    for line in lines:
        if line.startswith('>'):
            if sequence != '':
                headers_sequences[header] = sequence
            header = line.lstrip('>')
            sequence = ''
        else:
            sequence += line

    return headers_sequences
                