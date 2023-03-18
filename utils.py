from PyExp import sc_iter_filepath_folder

REVCOMP_DICTIONARY = dict(zip('ATCGNatcgn~[]', 'TAGCNtagcn~]['))

def get_revcomp(sequence):
    '''Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    '''
    return ''.join(REVCOMP_DICTIONARY.get(nucleotide, '') for nucleotide in reversed(sequence))

    
def sc_iter_fasta_brute(file_name, inmem=False, lower=False):
    """ Iter over fasta file."""
    
    header = None
    seq = []
    if file_name.endswith(".gz"):
        opener = gzip.open
        decoder = lambda x: x.decode("utf8")
    else:
        opener = open
        decoder = lambda x: x
    
    with opener(file_name) as fh:
        if inmem:
            data = [decoder(x) for x in fh.readlines()]
            decoder = lambda x: x
        else:
            data = fh
        for line in data:
            line = decoder(line)
            if line.startswith(">"):
                if seq or header:
                    sequence = "".join(seq)
                    if lower:
                        sequence = sequence.lower()
                    yield header, sequence
                header = line.strip()
                seq = []
                continue
            seq.append(line.strip())
        if seq or header:
            sequence = "".join(seq)
            if lower:
                sequence = sequence.lower()
            yield header, sequence

