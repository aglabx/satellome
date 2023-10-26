
REVCOMP_DICTIONARY = dict(zip('ATCGNatcgn~[]', 'TAGCNtagcn~]['))

def get_revcomp(sequence):
    '''Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    '''
    return ''.join(REVCOMP_DICTIONARY.get(nucleotide, '') for nucleotide in reversed(sequence))
