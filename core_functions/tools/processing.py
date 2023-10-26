#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 23.02.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

REVCOMP_DICTIONARY = dict(zip('ATCGNatcgn~[]', 'TAGCNtagcn~]['))

def get_revcomp(sequence):
    '''Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    '''
    return ''.join(REVCOMP_DICTIONARY.get(nucleotide, '') for nucleotide in reversed(sequence))
