import re

# These are patches for inconsistencies in RSAT and Microbes Online
KEGG_EXCEPTIONS = { 'Pseudomonas aeruginosa PAO1': 'Pseudomonas aeruginosa',
                    'Campylobacter jejuni NCTC11168': 'Campylobacter jejuni' }

def patch_gene(code, gene):
    """Microbes Online genes names that differ from RSAT names are renamed here"""    
    if code == 'bth':
        if gene.startswith('p'):
            # p5482_01 -> BT_p548201
            return 'BT_' + gene.replace('_', '')
        elif re.compile('BT\d+').match(gene):
            # BT1234 -> BT_1234
            return gene.replace('BT', 'BT_')
    return gene
