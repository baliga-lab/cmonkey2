import re

"""
These are patches for inconsistencies in KEGG, RSAT, Microbes Online and STRING
This can be regarded as a database of inconsistencies which was obtained by
creating concrete cMonkey runs on various microbial organisms
"""

KEGG_EXCEPTIONS = {'Pseudomonas aeruginosa PAO1': 'Pseudomonas aeruginosa',
                   'Campylobacter jejuni NCTC11168': 'Campylobacter jejuni'}


def patch_mo_gene(code, gene):
    """Microbes Online genes names that differ from RSAT names are
    renamed here"""
    if code == 'bth':
        if gene.startswith('p'):
            # p5482_01 -> BT_p548201
            return 'BT_' + gene.replace('_', '')
        elif re.compile('BT\d+').match(gene):
            # BT1234 -> BT_1234
            return gene.replace('BT', 'BT_')
    elif code == 'son':
        if re.compile('SO\d+').match(gene):
            # SO1234 -> SO_1234
            return gene.replace('SO', 'SO_')

    return gene


def patch_string_gene(code, gene):
    """STRING names unfortunately can differ from RSAT as well"""
    if code == 'cac':
        if gene.startswith('CA_'):
            return gene.replace('CA_', 'CA')
    return gene


def patch_ncbi_taxonomy(taxonomy_id):
    """patch for NCBI overrides:
    """
    # E.coli should be mapped to the 511145 code
    # to keep things consistent
    if taxonomy_id == '83333':
        return "511145"
    else:
        return taxonomy_id
