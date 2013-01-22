# These are patches for inconsistencies in RSAT and Microbes Online
KEGG_EXCEPTIONS = { 'Pseudomonas aeruginosa PAO1': 'Pseudomonas aeruginosa',
                    'Campylobacter jejuni NCTC11168': 'Campylobacter jejuni' }

def patch_synonyms(synonyms_dict, code):
    """exceptions: patch synonyms for specific organisms"""
    if code == 'bth':
        # bth in Microbes Online is inconsistent with RSAT - gene names contain
        # an underscore while RSAT does not. We fix this by adding an additional
        # entry in the synoyms without the underscore
        new_dict = synonyms_dict.copy()
        for key, value in synonyms_dict.items():
            if key.startswith('BT_p'):
                # BT_p548201 -> p5482_01
                new_key = (key[:-2] + '_' + key[-2:]).replace('BT_', '')
                new_dict[new_key] = value
            elif key.startswith('BT_'):
                # BT_1234 -> BT1234
                new_dict[key.replace('BT_', 'BT')] = value
        return new_dict
    else:
        return synonyms_dict
