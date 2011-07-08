"""The data_providers module contains functionality related to
integrating with 3rd party data providers"""


def get_organism_for_code(kegg_taxonomy_file, code):
    """using a KEGG taxonomy file, lookup the organism code to
    return the full name. taxonomy_file should be an instance of
    DelimitedFile"""
    for line in kegg_taxonomy_file.get_lines():
        if line[1] == code:
            return line[3]
    return None

__all__ = ['get_organism_for_code']
