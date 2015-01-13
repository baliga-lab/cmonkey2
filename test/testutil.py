"""common functionality for tests"""
import util
import rsat
import stringdb
import microbes_online
import organism as org

KEGG_FILE_PATH = 'config/KEGG_taxonomy'
GO_FILE_PATH = 'config/proteome2taxid'
CACHE_DIR = 'cache'
RSAT_BASE_URL = 'http://embnet.ccg.unam.mx/rsa-tools'

def make_halo(search_distances, scan_distances, ratios=None):
    """returns the organism object to work on"""
    keggfile = util.read_dfile(KEGG_FILE_PATH, comment='#')
    gofile = util.read_dfile(GO_FILE_PATH)
    rsatdb = rsat.RsatDatabase(RSAT_BASE_URL, CACHE_DIR,
                               'Halobacterium sp', 64091)
    mo_db = microbes_online.MicrobesOnline(CACHE_DIR)
    stringfile = 'testdata/string_links_64091.tab'

    nw_factories = []
    if stringfile != None:
        nw_factories.append(stringdb.get_network_factory('hal', stringfile, 0.5))
    else:
        logging.warn("no STRING file specified !")

    if ratios is not None:
        nw_factories.append(microbes_online.get_network_factory(
            mo_db, max_operon_size=ratios.num_rows / 20, weight=0.5))

    keggorg = util.make_dfile_map(keggfile, 1, 3)['hal']
    rsat_info = org.RsatSpeciesInfo(rsatdb, keggorg, 'Halobacterium_sp', 64091)
    gotax = util.make_dfile_map(gofile, 0, 1)[rsat_info.go_species()]
    return org.Microbe('hal', keggorg, rsat_info, gotax, mo_db, nw_factories,
                       search_distances, scan_distances, True, None)
