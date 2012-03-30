'''
Created on Mar 21, 2012

@author: frank
'''


from gitForks.cmonkey_python.cmonkey.stringdb import get_network_factory_FS


file = "/Users/Frank/Data/Resources/external/String/protein.links.detailed.v9.0.txt"
speciesL = ['10090']

result = get_network_factory_FS(file, speciesL)