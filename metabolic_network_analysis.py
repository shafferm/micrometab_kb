from __future__ import division

import networkx as nx
from collections import defaultdict
import gzip
from os import path
from picrust.util import get_picrust_project_dir, convert_precalc_to_biom
import warnings
import requests
from multiprocessing.pool import ThreadPool


def load_data_table(ids_to_load):
    """Stolen from https://github.com/picrust/picrust/blob/master/scripts/predict_metagenomes.py
    and modified.
    Load a data table, detecting gziped files and subset loading
    data_table_fp -- path to the input data table
    load_data_table_in_biom -- if True, load the data table as a BIOM table rather
    than as tab-delimited
    suppress_subset_loading -- if True, load the entire table, rather than just
    ids_of_interest
    ids_to_load -- a list of OTU ids for which data should be loaded
    gzipped files are detected based on the '.gz' suffix.
    """
    # first two lines adapted from determine_data_table_fp from
    # https://github.com/picrust/picrust/blob/master/scripts/predict_metagenomes.py
    # stolen setup from predict_metagenomes.py from PICRUSt
    precalc_data_dir = path.join(get_picrust_project_dir(), 'picrust', 'data')
    precalc_file_name = '_'.join(["ko", "13_5", 'precalculated.tab.gz'])
    data_table_fp = path.join(precalc_data_dir, precalc_file_name)

    if not path.exists(data_table_fp):
        raise IOError("File " + data_table_fp + " doesn't exist! Did you forget to download it?")

    genome_table_fh = gzip.open(data_table_fp, 'rb')
    genome_table = convert_precalc_to_biom(genome_table_fh, ids_to_load)
    return genome_table


genes_seen = {}
def get_kegg_rxns_from_gene(gene):
    global genes_seen
    if gene in genes_seen:
        return genes_seen[gene]
    else:
        r = requests.get('http://togows.org/entry/kegg-genes/%s/dblinks.json' % gene)
        if r.status_code == 200:
            if len(r.json()) == 0:
                warnings.warn("No gene found with id %s" % gene)
                genes_seen[gene] = list()
                return list()
            else:
                try:
                    genes_seen[gene] = r.json()[0]['RN']
                    return r.json()[0]['RN']
                except KeyError:
                    genes_seen[gene] = list()
                    return list()
        else:
            warnings.warn("Connection to kegg via togows not able to be established.")
            return list()


def get_kegg_rxns_from_gene2(gene):
    global genes_seen
    if gene in genes_seen:
        return genes_seen[gene]
    else:
        r = requests.get('http://rest.kegg.jp/get/%s' % gene)
        if r.status_code == 200:
            if len(r.text) == 0:
                warnings.warn("No gene found with id %s" % gene)
                return list()
            else:
                rxn_line = [i for i in r.text.split('\n') if "RN:" in i.strip().split()]
                if len(rxn_line) == 1:
                    rxns = rxn_line[0][12:].split()[1:]
                    genes_seen[gene] = rxns
                    return rxns
                elif len(rxn_line) == 0:
                    return list()
                else:
                    warnings.warn("More than one reaction DBLINKS RN for %s" % gene)
                    return list()
        elif r.status_code == 404:
            warnings.warn("Connection to kegg not able to be established.")
            return list()
        else:
            warnings.warn("No gene found with id %s" % gene)
            return list()


rxns_seen = {}
def get_kegg_rxn(rxn):
    global rxns_seen
    if rxn in rxns_seen:
        return rxns_seen[rxn]
    else:
        r = requests.get('http://togows.org/entry/kegg-reaction/%s/equation.json' % rxn)
        if r.status_code == 200:
            if len(r.json()) == 0:
                warnings.warn("No reaction found with id %s" % rxn)
                rxns_seen[rxn] = list(), list(), False
                return list(), list(), False
            else:
                eq_str = r.json()[0]
                equ = eq_str.split('=>')
                if equ[0][-1] == '<':
                    rev = True
                    equ[0] = equ[0][:-1]
                else:
                    rev = False
                reacts = list()
                for part in equ[0].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        reacts.append(part[:6])
                prods = list()
                for part in equ[1].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        prods.append(part[:6])
                rxns_seen[rxn] = reacts, prods, rev
                return reacts, prods, rev
        else:
            warnings.warn("Connection to kegg via togows not able to be established")
            return list(), list(), False


def get_kegg_rxn2(rxn):
    global rxns_seen
    if rxn in rxns_seen:
        return rxns_seen[rxn]
    else:
        r = requests.get('http://rest.kegg.jp/get/%s' % rxn)
        if r.status_code == 200:
            if len(r.text) == 0:
                warnings.warn("No reaction found with id %s" % rxn)
                rxns_seen[rxn] = list(), list(), False
                return list(), list(), False
            else:
                equ_line = [i for i in r.text.split('\n') if i.startswith('EQUATION')][0]
                equ_str = ' '.join(equ_line.split()[1:])
                equ = equ_str.split('=>')
                if equ[0][-1] == '<':
                    rev = True
                    equ[0] = equ[0][:-1]
                else:
                    rev = False
                reacts = list()
                for part in equ[0].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        reacts.append(part[:6])
                prods = list()
                for part in equ[1].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        prods.append(part[:6])
                rxns_seen[rxn] = reacts, prods, rev
                return reacts, prods, rev
        else:
            warnings.warn("Connection to kegg via togows not able to be established")
            return list(), list(), False


def get_reactome(genome, threads=20):
    reactome = list()
    pool = ThreadPool(processes=threads)
    pool.map_async(get_kegg_rxns_from_gene, genome, callback=reactome.extend)
    pool.close()
    pool.join()
    return set([j for i in reactome for j in i])


def get_reactome2(genome, threads=20):
    reactome = list()
    pool = ThreadPool(processes=threads)
    pool.map_async(get_kegg_rxns_from_gene2, genome, callback=reactome.extend)
    pool.close()
    pool.join()
    return set([j for i in reactome for j in i])


def get_rxns(reactome, threads=20):
    rxns = list()
    pool = ThreadPool(processes=threads)
    pool.map_async(get_kegg_rxn, reactome, callback=rxns.append)
    pool.close()
    pool.join()
    return [j for i in rxns for j in i]


def get_rxns2(reactome, threads=20):
    rxns = list()
    pool = ThreadPool(processes=threads)
    pool.map_async(get_kegg_rxn2, reactome, callback=rxns.append)
    pool.close()
    pool.join()
    return [j for i in rxns for j in i]


def make_metabolic_network(rxns, filter_very_common=True, filter_common=False, only_giant=False,
                           min_component_size=None):
    metab_net = nx.DiGraph()
    for reacts, prods, rev in rxns:
        for react in reacts:
            if react not in metab_net.nodes():
                metab_net.add_node(react)
            for prod in prods:
                if prod not in metab_net.nodes():
                    metab_net.add_node(prod)
                if (react, prod) not in metab_net.edges():
                    metab_net.add_edge(react, prod)

    if filter_very_common:
        try:
            f = open("cos_to_remove.txt")
            nodes_to_remove = set([i.strip() for i in f.readlines()])
            metab_net.remove_nodes_from(nodes_to_remove)
        except IOError:
            warnings.warn("cos_to_remove.txt not found. Run determine_cos_to_remove.py to create.")

    if filter_common:
        for node, degree in metab_net.degree_iter():
            if degree > 10:  # picked 20 by looking at degree distribution of some otus
                metab_net.remove_node(node)

    if only_giant:
        components = list(nx.weakly_connected_components(metab_net))
        giant_component = set()
        for component in components:
            if len(component) > len(giant_component):
                giant_component = component
        metab_net.remove_nodes_from(set(metab_net.nodes()) - giant_component)

    if not only_giant and type(min_component_size) == int:
        components = list(nx.weakly_connected_components(metab_net))
        nodes_to_remove = list()
        for component in components:
            if len(component) < min_component_size:
                nodes_to_remove += list(component)
        metab_net.remove_nodes_from(nodes_to_remove)

    return metab_net


def determine_seed_set(metab_net):
    sscs = list(nx.strongly_connected_components(metab_net))
    seed_sets = defaultdict(list)
    seed_group = 0
    for ssc in sscs:
        in_nodes = set([i[0] for i in metab_net.in_edges(ssc)])
        if len(in_nodes - set(ssc)) == 0:
            for co in ssc:
                metab_net.node[co]['Seed'] = 1
                metab_net.node[co]['SeedGroup'] = seed_group
                seed_sets[seed_group].append(co)
            seed_group += 1
        else:
            for co in ssc:
                metab_net.node[co]['Seed'] = 0
    return metab_net, seed_sets


def calculate_bss(network1, seeds1, network2, seeds2):
    # calculate bss for otu1 relative to otu2
    overlap_seed1net2 = 0
    network2_nodes = set(network2.nodes())
    for seeds in seeds1.values():
        if len(set(seeds) & network2_nodes) > 0:
            overlap_seed1net2 += 1
    seed1net2_bss = overlap_seed1net2 / len(seeds1)

    # calculate bss for otu2 relative to otu1
    overlap_seed2net1 = 0
    network1_nodes = set(network1.nodes())
    for seeds in seeds2.values():
        if len(set(seeds) & network1_nodes) > 0:
            overlap_seed2net1 += 1
    seed2net1_bss = overlap_seed2net1 / len(seeds2)

    return seed1net2_bss, seed2net1_bss


def calculate_mci(network1, seeds1, network2, seeds2):
    # calculate mci for otu1 relative to otu2
    otu2_inseeds = set.union(*[set(i) for i in seeds2.values()])
    otu2_nodes = set(network2.nodes())
    overlap_seed1seed2 = 0
    for seeds in seeds1.values():
        if len(set(seeds) & otu2_nodes) > 0:  # check if any of seed group in other network
            if len(set(seeds) & otu2_inseeds) == 0:  # check if any of seed group is seed in other network
                overlap_seed1seed2 += 1
    seed1net2_mci = overlap_seed1seed2 / len(seeds1)

    # calculate mci for otu2 relative to otu1
    otu1_inseeds = set.union(*[set(i) for i in seeds1.values()])
    otu1_nodes = set(network1.nodes())
    overlap_seed2seed1 = 0
    for seeds in seeds2.values():
        if len(set(seeds) & otu1_nodes) > 0:
            if len(set(seeds) & otu1_inseeds) == 0:
                overlap_seed2seed1 += 1
    seed2net1_mci = overlap_seed2seed1 / len(seeds2)

    return seed1net2_mci, seed2net1_mci
