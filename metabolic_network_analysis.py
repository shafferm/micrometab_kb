from __future__ import division

import networkx as nx
from collections import defaultdict


def determine_seed_set(metab_net):
    sscs = list(nx.strongly_connected_components(metab_net))
    seed_sets = defaultdict(list)
    seed_group = 0
    for ssc in sscs:
        in_nodes = set([i[0] for i in metab_net.in_edges(ssc)])
        if len(in_nodes - ssc) == 0:
            for co in ssc:
                metab_net.node[co]['Seed'] = "True"
                metab_net.node[co]['SeedGroup'] = seed_group
                seed_sets[seed_group].append(co)
            seed_group += 1
        else:
            for co in ssc:
                metab_net.node[co]['Seed'] = "False"
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
