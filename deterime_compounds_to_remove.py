import parse_KEGG
from collections import Counter
import networkx as nx

kegg = parse_KEGG.KEGG_Parser()

def make_metabolic_network(genome, filter_common=True, only_giant=False, min_component_size=None):
    # TODO: Add common co filter to make less dense networks
    metab_net = nx.DiGraph()
    for gene in genome:
        rxns = kegg.get_rxns_from_ko(gene)
        if rxns != set():
            for rxn in rxns:
                reacts, prods, rev = kegg.get_rxn(rxn)
                for react in reacts:
                    if react not in metab_net.nodes():
                        metab_net.add_node(react)
                    for prod in prods:
                        if prod not in metab_net.nodes():
                            metab_net.add_node(prod)
                        if (react, prod) not in metab_net.edges():
                            metab_net.add_edge(react, prod)
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

rxs = parse_KEGG.get_reactions()
all_cos = [k for i in rxs.values() for j in i[:1] for k in j]
co_counter = Counter(all_cos)
net = make_metabolic_network(parse_KEGG.get_ko_names().keys(), filter_common=False)
cos_to_remove = [node for node, deg in net.degree().iteritems() if deg >= 30]
with open("cos_to_remove.txt", 'w') as f:
    f.write('\n'.join(cos_to_remove))
    f.write('\n')
