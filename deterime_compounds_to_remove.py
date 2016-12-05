from collections import Counter

from micrometab_analysis import parse_KEGG
from micrometab_analysis.metabolic_network_analysis import make_metabolic_network

kegg = parse_KEGG.KEGG_Parser()

rxs = parse_KEGG.get_reactions()
all_cos = [k for i in list(rxs.values()) for j in i[:1] for k in j]
co_counter = Counter(all_cos)
net = make_metabolic_network(list(parse_KEGG.get_ko_names().keys()), filter_common=False)
cos_to_remove = [node for node, deg in net.degree().items() if deg >= 30]
with open("cos_to_remove.txt", 'w') as f:
    f.write('\n'.join(cos_to_remove))
    f.write('\n')
