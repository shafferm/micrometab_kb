from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from database_setup import Base, Genome
import multiprocessing
import networkx as nx
import gzip
from os import path
from picrust.util import get_picrust_project_dir, convert_precalc_to_biom
from parse_KEGG import KEGG_Parser
import json
from py2cytoscape import util as cy

engine = create_engine('sqlite:///gg_genomes.db')
Base.metadata.bind = engine
DBSession = sessionmaker(bind=engine)

session = DBSession()

kegg = KEGG_Parser()

GG_LOC = "/Users/shafferm/lab/HIV_5runs/qiime_files/99_otu_taxonomy.txt"
chunk_size = 2000
procs = 3


def load_data_table(ids_to_load):
    """Stolen from https://github.com/picrust/picrust/blob/master/scripts/predict_metagenomes.py
    and modified.
    Load a data table, detecting gziiped files and subset loading
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


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def generate_genome(chunk):
    if len(chunk) == 1:
        genome_table = load_data_table(list(chunk[0][0]))
    else:
        genome_table = load_data_table([i[0] for i in chunk])
    for otu in chunk:
        otu_id = otu[0]
        taxonomy = otu[1]
        nsti = genome_table.metadata(otu_id)['NSTI']
        genome = genome_table.ids(axis="observation")[genome_table.data(otu_id) > 0]
        genome = [str(i) for i in genome]
        metab_network = make_metabolic_network(genome, min_component_size=4)
        metab_network_json = cy.from_networkx(metab_network)
        # metab_network, seed_sets = mna.determine_seed_set(metab_network)
        # ggGenome = Genome(name=otu_id, nsti=NSTI, metab_net=str(metab_network), seed_sets=str(seed_sets),
        #                  genome=str(genome), taxonomy=taxonomy)
        genome = Genome(name=int(otu_id), nsti=float(nsti), metab_net=json.dumps(metab_network_json),
                        genome=str(genome), taxonomy=taxonomy)
        session.add(genome)
        session.commit()
        print otu_id


def main():
    gg_genomes = {i.strip().split('\t')[0]: i.strip().split('\t')[1] for i in open(GG_LOC).readlines()}

    # pool = multiprocessing.Pool(procs)

    # for chunk in chunks(gg_genomes.items(), chunk_size):
    #     pool.map(generate_genome, chunk)
    #     # potentially have return objects and add to database not concurrently with a loop here
    #     pool.close()
    #     pool.join()

    for chunk in chunks(gg_genomes.items(), chunk_size):
        generate_genome(chunk)


if __name__ == "__main__":
    main()
