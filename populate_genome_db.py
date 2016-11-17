import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from database_setup import Base, Genome
import multiprocessing
import json
from py2cytoscape import util as cy
from metabolic_network_analysis import load_data_table, make_metabolic_network

engine = create_engine('sqlite:///gg_genomes.db')
Base.metadata.bind = engine
DBSession = sessionmaker(bind=engine)

session = DBSession()


GG_LOC = "/Users/shafferm/lab/HIV_5runs/qiime_files/99_otu_taxonomy.txt"
chunk_size = 10
procs = 3


def breakup_list(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def generate_genome(otus):
    genome_table = load_data_table([i[0] for i in otus])
    genomes = list()
    for otu_id, taxonomy in otus:
        nsti = genome_table.metadata(otu_id)['NSTI']
        genome = genome_table.ids(axis="observation")[genome_table.data(otu_id) > 0]
        genome = [str(i) for i in genome]
        metab_network = make_metabolic_network(genome)
        metab_network_json = cy.from_networkx(metab_network)
        genome = Genome(name=int(otu_id), nsti=float(nsti), metab_net=json.dumps(metab_network_json),
                        genome=str(genome), taxonomy=taxonomy)
        genomes.append(genome)
    return genomes


def add_to_db(chunks):
    for chunk in chunks:
        for genome in chunk:
            session.add(genome)
            session.commit()
            print genome.name


def main(args):
    gg_genomes = {i.strip().split('\t')[0]: i.strip().split('\t')[1] for i in open(args.gg_loc).readlines()[:100]}

    chunks = breakup_list(gg_genomes.items(), chunk_size)
    pool = multiprocessing.Pool(args.nprocs)
    pool.map_async(generate_genome, chunks, callback=add_to_db)
    pool.close()
    pool.join()

    # for chunk in chunks(gg_genomes.items(), chunk_size):
    #     generate_genome(chunk)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--gg_loc", help="greengenes taxonomy location", default=GG_LOC)
    parser.add_argument("--nprocs", help="number of processors", default=procs)
    args = parser.parse_args()
    main(args)
