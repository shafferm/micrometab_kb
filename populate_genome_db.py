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
chunk_size = 2000
procs = 3


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
