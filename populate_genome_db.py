import argparse
import json
import multiprocessing
import os
from datetime import datetime

import requests
from py2cytoscape import util as cy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from database_setup import Base, Genome
from micrometab_analysis import metabolic_network_analysis as mna
from micrometab_analysis.picrust_util import load_data_table

engine = create_engine('sqlite:///gg_genomes.db')
Base.metadata.bind = engine
DBSession = sessionmaker(bind=engine)

session = DBSession()


GG_LOC = "/Users/shafferm/lab/HIV_5runs/qiime_files/99_otu_taxonomy.txt"
KEGG_LOC = "/Users/shafferm/KEGG_late_june2011_snapshot/"
chunk_size = 10000
procs = 3


# TODO: download greengenes and picrust files when called if not already present
def download_file(url, local_directory):
    """
    Adapted from http://stackoverflow.com/questions/16694907/how-to-download-large-file-in-python-with-requests-py
    Parameters
    ----------
    url: url to download from
    local_directory: location to download

    Returns
    -------
    local_path: full file path of downloaded file
    """
    local_filename = url.split('/')[-1]
    local_path = os.path.join(local_directory, local_filename)
    r = requests.get(url, stream=True)
    with open(local_path, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
    return local_path


def breakup_list(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def generate_genome_local(otus, loc=None):
    genome_table = load_data_table([i[0] for i in otus])
    genomes = list()
    for otu_id, taxonomy in otus:
        print(os.getpid(), otu_id)
        nsti = genome_table.metadata(otu_id)['NSTI']
        genome = genome_table.ids(axis="observation")[genome_table.data(otu_id) > 0]
        genome = [str(i) for i in genome]
        print(os.getpid(), otu_id, "genome length", len(genome))
        reactome = mna.get_reactome_local(genome, loc)
        print(os.getpid(), otu_id, "reactome length", len(reactome))
        rxns = mna.get_rxns_local(reactome, loc)
        print(os.getpid(), otu_id, "reaction count", len(rxns))
        metab_network = mna.make_metabolic_network(rxns, only_giant=True)
        print(os.getpid(), otu_id, "network made")
        metab_network_json = cy.from_networkx(metab_network)
        genome = Genome(name=int(otu_id), nsti=float(nsti), metab_net=json.dumps(metab_network_json),
                        genome=','.join(genome), taxonomy=taxonomy)
        genomes.append(genome)
    return genomes


def add_chunk_to_db(chunk):
    for genome in chunk:
        session.add(genome)
        session.commit()
        print(genome.name)


def add_chunks_to_db(chunks):
    for chunk in chunks:
        for genome in chunk:
            session.add(genome)
            session.commit()
            print(genome.name)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--gg_loc", help="greengenes taxonomy location", default=GG_LOC)
    parser.add_argument("--database_loc", help="location of kegg database", default=KEGG_LOC)
    parser.add_argument("--nprocs", help="number of processors", type=int, default=procs)
    parser.add_argument("--chunk_size", help="size of chunks to analyze per processor", type=int, default=chunk_size)
    parser.add_argument("--subset", help="size of subset of gg_genomes to analyze", type=int)
    args = parser.parse_args()

    start = datetime.now()

    if args.subset is None:
        gg_genomes = {i.strip().split('\t')[0]: i.strip().split('\t')[1] for i in open(args.gg_loc).readlines()}
    else:
        gg_genomes = {i.strip().split('\t')[0]: i.strip().split('\t')[1]
                      for i in open(args.gg_loc).readlines()[:args.subset]}

    chunks = breakup_list(list(gg_genomes.items()), args.chunk_size)
    pool = multiprocessing.Pool(args.nprocs)
    # pool.map_async(generate_genome, chunks, callback=add_chunks_to_db)
    for chunk in chunks:
        pool.apply_async(generate_genome_local, (chunk, args.database_loc), callback=add_chunk_to_db)
    pool.close()
    pool.join()

    finish = datetime.now()
    print(finish-start)


if __name__ == "__main__":
    main()
