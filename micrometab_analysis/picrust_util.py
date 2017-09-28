from os import path
from biom.table import Table
import numpy as np
from io import StringIO
import gzip
"""Functions stolen from picrust, adapted to python 3 and then adapted for my usage"""


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
    precalc_file_name = '_'.join(["ko", "13_5", 'precalculated.tab.gz'])
    data_table_fp = path.join(get_data_dir(), precalc_file_name)

    if not path.exists(data_table_fp):
        raise IOError("File " + data_table_fp + " doesn't exist! Did you forget to download it?")

    genome_table_fh = gzip.open(data_table_fp, 'rt')
    genome_table = convert_precalc_to_biom(genome_table_fh, ids_to_load)
    return genome_table


def get_data_dir():
    """ Returns the top-level PICRUST directory
    """
    # Get the full path of util.py
    current_file_path = path.abspath(__file__)
    # Get the directory containing util.py
    current_dir_path = path.dirname(current_file_path)
    # Return the directory containing the directory containing util.py
    return path.join(current_dir_path, 'data')


def determine_metadata_type(line):
    if ';' in line:
        if '|' in line:
            return 'list_of_lists'
        else:
            return 'list'
    else:
        return 'string'


def parse_metadata_field(metadata_str,metadata_format='string'):
    if metadata_format == 'string':
        return metadata_str
    elif metadata_format == 'list':
        return [e.strip() for e in metadata_str.split(';')]
    elif metadata_format == 'list_of_lists':
        return [[e.strip() for e in y.split(';')] for y in metadata_str.split('|')]


def convert_precalc_to_biom(precalc_in, ids_to_load=None,transpose=True,md_prefix='metadata_'):
    """Loads PICRUSTs tab-delimited version of the precalc file and outputs a BIOM object"""

    #if given a string convert to a filehandle
    if type(precalc_in) == str:
        fh = StringIO(precalc_in)
    else:
        fh=precalc_in

    #first line has to be header
    header_ids=fh.readline().strip().split('\t')

    col_meta_locs={}
    for idx,col_id in enumerate(header_ids):
        if col_id.startswith(md_prefix):
            col_meta_locs[col_id[len(md_prefix):]]=idx

    end_of_data=len(header_ids)-len(col_meta_locs)
    trait_ids = header_ids[1:end_of_data]

    col_meta=[]
    row_meta=[{} for i in trait_ids]

    if ids_to_load is not None and len(ids_to_load) > 0:
        ids_to_load=set(ids_to_load)
        load_all_ids=False
    else:
        load_all_ids=True

    matching=[]
    otu_ids=[]
    for line in fh:
        fields = line.strip().split('\t')
        row_id=fields[0]
        if(row_id.startswith(md_prefix)):
            #handle metadata

            #determine type of metadata (this may not be perfect)
            metadata_type=determine_metadata_type(line)
            for idx,trait_name in enumerate(trait_ids):
                row_meta[idx][row_id[len(md_prefix):]]=parse_metadata_field(fields[idx+1],metadata_type)

        elif load_all_ids or (row_id in set(ids_to_load)):
            otu_ids.append(row_id)
            matching.append(list(map(float,fields[1:end_of_data])))

            #add metadata
            col_meta_dict={}
            for meta_name in col_meta_locs:
                col_meta_dict[meta_name]=fields[col_meta_locs[meta_name]]
            col_meta.append(col_meta_dict)

            if not load_all_ids:
                ids_to_load.remove(row_id)

    if not otu_ids:
        raise ValueError("No OTUs match identifiers in precalculated file. PICRUSt requires an OTU table reference/closed picked against GreenGenes.\nExample of the first 5 OTU ids from your table: {0}".format(', '.join(list(ids_to_load)[:5])))

    if ids_to_load:
       raise ValueError("One or more OTU ids were not found in the precalculated file!\nAre you using the correct --gg_version?\nExample of (the {0}) unknown OTU ids: {1}".format(len(ids_to_load),', '.join(list(ids_to_load)[:5])))

    #note that we transpose the data before making biom obj
    matching = np.asarray(matching)
    if transpose:
        return Table(matching.T, trait_ids, otu_ids, row_meta, col_meta,
                     type='Gene table')
    else:
        return Table(matching, otu_ids, trait_ids, col_meta, row_meta,
                     type='Gene table')
