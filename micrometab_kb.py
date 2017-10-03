import json
from os import path
from itertools import zip_longest

from flask import Flask, render_template, request, redirect, url_for, flash, jsonify
from py2cytoscape import util as cy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from skbio import TreeNode

from database_setup import Genome, Base
from micrometab_analysis import metabolic_network_analysis as mna
from micrometab_analysis import picrust_util as pu

app = Flask(__name__)

engine = create_engine('sqlite:///gg_genomes.db')
Base.metadata.bind = engine

DBSession = sessionmaker(bind=engine)
session = DBSession()

# TODO: Figure out a way to add in tree distance between OTUs
# TODO: Pathway enrichment of seeds
# TODO: Product seed relationship
# TODO: Regression of 16S similarity vs competition and cooperation
# TODO: Is compound bacterial made only
# TODO: make a setup.py that downloads greengenes (and picrust?) and installs dependencies
# TODO: add human to database, other host organisms?
# TODO: add human vs bacterial or both generated
# TODO: extract micrometab_analysis to it's own package
# TODO: add tooltip explaining MCI and BSS


def pretty_taxa(taxa_str):
    taxa = [i[3:] for i in taxa_str.split('; ') if len(i[3:]) > 0]
    if len(taxa) > 5:
        return ' '.join(taxa[5:])
    elif len(taxa) == 5:
        return "%s family" % taxa[-1]
    elif len(taxa) == 4:
        return "%s order" % taxa[-1]
    elif len(taxa) == 3:
        return "%s class" % taxa[-1]
    elif len(taxa) == 2:
        return "%s phylum" % taxa[-1]
    elif len(taxa) == 1:
        return "%s kingdom" % taxa[-1]
    else:
        return "Unclassified"


# TODO: put this somewhere and have tree read in as a global
def get_tip2tip(otu1, otu2):
    tree_loc = path.join(pu.get_data_dir(), 'gg_13_8_otus', 'trees', '99_otus.tree')
    tree = TreeNode.read(tree_loc)
    tip_a = tree.find(str(otu1))
    tip_b = tree.find(str(otu2))
    return tip_a.distance(tip_b)


@app.route('/')
def welcome_page():
    return render_template('index.html')


@app.route('/result/single_otu/', methods=['GET', 'POST'])
def single_otu_result():
    if request.method == 'POST':
        if request.form['name']:
            try:
                genome = session.query(Genome).filter_by(name=request.form['name']).one()
            except NoResultFound:
                flash("OTU ID %s not in database." % request.form['name'])
                return redirect(url_for('welcome_page'))
            metab_net_json = json.loads(genome.metab_net)
            metab_net = cy.to_networkx(metab_net_json)
            metab_net, ss = mna.determine_seed_set(metab_net)
            seeds = [j for i in list(ss.values()) for j in i]
            return render_template('singleOTUResult.html', genome=genome, taxa_str=pretty_taxa(genome.taxonomy),
                                   seeds=sorted(seeds), eles=json.dumps(metab_net_json['elements']))
        else:
            flash("No OTU ID entered for single analysis.")
            return redirect(url_for('welcome_page'))
    else:
        flash("How did you get here to single?")
        return redirect(url_for('welcome_page'))


@app.route('/result/pair_otu/', methods=['GET', 'POST'])
def pair_otu_result():
    if request.method == 'POST':
        if request.form['name1'] and request.form['name2']:
            # get genomes from database
            exception = False
            genome1 = None
            try:
                genome1 = session.query(Genome).filter_by(name=request.form['name1']).one()
            except NoResultFound:
                flash("OTU %s not found in the database." % request.form['name1'])
                exception = True
            genome2 = None
            try:
                genome2 = session.query(Genome).filter_by(name=request.form['name2']).one()
            except NoResultFound:
                flash("OTU %s not found in the database." % request.form['name2'])
                exception = True
            if exception:
                return redirect(url_for('welcome_page'))

            # get tree data
            tip2tip = get_tip2tip(genome1.name, genome2.name)

            # get data and determine seeds
            metab_net1_json = json.loads(genome1.metab_net)
            metab_net1 = cy.to_networkx(metab_net1_json)
            metab_net2_json = json.loads(genome2.metab_net)
            metab_net2 = cy.to_networkx(metab_net2_json)
            metab_net1, ss1 = mna.determine_seed_set(metab_net1)
            metab_net2, ss2 = mna.determine_seed_set(metab_net2)
            seeds1 = set([j for i in list(ss1.values()) for j in i])
            seeds2 = set([j for i in list(ss2.values()) for j in i])

            # analyze seeds and determine metabolic characteristics
            seeds1_only = seeds1-seeds2
            if seeds1_only == set():
                seeds1_only = [None]
            seeds2_only = seeds2-seeds1
            if seeds2_only == set():
                seeds2_only = [None]
            shared_seeds = seeds1 & seeds2
            otu1_seeds_otu2_complement = seeds1_only & set(metab_net2.nodes())
            otu2_seeds_otu1_complement = seeds2_only & set(metab_net1.nodes())
            net1net2_bss, net2net1_bss = mna.calculate_bss(metab_net1, ss1, metab_net2, ss2)
            net1net2_mci, net2net1_mci = mna.calculate_mci(metab_net1, ss1, metab_net2, ss2)

            # render page
            return render_template('pairOTUResult.html', genome1=genome1, taxa_str1=pretty_taxa(genome1.taxonomy),
                                   seeds1=sorted(seeds1_only), eles1=json.dumps(metab_net1_json['elements']),
                                   genome2=genome2, taxa_str2=pretty_taxa(genome2.taxonomy), seeds2=sorted(seeds2_only),
                                   eles2=json.dumps(metab_net2_json['elements']), tip2tip=round(tip2tip, 2),
                                   shared_seeds=sorted(shared_seeds), net1net2_bss=round(net1net2_bss, 2),
                                   net2net1_bss=round(net2net1_bss, 2), net1net2_mci=round(net1net2_mci, 2),
                                   net2net1_mci=round(net2net1_mci, 2),
                                   otu1_seeds_otu2_complement=otu1_seeds_otu2_complement,
                                   otu2_seeds_otu1_complement=otu2_seeds_otu1_complement)
        else:
            flash("Need to enter two OTU ID's to compare OTUs")
            return redirect(url_for('welcome_page'))
    else:
        flash("How did you get to pairs then?")
        return redirect(url_for('welcome_page'))


@app.route('/get/<string:otu_ids>')
def get_otu_json(otu_ids):
    otu_ids = otu_ids.split(',')
    genome_dict = dict()
    for i in range(0, len(otu_ids), 995):
        chunk = otu_ids[i:i+995]
        genomes = session.query(Genome).filter(Genome.name.in_(chunk))
        if genomes.count() == len(chunk):
            genome_dict.update({genome.name: genome.serialize for genome in genomes})
        else:
            found_ids = set(otu_ids) - set([i.name for i in genomes])
            raise NoResultFound("Not all OTUs found. %s are missing" % ', '.join(found_ids))
    return jsonify(genome_dict)


if __name__ == '__main__':
    app.secret_key = 'super_secret_key'
    app.debug = True
    app.run()
