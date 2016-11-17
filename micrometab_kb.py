from flask import Flask, render_template, request, redirect, url_for, flash
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from database_setup import Genome, Base
import networkx as nx
import json
from py2cytoscape import util as cy
import metabolic_network_analysis as mna

app = Flask(__name__)

engine = create_engine('sqlite:///gg_genomes.db')
Base.metadata.bind = engine

DBSession = sessionmaker(bind=engine)
session = DBSession()

# TODO: Make taxa strings pretty
# TODO: Figure out a way to add in tree distance between OTUs
# TODO: Pathway enrichment of seeds
# TODO: Product seed relationship
# TODO: Regression of 16S similarity vs competition and cooperation
# TODO: Is compound bacterial made only
# TODO: show COs that go into the scores, what are complements and what are they competeing over
# TODO: add in top bar to link back


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
                flash("OTU id %s not in database." % request.form['name'])
                return redirect(url_for('welcome_page'))
            metab_net_json = json.loads(genome.metab_net)
            metab_net = cy.to_networkx(metab_net_json)
            metab_net, ss = mna.determine_seed_set(metab_net)
            seeds = [j for i in ss.values() for j in i]
            return render_template('singleOTUResult.html', genome=genome, taxa_str=pretty_taxa(genome.taxonomy), seeds=sorted(seeds),
                                   eles=json.dumps(metab_net_json['elements']))
        else:
            flash("No OTU id entered for single analysis.")
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

            # get data and render webpage
            metab_net1_json = json.loads(genome1.metab_net)
            metab_net1 = cy.to_networkx(metab_net1_json)
            metab_net2_json = json.loads(genome2.metab_net)
            metab_net2 = cy.to_networkx(metab_net2_json)
            metab_net1, ss1 = mna.determine_seed_set(metab_net1)
            metab_net2, ss2 = mna.determine_seed_set(metab_net2)
            seeds1 = set([j for i in ss1.values() for j in i])
            seeds2 = set([j for i in ss2.values() for j in i])
            seeds1_only = seeds1-seeds2
            seeds2_only = seeds2-seeds1
            shared_seeds = seeds1 & seeds2
            net1net2_bss, net2net1_bss = mna.calculate_bss(metab_net1, ss1, metab_net2, ss2)
            net1net2_mci, net2net1_mci = mna.calculate_mci(metab_net1, ss1, metab_net2, ss2)
            return render_template('pairOTUResult.html', genome1=genome1, taxa_str1=pretty_taxa(genome1.taxonomy),
                                   seeds1=sorted(seeds1_only), eles1=json.dumps(metab_net1_json['elements']),
                                   genome2=genome2, taxa_str2=pretty_taxa(genome2.taxonomy), seeds2=sorted(seeds2_only),
                                   eles2=json.dumps(metab_net2_json['elements']), shared_seeds=sorted(shared_seeds),
                                   net1net2_bss=round(net1net2_bss, 2), net2net1_bss=round(net2net1_bss, 2),
                                   net1net2_mci=round(net1net2_mci, 2), net2net1_mci=round(net2net1_mci, 2))
        else:
            flash("Need to enter two OTU id's to compare OTUs")
            return redirect(url_for('welcome_page'))
    else:
        flash("How did you get to pairs then?")
        return redirect(url_for('welcome_page'))


@app.route('/result/single_otu/<int:otu_id>/network/')
def generate_cyto_network(otu_id):
    genome = session.query(Genome).filter_by(name=otu_id).one()
    metab_net = json.loads(genome.metab_net)
    return render_template('simple_cyto.html', eles=json.dumps(metab_net['elements']))


@app.route('/simple/')
def generate_simple_network():
    house_net = nx.house_graph()
    house_net_json = cy.from_networkx(house_net)
    return render_template("simple_cyto.html", eles=json.dumps(house_net_json['elements']))

if __name__ == '__main__':
    app.secret_key = 'super_secret_key'
    app.debug = True
    app.run()
