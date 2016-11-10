from flask import Flask, render_template, request, redirect, url_for, flash, jsonify
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from database_setup import Genome, Base
import networkx as nx
import json
import uuid
from py2cytoscape import util as cy
import metabolic_network_analysis as mna

app = Flask(__name__)

engine = create_engine('sqlite:///gg_genomes.db')
Base.metadata.bind = engine

DBSession = sessionmaker(bind=engine)
session = DBSession()

# TODO: Make taxa strings pretty
# TODO: Figure out a way to add in tree distance between OTUs


@app.route('/')
def welcome_page():
    return render_template('index.html')


@app.route('/single_otu/', methods=['GET', 'POST'])
def retrive_single_otu_entry():
    if request.method == 'POST':
        otu_id = request.form['name']
        try:
            session.query(Genome).filter_by(name=otu_id).one()
            return redirect(url_for('single_otu_result', otu_id=otu_id))
        except NoResultFound:
            flash("GreenGenes OTU ID %s not found in database." % otu_id)
            return render_template('singleOTURequest.html')
    else:
        return render_template('singleOTURequest.html')


@app.route('/result/single_otu/<int:otu_id>/')
def single_otu_result(otu_id):
    genome = session.query(Genome).filter_by(name=otu_id).one()
    metab_net = cy.to_networkx(json.loads(genome.metab_net))
    metab_net, ss = mna.determine_seed_set(metab_net)
    seeds = [j for i in ss.values() for j in i]
    return render_template('singleOTUResult.html', genome=genome, seeds=sorted(seeds))


@app.route('/pair_otu/', methods=['GET', 'POST'])
def retrive_paired_otu_entry():
    if request.method == 'POST':
        exception_caught = False
        # check and make sure otu1 is in the database
        otu1_id = request.form['name1']
        try:
            session.query(Genome).filter_by(name=otu1_id).one()
        except NoResultFound:
            exception_caught = True
            flash("GreenGenes OTU ID %s not found in database." % otu1_id)

        # check and make sure otu2 is in the database
        otu2_id = request.form['name2']
        try:
            session.query(Genome).filter_by(name=otu2_id).one()
        except NoResultFound:
            exception_caught = True
            flash("GreenGenes OTU ID %s not found in database." % otu2_id)

        # if either OTU id entered is not in database then go back to request screen
        if exception_caught:
            return render_template('pairOTURequest.html')
        else:
            return redirect(url_for('pair_otu_result', otu1_id=otu1_id, otu2_id=otu2_id))

    else:
        return render_template('pairOTURequest.html')


@app.route('/result/pair_otu/<int:otu1_id>/<int:otu2_id>/')
def pair_otu_result(otu1_id, otu2_id):
    genome1 = session.query(Genome).filter_by(name=otu1_id).one()
    genome2 = session.query(Genome).filter_by(name=otu2_id).one()
    metab_net1 = cy.to_networkx(json.loads(genome1.metab_net))
    metab_net2 = cy.to_networkx(json.loads(genome2.metab_net))
    metab_net1, ss1 = mna.determine_seed_set(metab_net1)
    metab_net2, ss2 = mna.determine_seed_set(metab_net2)
    seeds1 = set([j for i in ss1.values() for j in i])
    seeds2 = set([j for i in ss2.values() for j in i])
    seeds1_only = seeds1-seeds2
    seeds2_only = seeds2-seeds1
    shared_seeds = seeds1 & seeds2
    net1net2_bss, net2net1_bss = mna.calculate_bss(metab_net1, ss1, metab_net2, ss2)
    net1net2_mci, net2net1_mci = mna.calculate_mci(metab_net1, ss1, metab_net2, ss2)
    return render_template('pairOTUResult.html', genome1=genome1, seeds1=sorted(seeds1_only), genome2=genome2,
                           seeds2=sorted(seeds2_only), shared_seeds=sorted(shared_seeds),
                           net1net2_bss=round(net1net2_bss, 2), net2net1_bss=round(net2net1_bss, 2),
                           net1net2_mci=round(net1net2_mci, 2), net2net1_mci=round(net2net1_mci, 2))


# Network testing
@app.route('/result/single_otu/<int:otu_id>/network_JSON/')
def generate_network_json(otu_id):
    genome = session.query(Genome).filter_by(name=otu_id).one()
    metab_network_json = json.loads(genome.metab_net)
    # metab_network_json['links'] = [{'source': metab_network_json['nodes'][link['source']]['id'],
    #                                 'target': metab_network_json['nodes'][link['target']]['id']}
    #                                for link in metab_network_json['links']]
    return jsonify(metab_network_json)


@app.route('/result/single_otu/<int:otu_id>/network/')
def generate_network(otu_id):
    return render_template('viewNetwork.html', otu_id=otu_id)


@app.route('/test_network/')
def generate_test_network():
    return render_template('viewTestNetwork.html')


@app.route('/house_network/')
def generate_house_network():
    house_net = nx.house_graph()
    house_net_json = cy.from_networkx(house_net)
    new_uuid = "cy" + str(uuid.uuid4())
    width = 1098
    height = 700
    nodes = house_net_json['elements']['nodes']
    edges = house_net_json['elements']['edges']
    return render_template('viewCytoscapeNetwork.html', nodes=json.dumps(nodes), edges=json.dumps(edges), uuid=new_uuid,
                           widget_width=str(width), widget_height=str(height))


@app.route('/result/single_otu/<int:otu_id>/cyto_network/')
def generate_cyto_network(otu_id):
    genome = session.query(Genome).filter_by(name=otu_id).one()
    metab_net = json.loads(genome.metab_net)
    nodes = metab_net['elements']['nodes']
    edges = metab_net['elements']['edges']
    new_uuid = "cy" + str(uuid.uuid4())
    width = 1098
    height = 700
    return render_template('viewCytoscapeNetwork.html', nodes=json.dumps(nodes), edges=json.dumps(edges), uuid=new_uuid,
                           widget_width=str(width), widget_height=str(height))


if __name__ == '__main__':
    app.secret_key = 'super_secret_key'
    app.debug = True
    app.run()
