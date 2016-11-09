from flask import Flask, render_template, request, redirect, url_for, flash, jsonify
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from database_setup import Genome, Base
import networkx as nx
import json
import uuid
from py2cytoscape import util as cy

app = Flask(__name__)

engine = create_engine('sqlite:///gg_genomes.db')
Base.metadata.bind = engine

DBSession = sessionmaker(bind=engine)
session = DBSession()


@app.route('/')
def welcome_page():
    return render_template('index.html')


@app.route('/single_otu/', methods=['GET', 'POST'])
def retrive_single_otu_entry():
    if request.method == 'POST':
        otu_id = request.form['name']
        return redirect(url_for('single_otu_result', otu_id=otu_id))
    else:
        return render_template('singleOTURequest.html', otu_ids=[str(i.name) for i in session.query(Genome).all()])


@app.route('/result/single_otu/<int:otu_id>/')
def single_otu_result(otu_id):
    genome = session.query(Genome).filter_by(name=otu_id).one()
    return render_template('singleOTUResult.html', genome=genome)


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

if __name__ == '__main__':
    app.debug = True
    app.run()
