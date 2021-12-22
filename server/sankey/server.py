from flask import Flask, request, jsonify
from flask_cors import CORS

import argparse
from jinja2 import Template

import scanpy as sc
import pandas as pd

app = Flask(__name__)
CORS(app)

@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route('/api/data')
def data():
   adata = app.config['adata']
   timepoint = app.config['timepoint']
   phenotype = app.config['phenotype']
   clone = app.config['clone']

   umap = pd.DataFrame(adata.obsm['X_umap'])
   umap.columns = ['UMAP_1', 'UMAP_2']

   df = adata.obs[[timepoint, phenotype, clone]]
   df = df.reset_index()
   df = df.rename(columns={'index': 'cell_id'})

   df = df.merge(umap, left_index=True, right_index=True)

   response = jsonify(df.to_dict(orient='records'))
   response.headers.add("Access-Control-Allow-Origin", "*")

   return response

@app.route('/api/cells')
def cells():
    timepoints = request.args.get('timepoint', []).split(',')
    phenotypes = request.args.get('phenotype', []).split(',')
    clones = request.args.get('clone', None)

    adata = app.config['adata']
    clone = app.config['clone']

    nodes = [f"{a[0]}_{a[1]}" for a in zip(timepoints, phenotypes)]

    df = adata.obs[adata.obs['viz_node'].isin(nodes)]
    
    if clones is not None:
        clones = clones.split(',')
        df = df[df[clone].isin(clones)]

    df = df.reset_index()
    df = df.rename(columns={'index': 'cell_id'})

    response = jsonify(df['cell_id'].tolist())
    response.headers.add("Access-Control-Allow-Origin", "*")

    return response


@app.route('/api/degenes')
def degenes():
    timepoints = request.args.get('timepoint', []).split(',')
    phenotypes = request.args.get('phenotype', []).split(',')
    clones = request.args.get('clone', None)

    adata = app.config['adata']
    clone = app.config['clone']

    nodes = [f"{a[0]}_{a[1]}" for a in zip(timepoints, phenotypes)]

    if clones is None:
        adata.obs['included'] = adata.obs.apply(lambda x: "included" if x['viz_node'] in nodes else "excluded", axis=1)
    else:
        clones = clones.split(',')
        adata.obs['included'] = adata.obs.apply(lambda x: "included" if x['viz_node'] in nodes and x[clone] in clones else "excluded", axis=1)

    sc.tl.rank_genes_groups(adata,"included")

    genes = pd.DataFrame(adata.uns['rank_genes_groups']["names"])[['included']].rename(columns={'included': 'gene'})
    adjpvals = pd.DataFrame(adata.uns['rank_genes_groups']["pvals_adj"])[['included']].rename(columns={'included': 'p'})
    logfc = pd.DataFrame(adata.uns['rank_genes_groups']["logfoldchanges"])[['included']].rename(columns={'included': 'fc'})
    df = logfc.merge(genes, left_index=True, right_index=True).merge(adjpvals, left_index=True, right_index=True)
    df = df[df['p'] < 0.05]
    df = df[df['fc'] > 0.25]
    final = df.to_dict(orient='records')

    final = sorted(final, key=lambda x: x['p'])

    response = jsonify(final)
    response.headers.add("Access-Control-Allow-Origin", "*")

    return response

def configure_app(path, timepoint, phenotype, clone, order, threshold):

    assert threshold > 0, 'Input threshold of at least 1'

    adata = open_data(path, timepoint, phenotype, clone)

    timepoints = get_timepoints(adata, order, timepoint)

    ## Filter data
    adata_filtered = filter_adata(adata, timepoint, phenotype, clone, timepoints, threshold)

    ## Add extra column for api
    adata_filtered.obs['viz_node'] = adata_filtered.obs[timepoint].astype(str) + "_" + adata_filtered.obs[phenotype].astype(str)


    app.config.update(
        adata=adata_filtered,
        timepoints=timepoints,
        threshold=threshold,
        timepoint=timepoint,
        phenotype=phenotype,
        clone=clone
    )

def open_data(filepath, timepoint, phenotype, clone_id):
    adata = sc.read(filepath)

    assert clone_id in adata.obs.columns, 'Missing field in obs: ' + clone_id
    assert timepoint in adata.obs.columns, 'Missing field in obs: ' + timepoint
    assert phenotype in adata.obs.columns, 'Missing field in obs: ' + phenotype


    return adata


def get_timepoints(adata, order, timepoint):
    if order is None:
        timepoints = adata.obs[timepoint].unique()
        assert len(timepoints) == 2, timepoint + ' column does not have 2 unique values, filter data or use order argument'

        return sorted(timepoints)
    return order


def filter_adata(adata, timepoint, phenotype, clone, timepoints, threshold):
    adata = adata[adata.obs[timepoint].isin(timepoints)]
    df = adata.obs

     ## Filter out all clones that are present in all timepoint values at threshold level

     ## !!! need to generalize to arb. timepoints

    pre = df[df[timepoint] == timepoints[0]]
    pre_counts = pre[clone].value_counts()
    pre_clones = set(pre[pre[clone].isin(pre_counts[pre_counts >= threshold].index)][clone])

    post = df[df[timepoint] == timepoints[1]]
    post_counts = post[clone].value_counts()
    post_clones = set(post[post[clone].isin(post_counts[post_counts >= threshold].index)][clone])

    clones = pre_clones.intersection(post_clones)

    return adata[df[clone].isin(clones)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Output sankey HTML with data file')
    parser.add_argument('dashboard_id', type=str)
    parser.add_argument('path', type=str)
    parser.add_argument('-t', '--threshold', type=int, help='Minimum number of cells present in a clone in each timepoint to be included', default=3)
    parser.add_argument('-o', '--order', type=str, nargs=2, help='Order of timepoints')
    parser.add_argument('--width', type=int, default=800, help='Pixel width of sankey plot')
    parser.add_argument('--height', type=int, default=700, help='Pixel height of sankey plot')
    parser.add_argument('--timepoint', type=str, default='timepoint', help='Column name for timepoint')
    parser.add_argument('--clone', type=str, default='clone_id', help='Column name for clone ID')
    parser.add_argument('--phenotype', type=str, default='cell_type', help='Column name for phenotype')

    args = parser.parse_args()

    configure_app(args.path, args.timepoint, args.phenotype, args.clone, args.order, args.threshold)

    app.run()