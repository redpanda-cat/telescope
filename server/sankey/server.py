from flask import Flask, request, jsonify
import scanpy as sc
import pandas as pd

app = Flask(__name__)

adata = sc.read('hacohen_compass.h5ad')
## !!! will need to generalize
adata.obs['viz_node'] = adata.obs['treatment'].astype(str) + "_" + adata.obs['cell_type'].astype(str)

@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"

@app.route('/cells/')
def cells():
    return []


@app.route('/degenes')
def degenes():
    timepoints = request.args.get('timepoint', []).split(',')
    phenotypes = request.args.get('phenotype', []).split(',')
    clones = request.args.get('clone', None)

    nodes = [f"{a[0]}_{a[1]}" for a in zip(timepoints, phenotypes)]

    if clones is None:
        adata.obs['included'] = adata.obs.apply(lambda x: "included" if x['viz_node'] in nodes else "excluded", axis=1)
    else:
        clones = clones.aplit(',')
        adata.obs['included'] = adata.obs.apply(lambda x: "included" if x['viz_node'] in nodes and x['clone_id'] in clones else "excluded", axis=1)

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

    return response