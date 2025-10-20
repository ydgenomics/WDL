### Date: 251014
### Image: scIB--

import scanpy as sc
import matplotlib.pyplot as plt
# import pandas
# import leidenalg
import click
import argparse

parser = argparse.ArgumentParser(description="Run unintegration analysis")
parser.add_argument("input_h5ad", type=str, default="/data/input/Files/yangdong/cotton/W202510110012842/03_integration/cotton_fibre.h5ad", 
                    help="Path to the input h5ad file")
parser.add_argument("out_h5ad", type=str, default="/data/work/cotton_fibre.h5ad",
                    help="Path to the output h5ad file")
parser.add_argument("out_umap", type=str, default="/data/work/cotton_fibre.umap.pdf",
                    help="Path to save the UMAP PDF output")
parser.add_argument("--batch_key", type=str, default="biosample",
                    help="Batch key in identifying HVG and harmony integration")
parser.add_argument("--sample_key", type=str, default="sample",
                    help="Sample key to distinguish sample")
parser.add_argument("--cluster_key", type=str, default="celltype",
                    help="Cluster key in samples")
parser.add_argument("--resolution_set", type=float, default=0.5,
                    help="Set for resolution, must be float")

args = parser.parse_args()

input_h5ad = args.input_h5ad
out_h5ad = args.out_h5ad
out_umap = args.out_umap
batch_key = args.batch_key
sample_key = args.sample_key
cluster_key = args.cluster_key
resolution_set = args.resolution_set

adata = sc.read_h5ad(input_h5ad)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key=batch_key)
# sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver="arpack")

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)
print("compute cluster using leiden, and save the account of clusters in celltype")
sc.tl.leiden(adata, resolution=resolution_set, key_added='celltype', flavor="igraph", n_iterations=2)
sc.tl.umap(adata)
sc.pl.umap(adata, color=[batch_key, sample_key, cluster_key], legend_loc='on data', ncols=1)
plt.savefig(out_umap, dpi=300, bbox_inches='tight') 
plt.close()
adata.write(filename=out_h5ad, compression="gzip")

sc.external.pp.bbknn(adata, batch_key=batch_key)
sc.tl.umap(adata)
sc.pl.umap(adata, color=[batch_key, sample_key, cluster_key, 'anno'], legend_loc='on data', ncols=1)
plt.savefig(out_umap, dpi=300, bbox_inches='tight') 
plt.close()