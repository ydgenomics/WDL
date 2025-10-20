# Date: 251014
# Image: /software/conda/Anaconda/bin/python 
# Coder: ydgenomics(yangdong@genomics.cn)

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import webbrowser
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import click

@click.command()
@click.option("--unintegrated_h5ad", type=str, default="/data/work/scIB/unintegrated.h5ad", help="Path to unintegrated h5ad file")
@click.option("--integrated_h5ad", type=str, default="/data/work/scIB/Scanorama.h5ad", help="Path to integrated h5ad file")
@click.option("--methods_list", type=str, default="Scanorama", help="List of methods separated by comma")
@click.option("--pcas_list", type=str, default="Scanorama", help="List of PCAs separated by comma")
@click.option('--batch_key', type=str, default="batch", help="Batch key")
@click.option('--label_key', type=str, default="cell_type", help="Information of biological cell name")
@click.option('--n_jobs', type=int, default=4, help="Number of jobs to use for parallelization of neighbor search")
@click.option("--prefix", type=click.Path(exists=False), default="lung", help="Prefix for output files")
def main(unintegrated_h5ad, integrated_h5ad, methods_list, pcas_list, batch_key, label_key, n_jobs, prefix):
    # Read files and split by comma
    files = integrated_h5ad.split(',')
    methods = methods_list.split(',')
    pcas = pcas_list.split(',')
    
    out_benchpdf=prefix+"_scIB.pdf"; print(out_benchpdf)
    out_benchcsv=prefix+"_scIB.csv"; print(out_benchcsv)
    out_h5ad=prefix+"_scIB.h5ad"; print(out_h5ad)

    # Process unintegrated data
    orig_ad = sc.read_h5ad(unintegrated_h5ad)
    # sc.set_figure_params(dpi_save=200, frameon=False, figsize=(10, 5))
    # sc.pp.normalize_total(orig_ad, target_sum=1e4)
    # sc.pp.log1p(orig_ad)
    # sc.pp.highly_variable_genes(orig_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # sc.pp.scale(orig_ad, max_value=10)
    # sc.tl.pca(orig_ad, svd_solver="arpack")
    # sc.pp.neighbors(orig_ad, n_neighbors=20, n_pcs=20)
    # sc.tl.leiden(orig_ad, resolution=0.5, key_added=cluster_key, flavor='igraph', n_iterations=2, directed=False)
    # sc.tl.umap(orig_ad, min_dist=0.3)
    # sc.pl.umap(orig_ad, color=[batch_key, label_key, cluster_key])
    # plt.savefig(out_rawpdf, dpi=300, bbox_inches='tight')
    orig_ad.obsm["Unintegrated"] = orig_ad.obsm["X_pca"]

    # Process integrated data and merge with unintegrated data
    for i in range(len(files)):
        adata = sc.read_h5ad(files[i])
        adata.obsm[methods[i]] = adata.obsm[pcas[i]]
        orig_ad.obsm[methods[i]] = adata.obsm[methods[i]]

    methods.append('Unintegrated')

    import time
    start = time.time()
    
    bm = Benchmarker(
        adata=orig_ad,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=methods,
        pre_integrated_embedding_obsm_key="X_pca",
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
        n_jobs=n_jobs,
    )
    bm.benchmark()
    end = time.time()
    orig_ad.write(out_h5ad, compression="gzip")
    print(f"Time: {int((end - start) / 60)} min {int((end - start) % 60)} sec")

    bm.plot_results_table()
    bm.plot_results_table(min_max_scale=False)
    plt.savefig(out_benchpdf, format='pdf', bbox_inches='tight')
    plt.close()
    df = bm.get_results(min_max_scale=False)
    print(df)
    df_transposed = df.transpose()
    df_transposed.to_csv(out_benchcsv, index=True)

if __name__ == "__main__":
    main()
