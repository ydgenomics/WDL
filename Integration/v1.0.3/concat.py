### Date: 250723
### Image: harmony-py--
### Coder: ydgenomics

import pandas as pd
import scanpy as sc
import anndata as ad
import sys

files_txt_path = sys.argv[1]
projects_txt_path = sys.argv[2]
species=sys.argv[3]
group_key=sys.argv[4] #"biosample"

with open(files_txt_path, 'r') as file:
    file_content = file.read().strip()

indataget = file_content.split(',')
print(indataget)

with open(projects_txt_path, 'r') as filen:
    file_content = filen.read().strip()

projects = file_content.split(',')
print(projects)

adatas={}
for i in range(len(indataget)): 
    key = projects[i]
    value = indataget[i]
    value = sc.read_h5ad(value)
    if 'biosample' in value.obs.columns:
        print('The raw key included `biosample` column, value of raw biosample named to biosample0')
        value.obs['biosample0']=value.obs['biosample']
    if 'counts' in value.layers:
        print("counts exist in raw h5ad")
    else:
        print("counts not exist in raw h5ad, .X as counts")
        value.layers["counts"] = value.X.copy()
    value.X = value.layers["counts"] #ensure concat used by raw data
    adatas[key] = value

adata = ad.concat(adatas, label=group_key, join="inner") # 'inner' or 'outer'
if 'celltype' in adata.obs.columns:
    print('The raw key included `celltype` column, value of raw celltype named to celltype0')
    adata.obs['celltype0']=adata.obs['celltype']

adata.obs_names_make_unique()
print(adata.obs[group_key].value_counts())
print(adata.obs.columns)
adata.write_h5ad(filename=species+'.h5ad',compression="gzip")