# ARBOLpy
python implementation of the R package ARBOL, scRNAseq iterative tiered clustering

![](https://github.com/jo-m-lab/ARBOLpy/blob/main/docs/ARBOLsmall.jpg?raw=true)

Iteratively cluster single cell datasets using a scanpy anndata object as input. Identifies and uses optimum 
clustering parameters at each tier of clustering. Current build includes SCtransform normalization. 
Outputs QC and visualization plots for each clustering event.  

## Install

By github:
```
pip install git+https://github.com/jo-m-lab/ARBOLpy.git
```

from PyPI
```
pip install arbolpy

import ARBOL
```

or clone the repository and source the functions directly from the script
```
git clone https://github.com/jo-m-lab/ARBOLpy.git

import "path/to/cloned/git/repo/ARBOLpy/ARBOL"
```

there is a docker image available with ARBOL and dependencies preinstalled
https://hub.docker.com/r/kkimler/arbolpy

## Recommended Usage

ARBOL was developed and used in the paper, "A treatment-naïve cellular atlas of pediatric Crohn’s disease predicts disease severity and therapeutic response"
Currently, a tutorial is only available for the R version, where the FGID atlas figure is reproduced: 
https://jo-m-lab.github.io/ARBOL/ARBOLtutorial.html

ARBOLpy is a stripped down version of ARBOL meant to perform iterative clustering with little overhead. 
Currently it does not include the two stop conditions that the R version uses to heuristically join similar clusters.
This results in the Python version overclustering data. Methods for merging the end clusters of the tree are available on the develop branch of the R version of ARBOL.

This package is meant as a starting point for the way that we approached clustering and and is meant to be edited/customized through community feedback through users such as yourself!

The main function of ARBOLpy is ARBOL() - here is an example call. 

```
import scanpy as sc
import ARBOL

adata = sc.datasets.pbmc3k()

tree = ARBOL.ARBOL(adata)

ARBOL.write_ARBOL_output(tree,output_csv='endclusts.csv')
```

The helper function write_ARBOL_output writes the anytree object's endclusters to a csv file.

**Note** This script can take a long time to run. Running on 20K cells could 
take >30 minutes. Running on 100k+ cells could take >3 hours. 

**Note** It has been tested up to 200k cells, and beyond 10k cells, maintains a linear relationship between resource usage and number of cells

**Python ARBOL resource usage**:  
	Pearson residuals normalization:  
 	- 1.2 GB RAM per 1000 cells  
 	- 2 minutes per 1000 cells  
 	TPM normalization:  
 	- 1.2 GB RAM per 1000 cells  
 	- 1:55 min per 1000 cells

**R ARBOL resource usage**:  
	Pearson residuals normalization (SCTransform):  
	- 1.2 GB RAM per 1000 cells  
	- 4 minutes per 1000 cells  

 The current RAM/time bottleneck is the silhouette analysis, which runs 30 rounds of clustering at different resolutions. 

## ARBOL() Parameters

* *adata* scanpy anndata object
* *normalize_method* normalization method, defaults to "Pearson", scanpy's experimental implementation of SCTransform. Also available: "TPM": as implemented in scanpy normalize_total()
* *tier* starting tier, defaults to 0
* *cluster* starting cluster, defaults to 0
* *min_cluster_size* minimum number of cells to allow further clustering
* *tree* anytree object to attach arbol to. Shouldn't be changed unless building onto a pre-existing tree.
* *parent* parent node of current clustering event, defaults to None. As with tree, shouldn't be changed unless building onto a pre-existing anytree object
* *max_tiers* maximum number of tiers to allow further clustering
* *min_silhouette_res* lower bound of silhouette analysis leiden clustering resolution parameter scan 
* *max_silhouette_res* upper bound
* *silhouette_subsampling_n* number of cells to subsample anndata for silhouette analysis (cluster resolution choice)
* *h5dir* where to save h5 objects for each tier and cluster, if None, does not save
* *figdir* where to save QC and viz figures for each tier and cluster, if None does not save

## Returns

* anytree object based on iterative tiered clustering
