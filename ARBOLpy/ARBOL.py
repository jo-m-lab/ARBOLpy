import numpy as np
import pandas as pd
import scanpy as sc
import os
from scipy.sparse import issparse
import seaborn as sb
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from anytree import Node, RenderTree, find_by_attr
from anytree.exporter import DictExporter
from collections import OrderedDict
import anytree
import math

sc.settings.verbosity = 4

#Functions used in wrapper

def QC_Plotting(adata):
    print('plotting QC')
    p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='pct_counts_mt',save='_QC_counts_genes.png')
    p2 = sc.pl.scatter(adata[adata.obs['n_counts']<5000], 'n_counts', 'n_genes', color='pct_counts_mt',save='_QC_counts_genes_lowEnd_zoom.png')
    t1 = sc.pl.violin(adata, keys=['n_counts','n_genes','pct_counts_mt'], size=2, log=True, cut=0,save='_QC_counts_genes_mt.png')

def TPM_log1p_Normalize(adata):
    print('starting normalization')
    print('copy X')
    adata.raw = adata
    print('normalizing')
    sc.pp.normalize_total(adata,target_sum=10**4)
    print('log transforming')
    sc.pp.log1p(adata)
    # returning adata feels unnecessary since everything scanpy does is inplace anyway
    # Interestingly though it does not print out the standard adata info you normally get when running 'adata'
    return adata

def Find_X_Variable_Genes(adata,num_genes):
    sc.pp.highly_variable_genes(adata,n_top_genes=num_genes)
    print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
    sc.pl.highly_variable_genes(adata)
    return adata

def Find_X_Variable_Genes_Pearson(adata,num_genes,figdir):
    sc.experimental.pp.highly_variable_genes(adata, flavor='pearson_residuals', n_top_genes=num_genes)
    print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
    #no residual variance plot currently. custom plot is difficult within experimental scanpy module
    return adata

def Pearson_Residuals_Normalize(adata):
    adata.raw = adata
    sc.experimental.pp.normalize_pearson_residuals(adata)
    return adata


def PC_choice_plot(adata,axis,diff):
    print('xcb')
    xsb = range(1,len(adata.uns['pca']['variance_ratio'])+1)
    print('plotting PCAs')
    ax0 = sb.scatterplot(y=adata.uns['pca']['variance_ratio'],x=np.array(xsb),ax=axis)
    ax0.set(title='variance explained per PC, cutoff when loss < 15%')
    ax0.axvline(max(max(np.where(diff > np.quantile(diff,0.85)))),linewidth=2, color='r')
    return ax0

def decide_PCs(adata):
    print('starting PCA calc')
    sc.pp.pca(adata,n_comps=50,use_highly_variable=True,svd_solver='arpack')
    t = adata.uns['pca']['variance_ratio']
    print('calculating PCA diff')
    diff = t[0:len(t)-1] - t[1:len(t)]
    print('calculating good PCs')
    goodPCs = np.where(diff > np.quantile(diff,0.85))[0]
    
    fig, ax = plt.subplots(1,2,figsize=(12,6))
    
    ax0 = PC_choice_plot(adata,axis=ax[0],diff=diff)
    ax1 = sc.pl.pca_scatter(adata, color='n_counts',ax=ax[1])
    
    return goodPCs,fig
    

def clustvres_plot(scores,axis,res_array,nclusts):
    ax3 = sb.scatterplot(y=nclusts,x=np.round(res_array, 1),ax=axis)
    ax3.set(title='Resolution vs. Number of Clusters')
    ax3.axvline(res_array[np.argmax(scores)],linewidth=2,color='r')
    return ax3

def silvres_plot(scores,axis,res_array):
    ax1 = sb.scatterplot(y=scores,x=np.round(res_array, 1),ax=axis)
    ax1.set(title='Resolution vs. Average Silhouette')
    ax1.axvline(res_array[np.argmax(scores)],linewidth=2, color='r')
    return ax1

def silvclust_plot(scores,axis,nclusts):
    ax2 = sb.scatterplot(y=scores,x=nclusts,ax=axis)
    ax2.set(title='Number of Clusters vs. Average Silhouette')
    ax2.axvline(nclusts[np.argmax(scores)],linewidth=2,color='r')
    return ax2

def silhouette_param_scan(adata,res_array):
    tmpadata = adata.copy()
    scores=[]
    clusts = []
    used_res = []
    print('testing out different resolutions of leiden clustering')
    for res in res_array:
        print(f'running leiden clustering at {res} resolution')
        sc.tl.leiden(tmpadata, resolution=res)
        if len(tmpadata.obs[f'leiden'].unique()) == 1:
            # silhoutte analysis does not work if it only finds one cluster, 
            # and such a resolution is useless anyway
            continue
        clusts.append(len(np.unique(tmpadata.obs[f'leiden'])))
        scores.append(silhouette_score(tmpadata.obsm['X_pca'],
                                  tmpadata.obs[f'leiden']
                                 ))
        used_res.append(res)
    
    fig, ax = plt.subplots(1,3,figsize=(16,6))
    if len(used_res) == 0:
        # aka all resolutions led to 1 cluster being produced
        return 0.01,fig
    ax1 = silvres_plot(scores,axis=ax[0],res_array=used_res)
    ax2 = silvclust_plot(scores,axis=ax[1],nclusts=clusts)
    ax3 = clustvres_plot(scores,axis=ax[2],res_array=used_res,nclusts=clusts)
    
    bestres = used_res[np.argmax(scores)]
    
    if max(scores) < 0:
        bestres = 0.01
    
    return bestres,fig


def PCA_Neighbors(adata,n_pcs):
    # didnt we calculate the PCs already? to set the number of PCs you should just be able to set n_pcs in the neighbors setting
    sc.pp.pca(adata,n_comps=n_pcs,use_highly_variable=True,svd_solver='arpack')

    print('calculate neighbors')
    # knn gets quite big for this number of cells (~400 for 700k cells), this might slow it down unneccesarrily
    # the lung cell atlas pre-print (2.2M cells) uses 30 neighbors for the full scale, then they do sub-clustering and use 15 neighbors https://www.biorxiv.org/content/10.1101/2022.03.10.483747v1.full
    # so, limit at 30 neighbors?
    # sc.pp.neighbors(adata,n_pcs=n_pcs,knn=math.ceil(0.5*math.sqrt(adata.X.shape[0])))
    sqrt_size = math.ceil(0.5*math.sqrt(adata.X.shape[0]))
    print(f'sqrt nr of cells (knn size) {sqrt_size}')
    if sqrt_size > 30:
        sqrt_size = int(30)
        print('knn > 30, setting to 30')
    sc.pp.neighbors(adata,n_pcs=n_pcs,n_neighbors=sqrt_size)
    print('finished neighbor calculation')
    return adata

def Cluster_UMAP(adata,bestres,node):
    print('calculating leiden')
    sc.tl.leiden(adata,resolution=bestres,key_added=node.rstrip('.'))
    print('calculating UMAP')
    sc.tl.umap(adata)
    # first umap here requires an obs column 'sample'
    sc.pl.umap(adata,color='sample',save='_sample.pdf')
    sc.pl.umap(adata,color=node.rstrip('.'),save='_cluster.pdf')
    return adata

def SplitAdataOnIdents(adata,ident):
    clusters = np.unique(adata.obs[ident])
    adatas=[adata[adata.obs[ident]==cluster] for cluster in clusters]
    return adatas

def Sort_cluster_list(cluster_list):
    cluster_list = sorted([int(i) for i in cluster_list])
    return [str(i) for i in cluster_list]

def Cluster_DE(adata,figdir,node,max_cells_per_cluster=2000):
    
    stripped_node = node.rstrip(".")
    cluster_sizes = adata.obs.groupby(stripped_node).size()
    
    if cluster_sizes.max() > max_cells_per_cluster:
        #adata.obs[f'subsampled_{stripped_node}'] = 'other'
        cluster_list = list(adata.obs[stripped_node].unique())
        cluster_list = Sort_cluster_list(cluster_list)
        print(cluster_list)
        obs_name_list = []
        
        for cluster in cluster_list:
            if cluster_sizes.loc[cluster] > max_cells_per_cluster:
                cluster_obs_names = adata[adata.obs[stripped_node]==cluster,:].obs_names
                sub_obs_names = np.random.choice(cluster_obs_names, size=max_cells_per_cluster, replace=False)
                obs_name_list+=list(sub_obs_names)
            else:
                cluster_obs_names = adata[adata.obs[stripped_node]==cluster,:].obs_names
                obs_name_list+=list(cluster_obs_names)
        DEGadata = adata[obs_name_list,:].copy()
        
        # run DE only on the subsampled 2000 cells per cluster through the groups= parameter 
        sc.tl.rank_genes_groups(DEGadata, stripped_node, method='wilcoxon')
        sc.pl.rank_genes_groups_heatmap(DEGadata,save=f'_wilcoxon_top10.pdf')
        sc.get.rank_genes_groups_df(DEGadata, group=None).to_csv(f'{figdir}/wilcoxon_1vAll_perCluster.csv')
        
        return DEGadata.uns['rank_genes_groups']
    else:
        sc.tl.rank_genes_groups(adata, stripped_node, method='wilcoxon')
        sc.pl.rank_genes_groups_heatmap(adata,save=f'_wilcoxon_top10.pdf')
        sc.get.rank_genes_groups_df(adata, group=None).to_csv(f'{figdir}/wilcoxon_1vAll_perCluster.csv')


### Wrapper function

def ARBOL(adata, tier = 0, cluster = 0,
          tree = '',
          parent = None,
          normalize_method='Pearson',
          max_tiers=10, min_cluster_size=100,
          silhoutte_subsampling_n=7500,
          min_silhouette_res = 0.004,
          max_silhouette_res = 3,
          figdir = 'figs', h5dir = 'h5s'):
    
    #create save folders in case they're not there
    if os.path.isdir(h5dir) == False:
        os.mkdir(h5dir)
    if os.path.isdir(figdir) == False:
        os.mkdir(figdir)
    
    #make sure QC metadata exists
    if tier == 0:
        adata.obs['n_counts'] = adata.X.sum(1)
        adata.obs['log_counts'] = np.log(adata.obs['n_counts'])[0]
        adata.obs['n_genes'] = (adata.X > 0).sum(1)
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    #keep track of position in tree for logging
    print("Starting tier: ",tier)
    print("cluster: ",cluster)
    print("Number of cells: ",adata.X.shape[0])
    if tier == 0:
        print("Number of genes: ",adata.X.shape[1])
    
    print(figdir)
    #QC plotting (done only once)
    if tier == 0 & cluster == 0:
        sc.settings.figdir = figdir
        QC_Plotting(adata)
    
    
    #Set node 
    if tier == 0:
        node = f'T{tier}C{cluster}.'
        cluster += 1
    else:
        node = f'{parent}T{tier}C{cluster}.'
    if parent==None:
        node=node.lstrip('None')
    
    #Write node to tree
    if node == 'T0C0.':
        tree = Node(node.rstrip('.'),
                    n=adata.X.shape[0],
                    tier=tier,
                    cells=adata.obs_names.tolist())
    else:
        Node(node.rstrip('.'),
             parent=find_by_attr(tree, parent.rstrip('.')),
             n=adata.X.shape[0],
             tier=tier,
             cells=adata.obs_names.tolist())
    
    #Quit if too few cells
    if adata.X.shape[0] < min_cluster_size:
        return tree
    
    # set figure output directory to hierarchical structure
    # create save folders in case they're not there
    
    #figdir = f"{figdir}/{node.split('.')[-2]}./"    # take only last element of node as figdir should contain hierachically clustered folders
    homedir = figdir
    figdir = f"{figdir}/{node.rstrip('.')}/" 
    if os.path.isdir(figdir) == False:
        os.mkdir(figdir)
    sc.settings.figdir = figdir

  
    #pre-processing
    
    if normalize_method == 'Pearson':
        #filter genes with zero counts to avoid errors in sctransform
        sc.pp.filter_genes(adata,min_counts=3)
        Find_X_Variable_Genes_Pearson(adata,2000,figdir=figdir) 
        Pearson_Residuals_Normalize(adata)

    if normalize_method == 'TPM':
        # log transform only once
        if tier == 0:
            TPM_log1p_Normalize(adata)  
        Find_X_Variable_Genes(adata,2000)    
    
    #find optimum PCs to use for clustering
    a,pcaplt = decide_PCs(adata)
    pcaplt.savefig(f'{figdir}/gene_choice_dispersion.png')
    
    min_res_power = np.log10(min_silhouette_res)
    max_res_power = np.log10(max_silhouette_res)
    
    if (adata.shape[0] > silhoutte_subsampling_n) & (silhoutte_subsampling_n < math.inf):
        print(f'''number of cells ({adata.shape[0]}) > {silhoutte_subsampling_n}, this is the limit used for the silhouette analysis
        subsampling to {silhoutte_subsampling_n} cells''')
        
        tmpadata = sc.pp.subsample(adata,n_obs=silhoutte_subsampling_n,random_state=42,copy=True)
        print(f'shape of subsampled adata: {tmpadata.shape}')
        
        print(f'calculating neighbors with {max(a)} PCs')
        PCA_Neighbors(tmpadata,n_pcs=max(a))

        #set preferred figure size for silhouette analysis
        figure(figsize=(12, 6), dpi=80)

        #find optimum cluster resolution by silhouette analysis
        print('starting silhouette analysis')
        bestres,fig = silhouette_param_scan(tmpadata,res_array = list(np.logspace(min_res_power, max_res_power, 30)))
        
        print(f'best resolution by silhouette analysis: {bestres}')

        fig.savefig(f"{figdir}/Silhouette_Analysis.png",dpi=199)
        
        # we can del tmpadata since we only needed bestres
        del tmpadata
        
        #now run at full scale with known best res
        print(f'calculating neighbors with {max(a)} PCs')
        PCA_Neighbors(adata,n_pcs=max(a))
        #Cluster
        print('clustering to UMAP space')
        Cluster_UMAP(adata,bestres=bestres,node=node)
    else:
        print(f'calculating neighbors with {max(a)} PCs')
        PCA_Neighbors(adata,n_pcs=max(a))

        #set preferred figure size for silhouette analysis
        figure(figsize=(12, 6), dpi=80)

        #find optimum cluster resolution by silhouette analysis
        print('starting silhouette analysis')
        bestres,fig = silhouette_param_scan(adata,res_array = list(np.logspace(min_res_power, max_res_power, 30)))

        print(f'best resolution by silhouette analysis: {bestres}')

        fig.savefig(f"{figdir}/Silhouette_Analysis.png",dpi=199)   

        #Cluster
        print('clustering to UMAP space')
        Cluster_UMAP(adata,bestres=bestres,node=node)
        
    #Quit if only one cluster found
    if len(np.unique(adata.obs[node.rstrip('.')])) == 1:
        return tree    
    
    #Write node to h5
    # f-strings are prettier i think: adata.write(f"{h5dir}/adata_{node.rstrip('.')}.h5ad")
    # also, people cant be trusted to add a / at the end of their directories, so make sure to add that
    # add compression to writing?
    adata.write_h5ad(f"{h5dir}/adata_{node.rstrip('.')}.h5ad",compression='gzip')
    
    #Plot and write Wilcoxon DE per cluster
    adata.uns['rank_genes_groups'] = Cluster_DE(adata,figdir=figdir,node=node,max_cells_per_cluster=2000)
    
    #Split anndata into subsets per cluster

    #reset anndata.X to raw counts if using pearson residual normalization
    #because need to re-calc for each subset. For TPM, no need to renorm.
    if normalize_method == 'Pearson':
        print('resetting anndata to raw anndata')
        adata = adata.raw.to_adata()

    subadatas = SplitAdataOnIdents(adata,ident=node.rstrip('.'))  

    #if we continue clustering, set cluster equal to 0, preserve parent and node for reset when coming back up tree
    cluster_hold = cluster
    node_hold = node
    cluster = 0
    newtier = tier + 1
    
    if newtier < max_tiers:

        for subadata in subadatas:
            #make sure h5dir and figdir are specified here as well
            tree = ARBOL(subadata,tier=newtier,
                         cluster=cluster,
                         tree=tree,parent=node,
                         normalize_method=normalize_method,
                         silhoutte_subsampling_n=silhoutte_subsampling_n,
                         min_silhouette_res=min_silhouette_res,
                         max_silhouette_res=max_silhouette_res,
                         h5dir=h5dir,
                         figdir=homedir
                        )
            cluster = cluster + 1
        #reset node and cluster
        node = node_hold
        cluster = cluster_hold
    
    return tree


def write_ARBOL_output(tree):
    for pre, fill, node in anytree.RenderTree(tree):
        print("%s%s n=%s tier=%s" % (pre, node.name, node.n, node.tier))

    endnodes = anytree.findall(tree, filter_=lambda node: len(node.children)==0)  
    cells = []
    tier = []
    endcluster = []
    for item in endnodes:
        cells += item.cells
        tier += [item.tier] * len(item.cells)
        endcluster += [item.name] * len(item.cells)
    endf = pd.DataFrame(index=cells)
    endf['tier']=tier
    endf['endcluster']=endcluster

    endf.to_csv(figdir + 'endclusters.csv')

