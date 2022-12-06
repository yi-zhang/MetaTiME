"""
Automatic annotator for tumor single cell data. 
# 
# """
### MetaTiME annotator
import pandas as pd
import anndata
from metatime import config
import scanpy as sc
from typing import Optional, Union
from metatime import mecs
import warnings



 
########### top1 cell state annotator #######
def annotator( dat: Union[pd.DataFrame, anndata.AnnData],
              mecnamedict: dict, 
              gcol ='overcluster', 
              #MINCELL = 5 
              ):  
    """ 
    Annotator marking top1st enriched cell states based on MeC scores.

    Parameters
    ----------
    dat
        dataframe with gene by cell scores and one column gcol (grouping column) e.g. leiden class
        or,
        anndata with gene by cell scores, and one column gcol in dat.obs 
        Malignant cells shall be removed to keep tumor microenvironmental cells only, such as immune cells, fibroblasts, endothelial cells
    mecnamedict
        for renaming dat columns into functional names
        can be loaded from pre-computed tumor MeC functional annotation
    gcol
        for grouping cells. e.g. a column for overclustered cluster assignment.

    Returns
    ----------
    projmat: pd.DataFrame
         Dataframe with mec projection scores and newly added predicted label column 'MetaTiME_'+gcol
    gpred
        median score of each gcol group
    gpreddict
        dictionary how is each gcol mapped to a label
    
    Examples
    ----------
    >>> projmat, gpred, gpreddict = annotator( projmat ,  mecnamedict)
    
    """
    if (isinstance( dat , anndata.AnnData)):
        # subselecting columns
        pdata = dat.copy()
        Xproj = pdata.to_df()
        #projmat = Xproj[Xproj.columns.intersection(config.scorectcols)].join(pdata.obs[[ gcol ]])
        projmat = Xproj.join(pdata.obs[[ gcol ]])
    else:
        projmat = dat.copy()

    projmat = projmat[ list( mecnamedict.keys() ) + [ gcol ]]

    # mean or median score for each group
    gpred = projmat.groupby( gcol ).median() # gpred: leiden cluster by mec.
    # top 1 score for each group, with filtering
    gpreddict = projmat.groupby( gcol ).median().idxmax(axis=1).apply(lambda t: mecnamedict[t])# cluster max quantile 75% larger than zscore-2 
    pass_maxmedian = (projmat.groupby( gcol ).quantile(q=0.7) >1 ).max(axis=1) # given median is max, top 30%cells shall pass zscore 1. Otherwise this may be a noisy group, setting to 'Others'
    gpreddict[ ~ pass_maxmedian] = 'Others'
    # generate top1 prediction
    projmat['MetaTiME_'+gcol] = projmat[ gcol ].apply(lambda t: gpreddict[t])
    
    # add to anndata
    if (isinstance( dat , anndata.AnnData)):
        if('MetaTiME_'+gcol in pdata.obs.columns):
            del pdata.obs['MetaTiME_'+gcol]
        pdata.obs['MetaTiME_'+gcol] = projmat['MetaTiME_'+gcol]
        pdata.obs['MetaTiME_'+gcol]=pd.Categorical(pdata.obs['MetaTiME_'+gcol])
        dat = pdata
    else:
        dat = projmat
    
    return( dat, gpred, gpreddict )


def overcluster(adata: anndata.AnnData,
                resolution : float=8, 
                random_state: int= 0, 
                clustercol :str = 'overcluster'):

    """
    Overcluster single cell data to get cluster level cell state annotation

    Parameters
    ----------
    adata
        scanpy object with adata.uns['neighbors'] computed.
        if adata.obsm['X_umap'] does not exist, recomputes umap coordinates.
        otherwise, keep the umap coordinates
    resolution
        clustering resolution
    random_state
        clustering random state
    clustercol:
        key to add to adata.obs that records cluster assignment

    Returns
    ----------
        scanpy object with clustering results.

    """
    sc.tl.leiden(adata, resolution=resolution, key_added = clustercol, random_state=random_state)
    if 'X_umap' not in adata.obsm.keys():
        sc.tl.umap(adata)
    return(adata)

def pdataToTable( pdata: anndata.AnnData, 
                    mectable: pd.DataFrame, 
                    gcol : str= 'overcluster'):

    """
    Convert projected scores to two simple pandas dataframes

    Parameters
    ----------
    pdata
        anndata with gene by mec scores, and one column gcol in pdata.obs 
    mectable
        for renaming dat columns into functional names
        can be loaded from pre-computed tumor MeC functional annotation
        Required columns: `['MeC_id', 'Annotation', 'UseForCellTypeAnno']`  
    gcol
        a column in pdata.obs for grouping cells. a column for overclustered cluster assignment.

    Returns
    ----------
    tuple
        projmat : a pandas dataframe with mec scores and a column for grouping cells. useful for annotating cell states. 
        mecscores: a pandas dataframe for per-cell mec scores. columns use functional annotation of mec ids


    """
    # projmat
    projmat = pdata.obs.join( pdata.to_df() )
    allmecids = mectable['MeC_id'].values.tolist()
    projmat = projmat[ allmecids  + [ gcol ] ]

    # mecscores
    mecscores = pdata.to_df()
    mecctnamedict = mecs.getmecnamedict_ct( mectable )
    newcols = []
    for t in projmat.columns:
        if(t in mecctnamedict.keys() ):
            newcols.append( mecctnamedict[t]  )
        else:
            newcols.append(t)
    mecscores = projmat.copy()
    mecscores.columns = newcols
    
    return(projmat, mecscores)


def saveToAdata( adata : anndata.AnnData, 
                projmat : pd.DataFrame, 
                gcol: str = 'overcluster',
                ANNOTATION_ONLY : bool= False
                ):
    """
    Save annotation to adata. 

    Parameters
    ----------
    adata: anndata.AnnData
        scRNA scanpy object
    projmat: pd.DataFrame
        Dataframe with mec projection scores and newly added predicted label column 'MetaTiME_'+gcol
    gcol
        A column in projmat for grouping cells. Typically a column for overclustered cluster assignment.
    ANNOTATION_ONLY
        Whether to only add cluster-wise annotation to anndata or Also append scores to adata.obs

    Returns
    ----------
    anndata.AnnData
        adata with per-cluster annotation column in adata.obs[[gcol]], and scores appended to adata.obs if ANNOTATION_ONLY==True

    """
    exist_cols = projmat.columns.intersection(adata.obs.columns)
    if( len(exist_cols ) >0 ):
            warnings.warn('Columns in projmat already overlap columns in adata.obs. Overwriting adata.obs. '
                )
            adata.obs = adata.obs[[ t for t in adata.obs.columns if t not in exist_cols] ]
    
    if( ANNOTATION_ONLY ):
        adata.obs['MetaTiME_'+gcol] = projmat['MetaTiME_'+gcol]
    else:
        adata.obs = adata.obs.merge( projmat , how = 'left', left_index = True, right_index = True)
    return( adata )
        

def saveToPdata( pdata : anndata.AnnData, 
                adata: anndata.AnnData,
                projmat : pd.DataFrame, 
                gcol: str = 'overcluster',
                BORROW_ADATA_EMBEDDING=True,
                ):
    """
    Save annotation to pdata.
    Borrow embedding from adata for easy visualization, including adata.obsm['X_pca'], adata.obsm['X_umap'], adata.obsm['X_pca_harmony']

    Parameters
    ----------
    pdata : anndata.AnnData
        scanpy object for per-cell projected score. 
    adata: anndata.AnnData
        scRNA scanpy object
    projmat: pd.DataFrame
        Dataframe with mec projection scores and newly added predicted label column 'MetaTiME_'+gcol
    gcol
        A column in projmat for grouping cells. Typically a column for overclustered cluster assignment.
    BORROW_ADATA_EMBEDDING
        Whether to borrow pca and umap embeddings from adata to write in pdata. For easy visualization.

    Returns
    ----------
    anndata.AnnData
        pdata with per-cluster annotation column in pdata.obs[[gcol]]. 

    """
    if(BORROW_ADATA_EMBEDDING):
        pdata.obsm['X_pca'] = adata.obsm['X_pca']
        pdata.obsm['X_umap'] = adata.obsm['X_umap']
        if('X_pca_harmony' in adata.obsm):
            pdata.obsm['X_pca_harmony'] = adata.obsm['X_pca_harmony']

    pdata.obs[gcol] = projmat[[gcol]]
    pdata.obs['MetaTiME_'+gcol] = projmat['MetaTiME_'+gcol]
    return( pdata )


