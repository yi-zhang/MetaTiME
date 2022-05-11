import pandas as pd
import os
import scanpy as sc
import numpy as np
import warnings
import anndata



############ project single cell to mec for creating new adata obj ###########

def projectMec(expr, mec, glstcol = 'TopGene20'):

    """ 
    Parameters
    ----------
    expr
        expression matrix, cell by gene. 
    mec
        mec table with one column, two types of format accepted
        format 1, each row is a mec with genes , comma seperated.
        format 2, each row is a genes, each column is a mec, value is float
    glstcol
        Used only when mec is the list format format1, and glstcol is the column name to record comma separated gene list for each mec.

    Returns
    ----------
    scorepd, cell by signature.

    Examples
    ---------
    >>> scorepd = mecmapper.projectMec( df, mec )
    
    """
    score_lst = []
   
    if(isinstance( mec.iloc[0,0], str) ) :
        # the mec format is format 1. no weighting.
        for i in mec.index:
            mec_i = str(i) 
            gs = mec.loc[i, glstcol ].split(',')
            gs_sub = [ t for t in expr.columns if t in gs ]
            score = pd.DataFrame( expr[ gs_sub ].mean(axis=1), columns= [ mec_i ] )
            score_lst.append( score  )
        scorepd = pd.concat(score_lst , axis=1).fillna(0)
    else:
        # the mec format is not format 1. using values as weight.
        shared_genes = mec.index.intersection( expr.columns )
        module_sharedgenes = mec.loc[ shared_genes ]
        expr_sharedgenes = expr[ shared_genes ]
        scorearray = expr_sharedgenes.values.dot( module_sharedgenes.values )
        scorepd = pd.DataFrame( scorearray, index = expr.index, columns = mec.columns)
        
    return( scorepd )



def annToDataFrame( adata_input, genescaling = False,  layer = 'norm_data'):
    """ 
    Extract expression matrix from scanpy object.
    Parameters
    ----------
    adata_input: anndata.AnnData
        Input scanpy object. 
    genescaling: bool
        Whether to z-scale extracted feature matrix.
    layer: str
        Layer of expression to extract from adata_input

    Returns
    ----------
    pd.dataframe
        Extracted expression matrix from scanpy object.
        
    """
    adata = adata_input.copy()
    try:
        df = adata.to_df(layer = layer )
    except:
        warnings.warn('Warning: no '+layer +' layer. using X.')
        df = adata.to_df()
    if(  genescaling ):
        df = (df - df.mean())/df.std()
    
    return( df )


def scale(df):
    """ standardize scaling feature """
    return(  (df-df.mean())/df.std()  )

def projectMecAnn( adata_input, 
                mec, 
                genescaling = False,  
                sigscaling = True, 
                addon = False, 
                layer = 'norm_data', 
                glstcol = 'TopGene20'):
    """
    Project single cell expression in AnnData to MeCs
    Calls: projectMeC,annToDataFrame

    Parameters
    ----------
    adata_input: anndata.AnnData
        Input scanpy object for gene expression
    mec: pd.DataFrame
        mec table with one column, two types of format both accepted
        format 1, each row is a mec with genes , comma seperated.
        format 2, each row is a genes, each column is a mec, value is float.
    genescaling: bool
        Whether to scale expression on gene level. Recommended to be False.
    sigscaling: bool
        Whether to scale projected scores across cells. Recomended and default is True.
    addon: bool
        Whether to keep the original adata and append the signautres in obs, or return an independent anndata (pdata) with only projected values (which saves memory).
    layer: str
        Layer of expression to extract from adata_input
    glstcol
        Used only when mec is the list format format 1, and glstcol is the column name to record comma separated gene list for each mec.


    Returns
    ----------
    anndata.AnnData
        If addon is False, return pdata where values are MeC-projected values. 
        If addon is True, return adata where values are same as in adata_input, but with extra obs columns.
    Examples
    ----------
    >>> pdata = mecmapper.projectMecAnn(adata, mec_score_topg, sigscaling=True, genescaling=False, addon=False)
    """

    adata = adata_input.copy()
    df = annToDataFrame( adata, genescaling = genescaling, layer = layer )    
    
    projected = projectMec(df, mec, glstcol = glstcol )
    projected.columns = [ str(t) for t in projected.columns ]
    if( sigscaling ):
        projected = scale(projected)
        
    if(addon):
        if( len(projected.columns.intersection(adata.obs.columns)) >0 ):
            warnings.warn('signature names already in adata.obs. cannot merge. If wish to remove in obs: \n \
                adata.obs = adata.obs[[t for t in adata.obs.columns if t not in mecnamedict.keys()]  \n or \
                    adata.obs = adata.obs[[t for t in adata.obs.columns if t[:6]!="score_" ]]]'
                )
        else:
            adata.obs = adata.obs.merge( projected , how = 'left', left_index = True, right_index = True)
        adataproj = adata.copy()
    else:
        adataproj = anndata.AnnData( projected, obs = adata.obs )
    return(adataproj)


def projectModuleAnn_aucell( adata, module, glstcol = 'TopGene20'):
    """ 
    Alternative function that projects scRNA data using top genes from MeCs and AUCell
    module has to be list mode.

    Parameters
    ----------
    adata: anndata.AnnData
        Input scanpy object for gene expression
    module: pd.DataFrame
        mec table with one column, two types of format accepted
        format 1, each row is a mec with genes , comma seperated.
    glstcol
        Used only when mec is the list format format 1, and glstcol is the column name to record comma separated gene list for each mec.

    Returns
    ----------
    anndata.AnnData
        adata with aucell score stored in extra columns starting with 'score_sig_' in adata.obs 

    
    """
    
    if(not glstcol): 
        glstcol = module.columns[0]
    for i in module.index:
        gene_list = module.loc[i,  glstcol ].split(',')
        use_gene_list = adata.var_names.intersection( gene_list )
        if(len(use_gene_list)>1):
            sc.tl.score_genes(adata, use_gene_list, score_name = 'score_sig_'+str(i) )
        else:
            adata.obs[ 'score_sig_'+str(i) ] = 0
    return(adata)

"""
### TODO. 
# Tutorial. mapper+annotator.
# Tutorial. differential signature.
# Tutorial functional extraction.
# metatime calling. (code to be wrapped in pipeline with a tutorial, that's it. )
# A class for projected data. this is per-cell class.
#init: score matrix, .obs features, cell location (adata.umap). function and category for features (mecs class.).
#fun:  plot single singature signature using mec name. plot all signature. plot all signature in
# 
#THEN, a seperate annotator function. takes projected class. takes adata input ( overclustering. )
# THEN, a separate differential signature comparison. 
# mecs class: print top gene, plot top gene, extract top tf, extract enrichr, extract lisa ranking, extract lisa-mec ranking plot. 
"""