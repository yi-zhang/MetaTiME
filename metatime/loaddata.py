""" Functions for new scRNA data processing """
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from typing import Optional, Union

###  function  ###



def isRawCount( ds ):
    """ 
    Tell whether the matrix is raw  count or float 
        if top 100 cells all mod 0, then determine it's raw count = integer

    Parameters
    ----------
    ds
        loom, or just a numpy matrix 2d. 
    
    Returns
    ----------
    int
        0, not raw count.
        1, is raw count.

    """
    if( ds[:100,:100].sum() - int( ds[:100,:100].sum() ) == 0 ) :
        datatype = 'rawcount'
        datatype = 1
    else:
        datatype = 'notcount'
        datatype = 0
    return(datatype)

### pre-processing adata

def adatapp(
    adata_input : anndata.AnnData, 
    mode: str = 'pp', # or  'umap'
    random_state = 42,
    MAX_MITO_PERCENT=5,
    MINGENES =500,
    MINCOUNTS = 1000,
    MINCELLS=5,
):

    """ 
    Pre-processing. 

    Parameters
    ----------
    adata_input
        scanpy scRNA object
        if adata.X is integer , or adata has 'counts' layer, like from 10X data.
        go through standard preprocessing.
        remove cells with <500 genes and <1000 counts. 
        remove genes with min cells < 5. 
        keep mitocondrial <5%. 
        normalize to 1e6. log. pca, umap. normalized value saved to 'norm_data' layer

        if adata.X is continuous, such as Smartseq data
        remove cells with <500 genes. remove genes with <5 cells. 
    mode
        choose from ['pp' , 'umap']
        pp: full pre-processing
        umap: preprocessing was done. only compute highly_varaible genes, pca, neighbors, umap.
    random_state
        umap random state
    MAX_MITO_PERCENT
    MINGENES 
    MINCOUNTS
    MINCELLS

    Returns
    ----------
        processed scanpy object

    Examples
    ----------
    >>> adata = adatapp(adata )

    """

    adata = adata_input.copy()

    if(mode=='pp'):
        isCount = isRawCount( adata.X )
        if( isCount ): 
            adata.layers["counts"] = adata.X.copy()
            print( '[Log] adata.X is count. Copy to .layers["counts"] . Pre-processing... ')
            adata.var_names_make_unique()
            adata.obs['n_counts'] = adata.X.sum(1).A1
            adata.obs['log_counts'] = np.log(adata.obs['n_counts']) #
            adata.obs['n_genes'] = (adata.X > 0).sum(1)
            print( '[Log] filtering ')
            sc.pp.filter_cells(adata, min_genes= MINGENES ) 
            sc.pp.filter_cells(adata, min_counts = MINCOUNTS )
            sc.pp.filter_genes(adata, min_cells = MINCELLS) 
            print('[Log] Human MT-percentage')
            adata.var['mt'] = adata.var_names.str.startswith('MT-') 
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            adata = adata[adata.obs.pct_counts_mt < MAX_MITO_PERCENT, :]
            sc.pp.normalize_total(adata, target_sum = 1e4)  # match NormalizeData in Seurat
            sc.pp.log1p(adata) # checked the values are logged.
            sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000 )
            sc.tl.pca(adata,  svd_solver='arpack')
            sc.pp.neighbors(adata )
            sc.tl.umap(adata, random_state = random_state)
            
        else: 
            # directly use the values if X not continuous: already normalized
            print( '[Log] adata.X is continuous . Pre-processing... ')
            adata.var_names_make_unique()
            adata.layers["norm_data"] = adata.X
            sc.pp.filter_cells(adata, min_genes=500) 
            sc.pp.filter_genes(adata, min_cells=5) 
            sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000 )
            sc.tl.pca(adata,  svd_solver='arpack')
            sc.pp.neighbors(adata )
            sc.tl.umap(adata, random_state = random_state)

    elif(mode=='umap'):
        
        sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000) 
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata, random_state = random_state)
    return(adata)


def load(
    file, 
    preprocessing = False,
    delimiter = ',',
):
    """ 
    Read in a expression matrix table file as scanpy object.

    If converting from an Seurat object, table can be saved in R:
        tmpdata = SeuratObject@assays$RNA@data
        write.csv( t( as.matrix(tmpdata) ), file = 'data.csv')

    Parameters
    ----------
    file
        file path for the expression file.
        if suffix is .h5 or .h5ad, data is loaded using scanpy.
        if not, load the file in format of a table text file
    preprocessing
        whether to do preprocessing for the loaded data
    delimiter
        delimiter of the expression matrix table, if file is a text table
    
    Returns
    ----------
    adata
        loaded scRNA data in scanpy object
        
    Examples
    ----------
    >>> adata = load( file, preprocessing=False )
    
    """
    adata = None

    suffix = file.split('.')[-1]
    if(suffix == 'h5' or suffix == 'h5ad'):
        annfile = file
        adata = sc.read(annfile)
    else:
        matfile = file
        adata = anndata.read_csv( matfile, delimiter = delimiter, first_column_names = True)
    if(not adata):
        return(None)

    if ( preprocessing == True ):
        #if 'norm_data' not in adata.layers:
        ## do preprocessing from raw
        print( ' Redo preprocessing ' )
        adata = adatapp( adata )
    else:
        print( ' Loaded ' + file)
    return(adata)



def batchharmonize(
    adata , 
    batchcols = [],
    random_state= 0
):
    """
    Harmonize batches such as patient, sample, using Harmony. Re-calculated neighbors, umap .

    Parameters
    ----------
    adata:
        scanpy object. 
    batchcols
        batch columns in adata.obs
    random_state
        umap random_state
    Returns
    ----------
    adata with batch correction and neighbors,umap calcualted

    """
    

    if((batchcols is not None) or (len(batchcols))>=1):
        print(f'[Log ] hamonize with batch ', batchcols)
        if('X_pca' not in adata.obsm.keys()):
            sc.tl.pca(adata,  svd_solver='arpack')
        try:
            sc.external.pp.harmony_integrate(adata, key = batchcols )
            sc.pp.neighbors(adata, use_rep = 'X_pca_harmony')
        except:
            print('[Log ] harmonization problem, skipping batch correction')
            sc.pp.neighbors(adata)
    else:
        print(f'[Log ]  no batch label detected. Computing neighbors using default')
        sc.pp.neighbors(adata)

    print(f'[Log ] Recomputing umap ' )
    sc.tl.umap(adata, random_state = random_state)
    return(adata)



"""
def isIndexNumberButNotGeneName( adata ):

    if( adata.var_names[0]=='0' ):
        return(True)
    else:
        return(False)
"""


def add_extra_metainfo( adata, meta ):
    """
    Add special extra metainformation for sample datasets

    Parameters
    ----------
    adata: anndata.AnnData
        scanpy object with scRNA expression
    meta: pd.DataFrame
        Dataframe with columns of metadata for each cell. index is cell barcode same as in adata.obs.index.
    
    Returns
    ----------
    anndata.AnnData
        scanpy object with new columns merged in obs. 

    Examples
    ----------
    >>> adata = testdata.add_extra_metainfo( adata, meta )
    """
    # add warning
    adata.obs = adata.obs.join( meta[[col for col in meta.columns if col not in adata.obs.columns]] )
    return(adata)