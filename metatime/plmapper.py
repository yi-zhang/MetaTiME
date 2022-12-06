"""
# function for plotting after mapper
"""
import numpy as np 
import matplotlib.pyplot as plt
import importlib as imp
import pandas as pd
import scanpy as sc
import anndata
from adjustText import adjust_text

# set the colormap and centre the colorbar
import matplotlib.colors as colors

class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar separately for above-midpoint and below-midpoint. Useful for visualizing the signature continuum with diverging colorbar. borrowed from [chris35wills](http://chris35wills.github.io/matplotlib_diverging_colorbar/)
    
    Parameters
    ----------
    colors.Normalize
    e.g. midpoint=0,vmin=vmin,vmax=vmax

    Examples 
    ----------
    # Pass kwargs to scatterplot
    >>> kwargs={'color_map':'RdBu_r', 'norm': MidpointNormalize(midpoint=0,vmin=vmin,vmax=vmax) }
    
    # note: 
    """
    
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def plot_umap_proj( adata : anndata.AnnData, 
                    Xproj : pd.DataFrame, 
                    mecnamedict, 
                   use_MeC_name = True, 
                   figfile = None,
                   N_col = 4 ,
                   N_MECS = 100,
                   ):
    """
    Plotting function: plot the signature continuum score for each funcntional MeC in one large figure.

    Parameters
    ----------
    adata: anndata.AnnData. 
        scanpy object for scRNA data.
        Two formats are accepted.
        1. scanpy object with gene by cell expression, with extra mec columns in adata.obs. 
        When adata.var_names has less than 100 features, go with format1.
        2. scanpy object with only gene by MeC scores in adata.X. projected style format. 

    Xproj: pd.DataFrame
        Dataframe for projected scores


    Returns
    ----------
        matplotlib figure

    
    """
    if(len(adata.var_names)< N_MECS ):
        # adata is a projected format
        Xproj = adata.to_df()
    elif([t for t in Xproj.columns if t not in adata.obs.columns]):
        adata.obs = adata.obs.merge( Xproj,left_index=True,right_index=True,how='left')
    N_col = N_col
    N_row = int(Xproj.shape[1] /N_col)+1
    fig, axes=plt.subplots(nrows= N_row , ncols= N_col ,figsize=(N_col*3,N_row*2), sharex=True,sharey=True, num=86)

    for i in range( Xproj.shape[1] ):
    #for i in range( 8 ):
        i_row = int(i/N_col)
        i_col = i%N_col
        vmin = min( min(Xproj[ Xproj.columns[i] ].values) , -3 )
        vmax = max( max(Xproj[ Xproj.columns[i] ].values) , 3 )
        kwargs={'color_map':'RdBu_r', 'norm': MidpointNormalize(midpoint=0,vmin=vmin,vmax=vmax) }
        # 
        if( use_MeC_name ): 
            #title = str(i)+' ' +mecnamedict[ Xproj.columns[i] ]
            title = Xproj.columns[i].replace('score_','') +' ' +mecnamedict[ Xproj.columns[i] ]
        else:
            title = Xproj.columns[i].replace('score_','') 
        sc.pl.umap(adata, color= Xproj.columns[i], 
                        title = title,
                       show=False, ax = axes[i_row, i_col], **kwargs,)
        axes[i_row, i_col].set_xlabel('')
        axes[i_row, i_col].set_ylabel('')
    if(figfile ):
        plt.savefig( figfile , dpi=200 )
    return( axes )

########## plotting full panel for pdata format ############

## todo: group and rank the components

def plot_umap_proj_pdata( pdata, mecnamedict, 
                   use_MeC_name = True, 
                   figfile = None,
                   N_col = 4 ,
                   ):
    """
    Save all projection images in one large figure. input is pdata format. Deprecated
    """
    pdata = pdata[ : , list(mecnamedict.keys())]
    Xproj = pdata.to_df()
    N_col = N_col
    N_row = int(Xproj.shape[1] /N_col)+1
    fig, axes=plt.subplots(nrows= N_row , ncols= N_col ,figsize=(N_col*3,N_row*2), sharex=True,sharey=True, num=86)

    for i in range( Xproj.shape[1] ):
    #for i in range( 8 ):
        i_row = int(i/N_col)
        i_col = i%N_col
        vmin = min( min(Xproj[ Xproj.columns[i] ].values) , -3 )
        vmax = max( max(Xproj[ Xproj.columns[i] ].values) , 3 )
        kwargs={'color_map':'RdBu_r', 'norm': MidpointNormalize(midpoint=0,vmin=vmin,vmax=vmax) }
        # 
        if( use_MeC_name ): 
            #title = str(i)+' ' +mecnamedict[ Xproj.columns[i] ]
            title = Xproj.columns[i].replace('score_','') +' ' +mecnamedict[ Xproj.columns[i] ]
        else:
            title = Xproj.columns[i].replace('score_','') 
        sc.pl.umap(pdata, color= Xproj.columns[i], 
                        title = title,
                       show=False, ax = axes[i_row, i_col], **kwargs,)
        axes[i_row, i_col].set_xlabel('')
        axes[i_row, i_col].set_ylabel('')
    if(figfile ):
        fig.savefig( figfile , dpi=200 )
    return( fig )


def plot_umap_mec(pdata : anndata.AnnData,
                 meccol: str, 
                 mecnamedict: dict, 
                 use_MeC_name: bool = True, 
                 figfile: str = None,
                 figsize=(3,3)
                 ):

    """
    Plot signature continuum for a specific MeC. 

    Parameters
    ----------
    pdata: anndata.AnnData 
        scanpy object with projection score matrix in pdata.X
    meccol
        MeC id
    mecnamedict


    Returns
    ----------

    """
    Xproj = pdata.to_df()

    fig,ax =plt.subplots(1,figsize= figsize )

    # keep value same range
    vmin = min( min(Xproj[meccol] .values) , -3 )
    vmax = max( max(Xproj[meccol] .values) , 3 )
    kwargs={'color_map':'RdBu_r', 'norm': MidpointNormalize(midpoint=0,vmin=vmin,vmax=vmax) }

    # figure title
    if( use_MeC_name ): 
        title = meccol.replace('score_','') +' ' +mecnamedict[ meccol ]
    else:
        title = meccol.replace('score_','') 

    sc.pl.umap(pdata, color= meccol, title=title,
                show=False, ax=ax, **kwargs,)
    if(figfile ):
        plt.savefig( figfile , dpi=200 )        
    return( fig )




def plot_umap_mecproj_2condition( pdata, meccol, mecnamedict, 
                    cellscond1,cellscond2,
                   use_MeC_name = True, 
                   figfile = '../tmp.png',
                  
                   ):
    """
    
    Pplotting comparison panel for pdata format, beta version.

    Examples
    ----------
    `Xproj = pdata.to_df()
    use_MeC_name = True
    mec_col ='score_1'
    mec_name = mecnamedict[ mec_col ]
    condcol = 'response'
    # cells from 2 conditions
    allcellscond1= adata[adata.obs[ condcol ].isin(['Non-responder']), ].obs.index
    allcellscond2= adata[adata.obs[ condcol ].isin(['Responder']), ].obs.index

    plmapper.plot_umap_mecproj_2condition( pdata, meccol = 'score_0', mecnamedict = mecnamedict, 
                                        cellscond1 = allcellscond1, cellscond2 = allcellscond2, 
                                        use_MeC_name = True, )
                                     
    `
    """
    Xproj = pdata.to_df()

    if((len(cellscond1)==0) or(len(cellscond2)==0)):
        return(None)
    # plot panel
    fig, axes=plt.subplots(nrows= 1 , ncols= 2 ,figsize=(2*3,1*3), sharex=True,sharey=True, num=2)

    # keep value same range
    vmin = min( min(Xproj[meccol] .values) , -3 )
    vmax = max( max(Xproj[meccol] .values) , 3 )
    kwargs={'color_map':'RdBu_r', 'norm': MidpointNormalize(midpoint=0,vmin=vmin,vmax=vmax) }
    """
    # figure title
    if( use_MeC_name ): 
        title = meccol.replace('score_','') +' ' +mecnamedict[ meccol ]
    else:
        title = meccol.replace('score_','') 
        """
    # condition1,2
    sc.pl.umap(pdata[ pdata.obs.index.intersection( cellscond1)  ], color= meccol, 
                    title = 'cond1',
                show=False, ax = axes[0], **kwargs,)

    sc.pl.umap(pdata[ pdata.obs.index.intersection( cellscond2) ], color= meccol, 
                    title = 'cond2',
                show=False, ax = axes[1], **kwargs,)
    axes[0].set_xlabel('')
    axes[1].set_ylabel('')
    if(figfile ):
        plt.savefig( figfile , dpi=200 )        
    return( fig )





############# annotator plots #############


def gen_mpl_labels(
    adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None
):
    """ 
    Get locations of cluster median . Borrowed from scanpy github forum.
    """
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)

def plot_annotation_on_data( 
    sdata, 
    COL='MetaTiME'  ,
    title = None, 
    fontsize = 8,
    MIN_CELL = 5,
    ):
    """ 
    Plot annotated cells with  non-overlapping fonts.
    Calls: gen_mpl_labels
    Occasionally, if sdata.uns['MetaTiME_overcluster_colors'] exists, color may not be updated. In that case, resetting color  `del pdata.uns['MetaTiME_overcluster_colors']`
    
    Parameters
	----------
	sdata : anndata.AnnData
        scanpy object for single cell , or scanpy object with projected signature. 
    COL : str
        column of cluster assignment in sdata.obs
    MIN_CELL: int
        Minimum number of cells to plot and mark

	Returns
	----------
        Matplotlib figure.

    Examples 
    ----------
    >>> plmapper.plot_annotation_on_data(pdata, title = 'MetaTiME')
    
    """
    #
    if not title:
        title = COL

    if( MIN_CELL >0 ):
        groupcounts = sdata.obs.groupby( COL ).count()
        groupcounts = groupcounts[groupcounts.columns[0]]
        group_with_good_counts = groupcounts[groupcounts>= MIN_CELL ].index.tolist()
        sdata = sdata[ sdata.obs[ COL ].isin( group_with_good_counts ) ]

    with plt.rc_context({"figure.figsize": (6, 6), "figure.dpi": 300, "figure.frameon": False}):
        #ax = sc.pl.umap(pdata, color="MetaTiME_overcluster", show=False, legend_loc=None, frameon=False, size=30)
        ax = sc.pl.umap(sdata, color= COL , show=False, legend_loc=None, add_outline=False, 
                #legend_loc='on data',legend_fontsize=6, legend_fontoutline=2,
                title= title, 
                palette=plt.cycler("color",plt.cm.tab20(np.linspace(0,1,20))), 
                #palette=plt.cycler("color",plt.cm.Set1(np.linspace(0,1,9))), 
                    )
        gen_mpl_labels(
            sdata,
            COL,
            exclude=("None",),  
            ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
            text_kwargs=dict(fontsize= fontsize ,weight='bold'),
        )
        fig = ax.get_figure()
        fig.tight_layout()
        #plt.show()
        return( fig,ax )
        #plt.savefig('MetaTiME_overcluster.pdf', dpi=200)


############## snippets ##########
"""
# save all projection images in single figures.

for i in range( Xproj.shape[1] ):
#for i in range( 2 ):
    vmin = min( min(Xproj[ Xproj.columns[i] ].values) , -2 )
    vmax = max( max(Xproj[ Xproj.columns[i] ].values) , 2 )
    kwargs={'color_map':'RdBu_r', 'norm': plmapper.MidpointNormalize(midpoint=0,vmin=vmin,vmax=vmax) }
    ax = sc.pl.umap(adata, color= Xproj.columns[i], **kwargs,
                    title = mecnamedict[ Xproj.columns[i] ],
                   show=False )
#    plt.savefig('../analysis/20211001_scICA/BRCA_Alex/BRCA_Alex.rawumap.mec_'+str(keyi)+'.png')
    plt.close()
    
    """
