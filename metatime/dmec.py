"""
## calculate differential mec signature comparing two conditions
"""
import numpy as np
import matplotlib.pyplot as plt
import importlib as imp
from typing import Optional, Union
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata
import scipy.stats as ss

def _get_score_pd( adata, meccols ):
    """ 
    Obtain scores from scanpy object.
    compatible with:
        scanpy object with MetaTiME projection score in .X .
        or,
        scanpy object with MetaTiME projection score saved in .obs[meccols]
    """
    if meccols[0] in adata.var_names:
        x = adata.to_df()[meccols]
    elif( meccols[0] in adata.obs.columns):
        x = adata.obs[ meccols ]
    else:
        return(None)
    return(x)

def count_cell_condition( obs,  countbycols : list):
    """
    Counting cells by condition 

    Params
    ----------
    obs: adata.obs,
        metadata dataframe
    countbycols: list, has to be columns in obs.

    Examples
    ----------
    condcol = 'treatment' 
    clustercol = 'MetaTiME' 
    samplecol = 'patient'
    counts = count_cell_condition( adata,obs, countbycols  = [condcol, clustercol, samplecol] )

    """
    cellcts = obs.groupby( countbycols ).count()['MetaTiME_overcluster'].fillna(0).unstack()
    return(cellcts)

def enrich_mec_cluster(  labeltype, adata, labelcol , meccols , mecnamedict ):
    """ 
    For one cluster, Output which MeC is enriched.
    criterion: MeC score cluster-mean passes zscore 1 (1-standard deviation)
    
    Parameters
    ----------
    adata: anndata.AnnData.
        compatible with:
        scanpy object with MetaTiME projection score in .X .
        or,
        scanpy object with MetaTiME projection score saved in .obs[meccols]

    """
    if meccols[0] in adata.var_names:
        x = adata[adata.obs[ labelcol ]==labeltype].to_df()[meccols]
    elif( meccols[0] in adata.obs.columns):
        x = adata.obs.loc[ adata.obs[ labelcol ]==labeltype,meccols ]
    else:
        return(False)
    tscores = x.mean(axis=0).sort_values(ascending=False)
    #tscores = x.median(axis=0).sort_values(ascending=False)
    tscores = tscores[tscores>=1]#  [:2] # top 2 pass 1std
    tscores = tscores.rename(mecnamedict)
    topmecs = tscores.index.tolist()
    
    return(tscores.to_dict())

def enrich_mec( adata , labelcol , mecnamedict ):
    """ 
    Determine enriched MeCs per group or cluster.

    Parameters
    ----------
    adata: anndata.AnnData.
        compatible with:
        scanpy object with MetaTiME projection score in .X .
        or,
        scanpy object with MetaTiME projection score saved in .obs[meccols]

    labelcol: str
        The column name of grouping cells
    mecnamedict
        for renaming dat columns into functional names
        can be loaded from pre-computed tumor MeC functional annotation

    Returns
    ----------
    pd.DataFrame
        A data frame with each cluster in a row, and a column 'topmec' listing MeCs with significant enrichment as candidate for testing differential signature.

    Example
    ----------
    `
        cluster_mec_enriched = dmec.enrich_mec( labeltypes, adata=adata, labelcol =labelcol, meccols = meccols, mecnamedict=mecnamedict)
    `
    """
    # columns for meta components
    meccols = list(mecnamedict.keys())
    # labeltypes: categories in the cell grouping column 
    labeltypes = adata.obs[labelcol].drop_duplicates().tolist()
    
    res1 = pd.DataFrame( index=labeltypes)
    res={}
    for labeltype in labeltypes:
        topmec = enrich_mec_cluster( labeltype ,  adata=adata, labelcol =labelcol, meccols = meccols , mecnamedict =mecnamedict) 
        res.update( {labeltype: topmec})
        res1.loc[labeltype,'topmec']= ','.join(topmec)
    return(res1)

def cohend( x1, x2 ):
    """ t  test effect size """
    cohens_d = (np.mean(x1) - np.mean(x2)) / (np.sqrt((np.std(x1) ** 2 + np.std(x2) ** 2) / 2))
    return(cohens_d)

def test_sig_smallcluster(Xproj, meci, cellscond1, cellscond2):
    """ 
    testing signature difference at each cell group, for certain mec signatures in meci.
    """
    
    x1 = Xproj.loc[cellscond1, meci].values
    x2 = Xproj.loc[cellscond2, meci].values

    #cut=2
    #y1=x1[x1>=cut]
    #y2=x2[x2>=cut]
    #[s,p]=ss.ranksums(x2, x1)
    [s,p]=ss.ttest_ind(x2, x1)
    #u,p = ss.mannwhitneyu(x1,x2)
    #logfc = np.log2((np.mean(y2)/len(x2))/(np.mean(y1)/len(x1)))# consider cell number differences
    #logfc = np.log2((np.median(y2)+1)/(np.median(y1)+1))# consider cell number differences
    #logfc = np.log2( (np.median(x2) +1 )/(np.median(x1) +1 ) )
    effect = (np.mean(x2)  - np.mean(x1) )  
    
    #logfc = cohend( x1, x2)
    #logfc = np.log2((np.mean(x2))/(np.mean(x1)))
    meanb=np.mean(x2); meana=np.mean(x1)
    #logctfc = np.log2((len(y2)/len(x2))/(len(y1)/len(x1)))# consider cell number differences
    #logctfc = np.log2((np.sum(y2)/len(x2))/(np.sum(y1)/len(x1)))# consider cell number differences
    #logctfc = np.log2( (np.sum(x2)/N2 +1)/(np.sum(x1)/N1 +1))# consider cell number differences
    effect1 = effect

    logp=-np.log10(p)
    return( logp, effect,effect1,meanb,meana )
### plot enriched 1st mec

def plot_top1mec( adata, mec_enriched , legend_loc ='on data'):
    mec_enriched['top1mec'] = mec_enriched['topmec'].str.split(',').str.get(0)
    adata.obs = adata.obs.merge( mec_enriched[['top1mec']], left_on = 'celltype_subset', right_index=True, how='left')
    sc.pl.umap(adata, color = ['top1mec'], legend_fontsize='5', legend_loc= legend_loc , alpha=0.4 )


def dmec( adata, 
          pdata,
          allcellscond1,
          allcellscond2,
          clustercol: str  ,
          cluster_mec_enriched: pd.DataFrame, 
          mecnamedict,
          test_clusters = None ,):
    """ 
    Differential signature analysis

    Parameters
    ----------
    adata
        scRNA object
    pdata
        Projection score table, cell by signature. 
    allcellscond1
        list of barcodes of cells from condition 1
    allcellscond2
        list of barcodes of cells from condition 2
    clustercol
        cell cluter column in adata.obs 
    cluster_mec_enriched
        what mec to test in each cluster in test_clusters
    mecnamedict
        for renaming dat columns into functional names
        can be loaded from pre-computed tumor MeC functional annotation
    test_clusters
        Which cell clusters to test based on adata.obs[clutsercol]

    Returns
    ----------
    diffmec
        dictonary containing per-cluster differential testing results, good for summarizing and generating plots
    diffmecsig
        dictonary containing per-cluster differential testing results filtered to contain significant ones
    diffmec_full
        dictonary containing per-cluster differential testing results with more details

    """
    if(not test_clusters or test_clusters == 'all'):
        test_clusters =adata.obs[clustercol].drop_duplicates().tolist()
    diffmec_full={}
    diffmec={}
    diffmecsig={}
    for test_cluster in test_clusters:
        cell_clusteri = adata[adata.obs[ clustercol ]== test_cluster ].obs.index
        cellscond1 = allcellscond1.intersection(cell_clusteri)
        cellscond2 = allcellscond2.intersection(cell_clusteri)
        if((len(cellscond1)==0) or(len(cellscond2)==0)):
            continue
        # signatures that has high (>1std) score for this cluster
        clusters_strong_sig1 = cluster_mec_enriched.loc[test_cluster,'topmec'].split(',')
        if(len(clusters_strong_sig1)==0 or clusters_strong_sig1[0]==''):
            continue


        Xproj = pdata.to_df()
        res = pd.DataFrame(index=Xproj.columns)
        #N2 = len(allcellscond2);N1 = len(allcellscond1)
        for meci in Xproj.columns:
            t = test_sig_smallcluster(Xproj, meci, cellscond1, cellscond2 )
            res.loc[meci, '-logp']= t[0]
            res.loc[meci, 'effect']= t[1]
            res.loc[meci, 'effect1']= t[2]
            res.loc[meci, 'meanb']= t[3]
            res.loc[meci, 'meana']= t[4]

        diff_mec = res.sort_values(ascending=False,by='effect')
        diff_mec['scorecol'] = diff_mec.index
        diff_mec = diff_mec.rename(index=mecnamedict)

        diff_mec_matchedmec = diff_mec.loc[clusters_strong_sig1]
        diff_mec_matchedmec_sig = diff_mec_matchedmec[diff_mec_matchedmec['-logp']>=2].sort_values(ascending=False, by='effect')
        tmp = {'diff_mec_matchedmec_sig': diff_mec_matchedmec_sig,
               'diff_mec': diff_mec,
               'full_diff': res}
        diffmec_full.update({test_cluster : tmp })
        if(len(diff_mec_matchedmec_sig)>0):
            diffmecsig.update({ test_cluster: diff_mec_matchedmec_sig })
        if(len(diff_mec)>0):
            diffmec.update({ test_cluster: diff_mec_matchedmec })
    return( diffmec, diffmecsig, diffmec_full)


def plotdata_diffmec( diffmec, mecnamedict):
    """
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    from plmapper import MidpointNormalize
    mat_p = pd.DataFrame(index = diffmec.keys(), columns = mecnamedict.values())
    mat_eff = pd.DataFrame(index = diffmec.keys(), columns = mecnamedict.values())
    mat_b = pd.DataFrame(index = diffmec.keys(), columns = mecnamedict.values())
    mat_a = pd.DataFrame(index = diffmec.keys(), columns = mecnamedict.values())
    for i in diffmec.keys():
        #print(i)
        mat_p.loc[i, diffmec[i].index] = diffmec[i]['-logp'].values
        mat_b.loc[i, diffmec[i].index] = diffmec[i]['meanb'].values
        mat_a.loc[i, diffmec[i].index] = diffmec[i]['meana'].values
        mat_eff.loc[i, diffmec[i].index] = diffmec[i]['logctfc'].values
    cluster_order = mat_p.max(axis=1).sort_values(ascending=False).index
    
    mec_order = mat_p.max(axis=0).sort_values(ascending=False).index
    tmp = pd.DataFrame(mat_p.loc[mat_p.index.intersection(cluster_order), mat_p.columns.intersection( mec_order) ].stack().rename('-logp'))
    #tmp = pd.DataFrame(mat_p.loc[cluster_order, mec_order].stack().rename('-logp'))
    tmp=tmp.merge( pd.DataFrame(mat_eff.stack().rename('effect')) , left_index=True,right_index=True, how='left')
    tmp=tmp.merge( pd.DataFrame(mat_b.stack().rename('meanb') ) , left_index=True,right_index=True, how='left')
    tmp=tmp.merge( pd.DataFrame(mat_a.stack().rename('meana') ) , left_index=True,right_index=True, how='left')
    tmp = tmp.reset_index()
    return(tmp)

def plot_diffmec( diffmec, mecnamedict):
    """
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    from metatime import plmapper
    from metatime.plmapper import MidpointNormalize
    mat_p = pd.DataFrame(index = diffmec.keys(), columns = mecnamedict.values())
    mat_eff = pd.DataFrame(index = diffmec.keys(), columns = mecnamedict.values())
    mat_b = pd.DataFrame(index = diffmec.keys(), columns = mecnamedict.values())
    mat_a = pd.DataFrame(index = diffmec.keys(), columns = mecnamedict.values())
    for i in diffmec.keys():
        #print(i)
        mat_p.loc[i, diffmec[i].index] = diffmec[i]['-logp'].values
        mat_b.loc[i, diffmec[i].index] = diffmec[i]['meanb'].values
        mat_a.loc[i, diffmec[i].index] = diffmec[i]['meana'].values
        mat_eff.loc[i, diffmec[i].index] = diffmec[i]['effect'].values
    cluster_order = mat_p.max(axis=1).sort_values(ascending=False).index
    
    mec_order = mat_p.max(axis=0).sort_values(ascending=False).index
    tmp = pd.DataFrame(mat_p.loc[mat_p.index.intersection(cluster_order), mat_p.columns.intersection( mec_order) ].stack().rename('-logp'))
    #tmp = pd.DataFrame(mat_p.loc[cluster_order, mec_order].stack().rename('-logp'))
    tmp=tmp.merge( pd.DataFrame(mat_eff.stack().rename('effect')) , left_index=True,right_index=True, how='left')
    tmp=tmp.merge( pd.DataFrame(mat_b.stack().rename('meanb') ) , left_index=True,right_index=True, how='left')
    tmp=tmp.merge( pd.DataFrame(mat_a.stack().rename('meana') ) , left_index=True,right_index=True, how='left')
    tmp = tmp.reset_index()

    sns.set_theme(style="whitegrid")
    kwargs={'hue_norm': MidpointNormalize(midpoint=0,vmin=min(tmp['effect']), vmax = max(tmp['effect']) ) ,
            }
    g = sns.relplot(data = tmp, 
                x = 'level_0', y='level_1', hue = 'effect', size = '-logp',
                palette="RdBu_r",
                #hue_norm=(-1, 1), 
                edgecolor=".7",
        height=8, aspect = 1.2, 
        #width=8, 
        sizes=(20, 100), size_norm=(-.2, .8),
                **kwargs
    )
    g.set(xlabel="", ylabel="", aspect="equal")
    g.despine(left=True, bottom=True)
    g.ax.margins(.02)
    for label in g.ax.get_xticklabels():
        label.set_rotation(90)
    for artist in g.legend.legendHandles:
        artist.set_edgecolor(".7")
    return(g)


def topdiff_df(diffmec,
             CUT_LOGPPOS = 1.3,
             CUT_EFFECTPOS = 1,
             CUT_EFFECTNEG = -1,
             FCCOL = 'effect'
             ):
    """
    Organize per-cluster differential signature results into tables

    Parameters
    ----------
    diffmec
        differential expression dict generated by dmec.dmec
    CUT_LOGPPOS 
        negative log p value cutoff to determin significance, default: 1.3 (#p=0.5)
    CUT_EFFECTPOS
        positive effect size to cut . default: 1 ( 1 std shift)
    CUT_EFFECTNEG
        negative effect size to cut . default: -1 ( 1 std shift)
    FCCOL
        the 'effect' to use (fold change)

    Example
    ----------
    `diff, diff_sig = dmec.topdiff_df(diffmec)`

    """
    tlst = []
    for cluster1 in diffmec.keys():
        df = diffmec[ cluster1 ].copy()
        df['cluster1']=cluster1
        df['mec']=df.index
        df.index= df.index + '@' + cluster1
        # MHCI@Mac_IL1B; Mac_IL1B (not Mac_IL1B@Mac_IL1B)
        df.index=df.index.str.replace( cluster1+'@'+cluster1 , cluster1 )
        tlst.append(df)
    dat = pd.concat(tlst)
    sigpos = dat[( dat[ FCCOL ]>= CUT_EFFECTPOS ) & ( dat['-logp']>= CUT_LOGPPOS ) & ( (dat['meanb']>=1) )]
    #sigpos = dat[( dat[ FCCOL ]>= CUT_EFFECTPOS ) & ( dat['-logp']>= CUT_LOGPPOS ) ]
    signeg = dat[( dat[ FCCOL ]<= CUT_EFFECTNEG ) & ( dat['-logp']>= CUT_LOGPPOS ) & ( (dat['meana']>=1) )]
    #signeg = dat[( dat[ FCCOL ]<= CUT_EFFECTNEG ) & ( dat['-logp']>= CUT_LOGPPOS ) ]
    markpos = sigpos.sort_values(ascending=False, by= FCCOL ).index
    markneg = signeg.sort_values(ascending=True, by= FCCOL ).index
    dat_sig = dat.loc[markpos.tolist() + markneg.tolist()].sort_values(ascending=False, by ='-logp')
    return( dat, dat_sig )


def plot_topdiff( diffmec, fontsize = 5 ):
    """
    Parameters
    ----------
    diffmec
        dictonary containing per-cluster differential testing results
    fontsize
        fontsize on the differential scatter plot

    Example
    ----------
    `sns.set_theme(style='white')
    fig_topdiff = dmec.plot_topdiff(diffmec, fontsize=6)
    `
    """
    CUT_LOGPPOS = 1.3 #p=0.5
    CUT_EFFECTPOS = 1 # 1 std shift
    CUT_EFFECTNEG = -1 # 1 std shift
    #FCCOL = 'logctft'
    FCCOL = 'effect'
    tlst = []
    for cluster1 in diffmec.keys():
        df = diffmec[ cluster1 ].copy()
        df['cluster1']=cluster1
        df['mec']=df.index
        df.index= df.index + '@' + cluster1
        # MHCI@Mac_IL1B; Mac_IL1B (not Mac_IL1B@Mac_IL1B)
        df.index=df.index.str.replace( cluster1+'@'+cluster1 , cluster1 )
        tlst.append(df)
    dat = pd.concat(tlst)

    fig,ax = plt.subplots(1, figsize=(6,5))
    sns.set_theme(style="white")
    sns.scatterplot(data=dat, x= FCCOL ,y='-logp', color='gray',legend=False,alpha=0.3)
    sigpos = dat[( dat[ FCCOL ]>= CUT_EFFECTPOS ) & ( dat['-logp']>= CUT_LOGPPOS ) & ( (dat['meanb']>=1) )]
    signeg = dat[( dat[ FCCOL ]<= CUT_EFFECTNEG ) & ( dat['-logp']>= CUT_LOGPPOS ) & ( (dat['meana']>=1) )]
    markpos = sigpos.sort_values(ascending=False, by= FCCOL ).index
    markneg = signeg.sort_values(ascending=True, by= FCCOL ).index

    # tuning figure
    lim = np.max([abs(t) for t in ax.get_xlim()])
    ax.set_xlim([-lim,lim])
    #ax.set_ylim([-30, ax.get_ylim()[1]])
    ax.axvline( CUT_EFFECTPOS ,linestyle='--', color='gray', alpha=0.4)
    ax.axvline( CUT_EFFECTNEG,linestyle='--', color='gray', alpha=0.4)
    ax.axvline( 0,linestyle='--', color='gray', alpha=0.6)
    ax.axhline( CUT_LOGPPOS ,linestyle='--', color='gray', alpha=0.4)
    # color
    sns.scatterplot(data=sigpos, x= FCCOL ,y='-logp', color='brown', size='meanb',legend=False, ax=ax)
    sns.scatterplot(data=signeg, x= FCCOL ,y='-logp', color='steelblue', size='meana',legend=False, ax=ax)
    # text
    texts=[]
    for meci in markpos:
        texts.append( plt.text( sigpos.loc[meci, FCCOL ],  sigpos.loc[meci,'-logp'], meci.replace('@','@') , fontsize= fontsize, color='red'))
    for meci in markneg:
        texts.append( plt.text( signeg.loc[meci, FCCOL ],  signeg.loc[meci,'-logp'], meci.replace('@','@') , fontsize= fontsize, color='blue'))

    from adjustText import adjust_text
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'), autoalign='xy')
    ax.set_xlabel('Effect')
    return( fig, ax )

