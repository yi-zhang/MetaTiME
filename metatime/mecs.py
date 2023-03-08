"""
## class and functions for meta components
"""
import numpy as np
import pandas as pd
import os
import sys

from metatime import config

################  Special manual annotation of MetaTiME precomputed MeCs #############

def load_mecname( mode='table', 
                mecDIR= config.SCMECDIR
                ):
    """
    Load functional annotation for pre-computed tumor microenvironemnt MeCs

    Parameters
    ----------
    mode
        choose from ['mecnamedict', 'table', 'meciddict']
        load manual assigned name for easy understanding of assigned names
        from file: MeC_anno_name.tsv under mecDIR.
        Required columns: `['MeC_id', 'Annotation', 'UseForCellStateAnno']` 
        Required seperator: tab
        Annotation column NA will be filtered.

    Returns
    ----------
    functional annotation in desided format

    """
    mecname = pd.read_table(os.path.join( mecDIR , 'MeC_anno_name.tsv'))
    if(mode=='table'):
        mectable=mecname.copy()
        #mectable['color_level0'] = mectable['MajorLineage_level0'].apply(lambda t: config.level0colordict[t])
        return(mectable)        
    else:
        mask = (~mecname['Annotation'].isna())
        #mecname.loc[mask,'Annotation']=mecname.loc[mask,'MeC_id' ]
        mecnamedict = mecname[mask].set_index('MeC_id')[['Annotation']].to_dict()['Annotation']

    if(mode=='mecnamedict'):
        return(mecnamedict)

def getmecnamedict_ct( mectable, only_include_mecs_UseForCellStateAnno = True ):
    """ 
    Collect list of meta components to be used for cell state annotation. 

    Parameters
    ----------
    mectable
        must have `UseForCellStateAnno` column with 0 or 1. 1:used in cell state annotation. this is helpful to remove pan-cell cell state like general mitochondrial activity component.
    only_include_mecs_UseForCellStateAnno
        filter mecs used for cell state annotation. 
        Default: True 
    Returns
    ----------
    dict 
        a subsetted dictionary only containing MeCs used for cell state enrichment.

    """
    if(only_include_mecs_UseForCellStateAnno):
        mecnamedict = mectable[mectable['UseForCellStateAnno']==1].set_index('MeC_id')['Annotation'].to_dict()
    else:
        mecnamedict = mectable.set_index('MeC_id')['Annotation'].to_dict()
    return(mecnamedict)



class MetatimeMecs():
    """
    Class for MetaTiME Mecs.
    
    Parameters
    ----------
    modelpath
        directory of model files

    Attributes
    ----------
    mecscore
        pandas dataframe, z-weights  gene by component
    mectopg
        pandas dataframe, list for top gene only
    mecanno
        pandas dataframe, annotation for each mec. 

    """
    def __init__(self, mec_score, mec_topg, mec_anno ) :
        self.mec_score = mec_score
        self.mec_topg = mec_topg
        self.mec_anno = mec_anno

    @staticmethod
    def load_mec_precomputed(  mecDIR = config.SCMECDIR, ):
        """
        Load pre-computed Meta-component matrix and ordered list.
        Look for precomputed files in mecDIR

            MeC_allgene_average-weights.tsv
            MeC_topgene.tsv

        Parameters
        ----------
        mecDIR
            Path for model files

        Returns
        ----------
        MetatimeMecs
            Loaded MeCs    

        Examples 
        ----------
        >>> mecmodel = mecs.MetatimeMecs.load_mec_precomputed()

        """
        try:
            mec_score = pd.read_table( os.path.join( mecDIR , 'MeC_allgene_average-weights.tsv') , index_col = 0)
            mec_topg = pd.read_table( os.path.join( mecDIR, 'MeC_topgene.tsv') , index_col = 0)
            mec_anno = pd.read_table( os.path.join( mecDIR, 'MeC_anno.tsv') , index_col = 0)
            # use 'MeC' name
            mec_score.columns = 'MeC_'+mec_score.columns
            return( MetatimeMecs( mec_score, mec_topg, mec_anno ) )

        except Exception as exception:
            raise Exception(f' Error loading model. check files in {mecDIR}. MeC_allgene_average-weights.tsv, MeC_topgene.tsv, MeC_anno.tsv. {exception}')

    @property
    def feature(self) -> np.ndarray:
        """ get genes covered in the mec z-weight matrix """
        return( self.mec_score.index)
    
    @property
    def nmecs(self) -> np.float:
        """ get number of meta-components """
        return( self.mec_score.columns.shape[0])


    # interpretation functions
