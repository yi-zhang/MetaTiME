## function for meta components
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
        from file: MeC_anno_name.txt under mecDIR.
        Required columns: `['MeC_id', 'Annotation', 'UseForCellTypeAnno']`
        Annotation column NA will be filtered.

    Returns
    ----------
    functional annotation in desided format

    """
    mecname = pd.read_table(os.path.join( mecDIR , 'MeC_anno_name.txt'))
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

def getmecnamedict_ct( mectable ):
    """ 
    Collect list of meta components to be used for cell state annotation. 

    Parameters
    ----------
    mectable
        must have `UseForCellTypeAnn` column with 0 or 1. 1:used in cell state annotation. this is helpful to remove pan-cell cell state like general mitochondrial activity component.
        
    Returns
    ----------
    dict 
        a subsetted dictionary only containing MeCs used for cell state enrichment.

    """
    mecnamedict = mectable[mectable['UseForCellTypeAnno']==1].set_index('MeC_id')['Annotation'].to_dict()
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

            MeC_allgene_average-weights.txt
            MeC_topgene.txt

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
            mec_score = pd.read_table( os.path.join( mecDIR , 'MeC_allgene_average-weights.txt') , index_col = 0)
            mec_topg = pd.read_table( os.path.join( mecDIR, 'MeC_topgene.txt') , index_col = 0)
            mec_anno = pd.read_table( os.path.join( mecDIR, 'MeC_anno.txt') , index_col = 0)
            # use 'MeC' name
            mec_score.columns = 'MeC_'+mec_score.columns
            return( MetatimeMecs( mec_score, mec_topg, mec_anno ) )

        except Exception as exception:
            raise Exception(f' Error loading model. check files in {mecDIR}. MeC_allgene_average-weights.txt, MeC_topgene.txt, MeC_anno.txt. {exception}')

    @property
    def feature(self) -> np.ndarray:
        """ get genes covered in the mec z-weight matrix """
        return( self.mec_score.index)
    
    @property
    def nmecs(self) -> np.float:
        """ get number of meta-components """
        return( self.mec_score.columns.shape[0])


    # interpretation functions
