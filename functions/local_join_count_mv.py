import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.base import BaseEstimator

class Local_Join_Count_MV(BaseEstimator):

    """Global Spatial Pearson Statistic"""

    def __init__(self, connectivity=None):
        """
        Initialize a Join_Counts_Local estimator
        Arguments
        ---------
        connectivity:   scipy.sparse matrix object
                        the connectivity structure describing the relationships
                        between observed units. Will be row-standardized.
        Attributes
        ----------
        MJC_:  numpy.ndarray (1,)
               array containing the Multivariate Local Join Counts ...
        """

        self.connectivity = connectivity

    def fit(self, variables):
        """
        Arguments
        ---------
        y       :   numpy.ndarray
                    array containing binary (0/1) data
        Returns
        -------
        the fitted estimator.
        Notes
        -----
        Technical details and derivations can be found in :cite:`AnselinLi2019`.
        """

        # Need not be flattened?

        w = self.connectivity
        w.transformation = 'b'

        self.MJC_ = self._statistic(variables, w)

        return self

    @staticmethod
    def _statistic(variables, w):
        # Create adjacency list. Note that remove_symmetric=False -
        # different from the esda.Join_Counts() function.
        adj_list = w.to_adjlist(remove_symmetric=False)

        # The zseries
        zseries = [pd.Series(i, index=w.id_order) for i in variables]
        # The focal values
        focal = [zseries[i].loc[adj_list.focal].values for
                 i in range(len(variables))]
        # The neighbor values
        neighbor = [zseries[i].loc[adj_list.neighbor].values for
                    i in range(len(variables))]

        # Find instances where all surrounding 
        # focal and neighbor values == 1
        focal_all = np.array(np.all(np.dstack(focal)==1, 
                                    axis=2))
        neighbor_all = np.array(np.all(np.dstack(neighbor)==1, 
                                       axis=2))
        MCLC = (focal_all == True) & (neighbor_all == True)
        # Convert list of True/False to boolean array 
        # and unlist (necessary for building pd.DF)
        MCLC = list(MCLC*1)
        
        # Create a df that uses the adjacency list
        # focal values and the BBs counts
        adj_list_MCLC = pd.DataFrame(adj_list.focal.values,
                                     MCLC).reset_index()
        # Temporarily rename the columns
        adj_list_MCLC.columns = ['MCLC', 'ID']
        adj_list_MCLC = adj_list_MCLC.groupby(by='ID').sum()

        return (adj_list_MCLC.MCLC.values)