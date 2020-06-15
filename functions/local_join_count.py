import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
import libpysal

PERMUTATIONS = 999

class Local_Join_Count(BaseEstimator):

    """Local Join Count Statistic"""

    def __init__(self, connectivity=None, permutations=PERMUTATIONS):
        """
        Initialize a Join_Counts_Local estimator
        Arguments
        ---------
        connectivity:   scipy.sparse matrix object
                        the connectivity structure describing the relationships
                        between observed units. Need not be row-standardized.
        Attributes
        ----------
        BB_:  numpy.ndarray (1,)
              array containing the estimated Local Join Count coefficients,
              where element [0,0] is the number of Local Join Counts, ...
        """

        self.connectivity = connectivity
        self.permutations = permutations

    def fit(self, y, permutations):
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
        Technical details and derivations found in :cite:`AnselinLi2019`.
        """
        y = np.asarray(y).flatten()

        w = self.connectivity
        # Binary weights are needed for this statistic
        w.transformation = 'b'

        self.BB_ = self._statistic(y, w)
        
        if permutations:
            print(PERMUTATIONS)

        # Need the >>> return self to get the associated .BB_ attribute
        # (significance in future, i.e. self.reference_distribution_ in lee.py)
        return self

    @staticmethod
    def _statistic(y, w):
        # Create adjacency list. Note that remove_symmetric=False - this is
        # different from the esda.Join_Counts() function.
        adj_list = w.to_adjlist(remove_symmetric=False)
        zseries = pd.Series(y, index=w.id_order)
        focal = zseries.loc[adj_list.focal].values
        neighbor = zseries.loc[adj_list.neighbor].values
        BB = (focal == 1) & (neighbor == 1)
        adj_list_BB = pd.DataFrame(adj_list.focal.values,
                                   BB.astype('uint8')).reset_index()
        adj_list_BB.columns = ['BB', 'ID']
        adj_list_BB = adj_list_BB.groupby(by='ID').sum()
        BB = adj_list_BB.BB.values
        return (BB)