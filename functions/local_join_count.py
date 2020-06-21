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
        BB:  numpy.ndarray (1,)
             array containing the estimated Local Join Count coefficients,
             where element [0,0] is the number of Local Join Counts, ...
        """

        self.connectivity = connectivity
        self.permutations = permutations

    def fit(self, y, permutations=999):
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
        
        self.y = y
        self.n = len(y)
        self.w = w
        
        self.BB = self._statistic(y, w)
        
        if permutations:
            self._crand()
            sim = np.transpose(self.rjoins)
            above = sim >= self.BB
            larger = above.sum(0)
            low_extreme = (self.permutations - larger) < larger
            larger[low_extreme] = self.permutations - larger[low_extreme]
            # 1 - simulated p-value? or just the simulated p-value?
            # values of 0.001 seem to be NA or error?
            self.p_sim = (larger + 1.0) / (permutations + 1.0)

        # Need the >>> return self to get the associated .BB attribute
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
    
    def _crand(self):
        """
        conditional randomization

        for observation i with ni neighbors,  the candidate set cannot include
        i (we don't want i being a neighbor of i). we have to sample without
        replacement from a set of ids that doesn't include i. numpy doesn't
        directly support sampling wo replacement and it is expensive to
        implement this. instead we omit i from the original ids,  permute the
        ids and take the first ni elements of the permuted ids as the
        neighbors to i in each randomization.

        """
        # converted z to y
        # renamed lisas to joins
        y = self.y
        n = len(y)
        joins = np.zeros((self.n, self.permutations))
        n_1 = self.n - 1
        prange = list(range(self.permutations))
        k = self.w.max_neighbors + 1
        nn = self.n - 1
        rids = np.array([np.random.permutation(nn)[0:k] for i in prange])
        ids = np.arange(self.w.n)
        ido = self.w.id_order
        w = [self.w.weights[ido[i]] for i in ids]
        wc = [self.w.cardinalities[ido[i]] for i in ids]

        for i in range(self.w.n):
            idsi = ids[ids != i]
            np.random.shuffle(idsi)
            tmp = y[idsi[rids[:, 0:wc[i]]]]
            joins[i] = y[i] * (w[i] * tmp).sum(1)
        self.rjoins = joins