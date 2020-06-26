import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.base import BaseEstimator
from libpysal import weights

PERMUTATIONS = 999

class Local_Join_Count_MV(BaseEstimator):

    """Multivariate local join counts"""

    def __init__(self, connectivity=None, permutations=PERMUTATIONS):
        """
        Initialize a Local_Join_Count_MV estimator
        Arguments
        ---------
        connectivity:   scipy.sparse matrix object
                        the connectivity structure describing the relationships
                        between observed units. Will be row-standardized.
        Attributes
        ----------
        LJC       :   numpy.ndarray
                      array containing the estimated Multivariate Local Join Counts.
        p_sim       :   numpy.ndarray
                        array containing the simulated p-values for each unit.
        """

        self.connectivity = connectivity
        self.permutations = permutations

    def fit(self, variables, permutations=999):
        """
        Arguments
        ---------
        variables     :   numpy.ndarray
                          array(s) containing binary (0/1) data
        Returns
        -------
        the fitted estimator.
        Notes
        -----
        Technical details and derivations can be found in :cite:`AnselinLi2019`.
        """

        w = self.connectivity
        # Fill the diagonal with 0s
        w = weights.util.fill_diagonal(w, val=0)
        w.transform = 'b'
        
        self.n = len(variables[0])
        self.w = w
        
        self.variables = variables
        
        self.ext = np.prod(np.vstack(variables), axis=0)

        self.LJC = self._statistic(variables, w)
        
        if permutations:
            self._crand()
            sim = np.transpose(self.rjoins)
            above = sim >= self.LJC
            larger = above.sum(0)
            low_extreme = (self.permutations - larger) < larger
            larger[low_extreme] = self.permutations - larger[low_extreme]
            self.p_sim = (larger + 1.0) / (permutations + 1.0)
            # Set p-values for those with LJC of 0 to NaN
            self.p_sim[self.LJC==0] = 'NaN'

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
        # converted y to z
        # renamed lisas to joins
        ext = self.ext
        # Get length based on first variable
        n = len(ext)
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
            # Mirroring moran_local_bv()
            tmp = ext[idsi[rids[:, 0:wc[i]]]]
            joins[i] = ext[i] * (w[i] * tmp).sum(1)
        self.rjoins = joins