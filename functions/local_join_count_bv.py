import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.base import BaseEstimator

PERMUTATIONS = 999

class Local_Join_Count_BV(BaseEstimator):

    """Global Spatial Pearson Statistic"""

    def __init__(self, connectivity=None, permutations=PERMUTATIONS):
        """
        Initialize a Join_Counts_Local estimator
        Arguments
        ---------
        connectivity:   scipy.sparse matrix object
                        the connectivity structure describing the relationships
                        between observed units. Will be row-standardized.
        Attributes
        ----------
        LJC :  numpy.ndarray (1,)
               array containing the estimated Bivariate Local Join Counts
        """

        self.connectivity = connectivity
        self.permutations = permutations

    def fit(self, x, z, case=None, permutations=999):
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
        x = np.asarray(x).flatten()
        z = np.asarray(z).flatten()

        w = self.connectivity
        w.transformation = 'b'

        self.x = x
        self.z = z
        self.n = len(x)
        self.w = w

        self.LJC = self._statistic(x, z, w, case=case)

        if permutations:
            self._crand()
            sim = np.transpose(self.rjoins)
            above = sim >= self.LJC
            larger = above.sum(0)
            low_extreme = (self.permutations - larger) < larger
            larger[low_extreme] = self.permutations - larger[low_extreme]
            # 1 - simulated p-value? or just the simulated p-value?
            # values of 0.001 seem to be NA or error?
            self.p_sim = (larger + 1.0) / (permutations + 1.0)

        return self

    @staticmethod
    def _statistic(x, z, w, case=None):
        # Create adjacency list. Note that remove_symmetric=False - this is
        # different from the esda.Join_Counts() function.
        adj_list = w.to_adjlist(remove_symmetric=False)

        # First, set up a series that maps the y values
        # (input as self.y) to the weights table
        zseries_x = pd.Series(x, index=w.id_order)
        zseries_z = pd.Series(z, index=w.id_order)

        # Next, map the y values to the focal (i) values
        focal_x = zseries_x.loc[adj_list.focal].values
        focal_z = zseries_z.loc[adj_list.focal].values

        # Repeat the mapping but for the neighbor (j) values
        neighbor_x = zseries_x.loc[adj_list.neighbor].values
        neighbor_z = zseries_z.loc[adj_list.neighbor].values

        if case == "BJC":
            BJC = (focal_x == 1) & (focal_z == 0) & \
                  (neighbor_x == 0) & (neighbor_z == 1)
            adj_list_BJC = pd.DataFrame(adj_list.focal.values,
                                        BJC.astype('uint8')).reset_index()
            adj_list_BJC.columns = ['BJC', 'ID']
            adj_list_BJC = adj_list_BJC.groupby(by='ID').sum()
            return (adj_list_BJC.BJC.values)
        elif case == "CLC":
            CLC = (focal_x == 1) & (focal_z == 1) & \
                  (neighbor_x == 1) & (neighbor_z == 1)
            adj_list_CLC = pd.DataFrame(adj_list.focal.values,
                                        CLC.astype('uint8')).reset_index()
            adj_list_CLC.columns = ['CLC', 'ID']
            adj_list_CLC = adj_list_CLC.groupby(by='ID').sum()
            return (adj_list_CLC.CLC.values)
        else:
            print("Please specify which type of bivariate Local Join Count \
            you would like to calculate (either 'BJC' or 'CLC'). See Anselin \
            and Li 2019 p. 9-10 for more information")
            
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
        y = self.x
        z = self.z
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
            # Add z[i]?
            joins[i] = y[i] * z[i] * (w[i] * tmp).sum(1)
        self.rjoins = joins