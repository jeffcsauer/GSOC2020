import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from libpysal import weights


PERMUTATIONS = 999


class Local_Join_Count(BaseEstimator):

    """Univariate Local Join Count Statistic"""

    def __init__(self, connectivity=None, permutations=PERMUTATIONS):
        """
        Initialize a Local_Join_Count estimator
        Arguments
        ---------
        connectivity     : scipy.sparse matrix object
                           the connectivity structure describing
                           the relationships between observed units.
                           Need not be row-standardized.
        Attributes
        ----------
        LJC             : numpy array
                          array containing the univariate
                          Local Join Count (LJC).
        p_sim           : numpy array
                          array containing the simulated
                          p-values for each unit.

        """

        self.connectivity = connectivity
        self.permutations = permutations

    def fit(self, y, permutations=999):
        """
        Arguments
        ---------
        y               : numpy.ndarray
                          array containing binary (0/1) data
        Returns
        -------
        the fitted estimator.

        Notes
        -----
        Technical details and derivations found in :cite:`AnselinLi2019`.

        Examples
        --------
        >>> import libpysal
        >>> w = libpysal.weights.lat2W(4, 4)
        >>> y = np.ones(16)
        >>> y[0:8] = 0
        >>> LJC_uni = Local_Join_Count(connectivity=w).fit(y)
        >>> LJC_uni.LJC
        >>> LJC_uni.p_sim

        Guerry data replicating GeoDa tutorial
        >>> import libpysal
        >>> import geopandas as gpd
        >>> guerry = libpysal.examples.load_example('Guerry')
        >>> guerry_ds = gpd.read_file(guerry.get_path('Guerry.shp'))
        >>> guerry_ds['SELECTED'] = 0
        >>> guerry_ds.loc[(guerry_ds['Donatns'] > 10997), 'SELECTED'] = 1
        >>> w = libpysal.weights.Queen.from_dataframe(guerry_ds)
        >>> LJC_uni = Local_Join_Count(connectivity=w).fit(guerry_ds['SELECTED'])
        >>> LJC_uni.LJC
        >>> LJC_uni.p_sim
        """
        y = np.asarray(y).flatten()

        w = self.connectivity
        # Fill the diagonal with 0s
        w = weights.util.fill_diagonal(w, val=0)
        w.transform = 'b'

        self.y = y
        self.n = len(y)
        self.w = w

        self.LJC = self._statistic(y, w)

        if permutations:
            self._crand()
            sim = np.transpose(self.rjoins)
            above = sim >= self.LJC
            larger = above.sum(0)
            low_extreme = (self.permutations - larger) < larger
            larger[low_extreme] = self.permutations - larger[low_extreme]
            self.p_sim = (larger + 1.0) / (permutations + 1.0)
            # Set p-values for those with LJC of 0 to NaN
            self.p_sim[self.LJC == 0] = 'NaN'

        return self

    @staticmethod
    def _statistic(y, w):
        # Create adjacency list. Note that remove_symmetric=False - this is
        # different from the esda.Join_Counts() function.
        adj_list = w.to_adjlist(remove_symmetric=False)
        zseries = pd.Series(y, index=w.id_order)
        focal = zseries.loc[adj_list.focal].values
        neighbor = zseries.loc[adj_list.neighbor].values
        LJC = (focal == 1) & (neighbor == 1)
        adj_list_LJC = pd.DataFrame(adj_list.focal.values,
                                    LJC.astype('uint8')).reset_index()
        adj_list_LJC.columns = ['LJC', 'ID']
        adj_list_LJC = adj_list_LJC.groupby(by='ID').sum()
        LJC = adj_list_LJC.LJC.values
        return (LJC)

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