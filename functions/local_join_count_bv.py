import numpy as np
import pandas as pd
import warnings
from scipy import sparse
from sklearn.base import BaseEstimator
from libpysal import weights


PERMUTATIONS = 999


class Local_Join_Count_BV(BaseEstimator):

    """Univariate Local Join Count Statistic"""

    def __init__(self, connectivity=None, permutations=PERMUTATIONS):
        """
        Initialize a Local_Join_Count_BV estimator
        Arguments
        ---------
        connectivity     : scipy.sparse matrix object
                           the connectivity structure describing
                           the relationships between observed units.
                           Need not be row-standardized.
        Attributes
        ----------
        LJC              : numpy.ndarray
                           array containing the estimated
                           Bivariate Local Join Counts
        p_sim            : numpy.ndarray
                           array containing the simulated
                           p-values for each unit.
        """

        self.connectivity = connectivity
        self.permutations = permutations

    def fit(self, x, z, case="CLC", permutations=999):
        """
        Arguments
        ---------
        x                : numpy.ndarray
                           array containing binary (0/1) data
        z                : numpy.ndarray
                           array containing binary (0/1) data
        Returns
        -------
        the fitted estimator.

        Notes
        -----
        Technical details and derivations can be found in :cite:`AnselinLi2019`.

        Examples
        --------
        >>> import libpysal
        >>> w = libpysal.weights.lat2W(4, 4)
        >>> x = np.ones(16)
        >>> x[0:8] = 0
        >>> z = [0,1,0,1,1,1,1,1,0,0,1,1,0,0,1,1]
        >>> LJC_BV_C1 = Local_Join_Count_BV(connectivity=w).fit(x, z, case="BJC")
        >>> LJC_BV_C2 = Local_Join_Count_BV(connectivity=w).fit(x, z, case="CLC")
        >>> LJC_BV_C1.LJC
        >>> LJC_BV_C1.p_sim
        >>> LJC_BV_C2.LJC
        >>> LJC_BV_C2.p_sim

        Commpop data replicating GeoDa tutorial (Case 1)
        >>> import libpysal
        >>> import geopandas as gpd
        >>> commpop = gpd.read_file("https://github.com/jeffcsauer/GSOC2020/raw/master/validation/data/commpop.gpkg")
        >>> w = libpysal.weights.Queen.from_dataframe(commpop)
        >>> LJC_BV_Case1 = Local_Join_Count_BV(connectivity=w).fit(commpop['popneg'], commpop['popplus'], case='BJC')
        >>> LJC_BV_Case1.LJC
        >>> LJC_BV_Case1.p_sim

        Guerry data replicating GeoDa tutorial (Case 2)
        >>> import geopandas as gpd
        >>> import libpysal
        >>> guerry = gpd.read_file('https://github.com/jeffcsauer/GSOC2020/raw/master/validation/data/guerry/guerry_geodavalues.gpkg')
        >>> guerry['infq5'] = 0
        >>> guerry['donq5'] = 0
        >>> guerry.loc[(guerry['Infants'] > 23574), 'infq5'] = 1
        >>> guerry.loc[(guerry['Donatns'] > 10973), 'donq5'] = 1
        >>> w = libpysal.weights.Queen.from_dataframe(guerry)
        >>> LJC_BV_Case2 = Local_Join_Count_BV(connectivity=w).fit(guerry['infq5'], guerry['donq5'], case='CLC')
        >>> LJC_BV_Case2.LJC
        >>> LJC_BV_Case2.p_sim

        """
        x = np.asarray(x).flatten()
        z = np.asarray(z).flatten()

        w = self.connectivity
        # Fill the diagonal with 0s
        w = weights.util.fill_diagonal(w, val=0)
        w.transform = 'b'

        self.x = x
        self.z = z
        self.n = len(x)
        self.w = w
        self.case = case

        self.LJC = self._statistic(x, z, w, case=case)

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
    def _statistic(x, z, w, case):
        # Create adjacency list. Note that remove_symmetric=False - this is
        # different from the esda.Join_Counts() function.
        adj_list = w.to_adjlist(remove_symmetric=False)

        # First, set up a series that maps the values
        # to the weights table
        zseries_x = pd.Series(x, index=w.id_order)
        zseries_z = pd.Series(z, index=w.id_order)

        # Map the values to the focal (i) values
        focal_x = zseries_x.loc[adj_list.focal].values
        focal_z = zseries_z.loc[adj_list.focal].values

        # Map the values to the neighbor (j) values
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
            raise NotImplementedError(f'The requested LJC method ({case}) \
            is not currently supported!')

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
        x = self.x
        z = self.z
        case = self.case

        n = len(x)
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
            tmp_x = x[idsi[rids[:, 0:wc[i]]]]
            tmp_z = z[idsi[rids[:, 0:wc[i]]]]
            if case == "BJC":
                joins[i] = x[i] * (w[i] * tmp_z).sum(1)
            elif case == "CLC":
                joins[i] = z[i] * (w[i] * tmp_z * tmp_x).sum(1)
            else:
                raise NotImplementedError(f'The requested LJC method \
                ({case}) is not currently supported!')
        self.rjoins = joins