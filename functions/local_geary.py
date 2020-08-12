import numpy as np
import pandas as pd
import warnings
from scipy import sparse
from scipy import stats
from sklearn.base import BaseEstimator
import libpysal as lp
from esda.crand import (
    crand as _crand_plus,
    njit as _njit,
    _prepare_univariate
)


PERMUTATIONS = 999
SIG = 0.05


class Local_Geary(BaseEstimator):
    """Local Geary - Univariate"""

    def __init__(self, connectivity=None, autocorr=False, sig=SIG,
                 permutations=PERMUTATIONS, n_jobs=1, keep_simulations=True,
                 seed=None):
        """
        connectivity     : scipy.sparse matrix object
                           the connectivity structure describing
                           the relationships between observed units.
                           Need not be row-standardized.
        autocorr         : boolean
                           (default=False)
                           If True use, indicate if an observation 
                           belongs to an outlier, cluster, undefined,
                           or non-significant group. 1 = outlier, 
                           2 = cluster, 3 = undefined, 4 = non-significant.
        sig              : float
                           (default=0.05)
                           Default significance threshold used for 
                           creation of autocorr groups.
        permutations     : int
                           number of random permutations for calculation
                           of pseudo p_values
        n_jobs           : int
                           Number of cores to be used in the conditional
                           randomisation. If -1, all available cores are used.
        keep_simulations : Boolean
                           (default=True)
                           If True, the entire matrix of replications under
                           the null is stored in memory and accessible;
                           otherwise, replications are not saved
        seed             : None/int
                           Seed to ensure reproducibility of conditional
                           randomizations. Must be set here, and not outside
                           of the function, since numba does not correctly
                           interpret external seeds nor
                           numpy.random.RandomState instances.

        Attributes
        ----------
        localG          : numpy array
                          array containing the observed univariate
                          Local Geary values.
        p_sim           : numpy array
                          array containing the simulated
                          p-values for each unit.
        """

        self.connectivity = connectivity
        self.autocorr = autocorr
        self.sig = sig
        self.permutations = permutations
        self.n_jobs = n_jobs
        self.keep_simulations = keep_simulations
        self.seed = seed

    def fit(self, x, n_jobs=1, permutations=999):
        """
        Arguments
        ---------
        x                : numpy.ndarray
                           array containing continuous data

        Returns
        -------
        the fitted estimator.

        Notes
        -----
        Technical details and derivations can be found in :cite:`Anselin1995`.

        Examples
        --------
        Guerry data replication GeoDa tutorial
        >>> import libpysal as lp
        >>> import geopandas as gpd
        >>> guerry = lp.examples.load_example('Guerry')
        >>> guerry_ds = gpd.read_file(guerry.get_path('Guerry.shp'))
        >>> w = libpysal.weights.Queen.from_dataframe(guerry_ds)
        >>> y = guerry_ds['Donatns']
        >>> lG = Local_Geary(connectivity=w).fit(y)
        >>> lG.localG[0:5]
        >>> lG.p_sim[0:5]
        """
        x = np.asarray(x).flatten()

        w = self.connectivity
        w.transform = 'r'

        self.localG = self._statistic(x, w)

        if self.permutations:
            self.p_sim, self.rlocalG = _crand_plus(
                z=(x - np.mean(x))/np.std(x),
                w=w,
                observed=self.localG,
                permutations=permutations,
                keep=True,
                n_jobs=n_jobs,
                stat_func=_local_geary
            )
            
        if self.autocorr:    
            Eij_mean = np.mean(self.localG)
            x_mean = np.mean(x)
            # Create empty vector where default is undefined.
            self.q = np.ones(len(x))*3
            # Outliers
            self.q[(self.localG < Eij_mean) & \
                   (y > y_mean) & \
                   (self.p_sim<=self.sig)] = 1
            # Clusters
            self.q[(self.localG < Eij_mean) & \
                   (y < y_mean) & \
                   (self.p_sim<=self.sig)] = 2
            # Undefined: default value in empty vector
            # Non-significant
            self.q[self.p_sim > self.sig] = 4

        del (self.keep_simulations, self.n_jobs,
             self.permutations, self.seed, self.rlocalG,
             self.connectivity)

        return self

    @staticmethod
    def _statistic(x, w):
        # Caclulate z-scores for x
        zscore_x = (x - np.mean(x))/np.std(x)
        # Create focal (xi) and neighbor (zi) values
        adj_list = w.to_adjlist(remove_symmetric=False)
        zseries = pd.Series(zscore_x, index=w.id_order)
        zi = zseries.loc[adj_list.focal].values
        zj = zseries.loc[adj_list.neighbor].values
        # Carry out local Geary calculation
        gs = sum(list(w.weights.values()), []) * (zi-zj)**2
        # Reorganize data
        adj_list_gs = pd.DataFrame(adj_list.focal.values, gs).reset_index()
        adj_list_gs.columns = ['gs', 'ID']
        adj_list_gs = adj_list_gs.groupby(by='ID').sum()

        localG = adj_list_gs.gs.values

        return (localG)

# --------------------------------------------------------------
# Conditional Randomization Function Implementations
# --------------------------------------------------------------

# Note: does not using the scaling parameter


@_njit(fastmath=True)
def _local_geary(i, z, permuted_ids, weights_i, scaling):
    zi, zrand = _prepare_univariate(i, z, permuted_ids, weights_i)
    return (zi-zrand)**2 @ weights_i