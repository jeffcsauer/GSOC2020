# Introduction

*A huge thank you to to @sjsrey, @ljwolf, @darribas, @slumnitz, @TaylorOshan, and the larger [PySAL group](https://github.com/orgs/pysal/people) for acting as mentors throughout the summer and welcoming me to the PySAL community. This was a wonderful exercise in developing open-source software and I would highly recommend future students consider contributing to PySAL!*

Welcome to the summary document for the Google Summer of Code 2020 project entitled [*PySAL ESDA Enhancements: Local join count and LOSH statistics*](https://docs.google.com/document/d/1WjHjy5Eyk4WG5QWfnsnhWg1r4-e09JXXCx0iaPphg6c/edit?usp=sharing). In short, the original objectives of the project were as follows:

- Implement Local Spatial Heteroskedasticity (LOSH) estimators
- Implement univariate, bivariate, and multivariate Local Join Count (LJC) estimators

Each estimator includes docstrings, doctests, tests, and an example notebook demonstrating its application. After completing the above objectives towards the end of July, additional objectives were outlined. These additional objectives were to be completed as much as possible by the end of Google Summer of Code 2020 and serve as a bridge to continue contributing to PySAL in the future. 

- Implement univariate and multivariate local Geary estimators
- Experiment with conditional randomization of the LOSH estimator

As of the writing of this document and the end of GSOC, the following progress has been made across all of the above objectives:

| Function              | Generating correct values | Generating correct inference | Documentation | Pull request | Overall, complete?      |
|-----------------------|---------------------------|------------------------------|---------------|--------------|-------------------------|
| `LOSH`                | Yes                       | Yes                          | Yes           | [PR#139](https://github.com/pysal/esda/pull/139)     | Yes |
| `Local_Join_Count` <br> (univariate)      | Yes        | Yes                          | Yes           | [PR#139](https://github.com/pysal/esda/pull/139)     | Yes |
| `Local_Join_Count_BV` <br> (bivariate)    | Yes        | Yes                          | Yes           | [PR#139](https://github.com/pysal/esda/pull/139)     | Yes |
| `Local_Join_Count_MV` <br> (multivariate) | Yes        | Yes                          | Yes           | [PR#139](https://github.com/pysal/esda/pull/139)     | Yes |
| `Local_Geary` <br> (univariate)           | Yes        | Yes                          | Yes   | TBD | Yes |
| `Local_Geary_MV` <br> (multivariate)      | Yes        | No                           | In progress   | TBD | No |

For ease of access and perpetuity, each of the above functions are copied below. Additional links to documentation and where the most recent version of the estimator may be found is provided with each function. If you are interested in viewing the development history of these estimators, I recommend you visit the [Github repository](https://github.com/jeffcsauer/GSOC2020) where the majority of the work was carried out. 

## Local Spatial Heteroskedasticity (LOSH)

Links to:
- [Documentation](https://github.com/jeffcsauer/GSOC2020/blob/master/docs/LOSH.ipynb)
- [Location in PySAL `esda`](https://github.com/pysal/esda/blob/master/esda/losh.py)
- [Location in GSOC 2020 workbook](https://github.com/jeffcsauer/GSOC2020/blob/master/functions/losh.py)

Stable release version:

<details>
  <summary>Click to expand code</summary>
  
``` python

import numpy as np
import warnings
from scipy import sparse
from scipy import stats
from sklearn.base import BaseEstimator
import libpysal as lp


class losh(BaseEstimator):
    """Local spatial heteroscedasticity (LOSH)"""

    def __init__(self, connectivity=None, inference=None):
        """
        Initialize a losh estimator

        Arguments
        ---------
        connectivity     : scipy.sparse matrix object
                           the connectivity structure describing the
                           relationships between observed units.
        inference        : str
                           describes type of inference to be used. options are
                           "chi-square" or "permutation" methods.

        Attributes
        ----------
        Hi               : numpy array
                           Array of LOSH values for each spatial unit.
        ylag             : numpy array
                           Spatially lagged y values.
        yresid           : numpy array
                           Spatially lagged residual values.
        VarHi            : numpy array
                           Variance of Hi.
        pval             : numpy array
                           P-values for inference based on either
                           "chi-square" or "permutation" methods.
        """

        self.connectivity = connectivity
        self.inference = inference

    def fit(self, y, a=2):
        """
        Arguments
        ---------
        y                : numpy.ndarray
                           array containing continuous data
        a                : int
                           residual multiplier. Default is 2 in order
                           to generate a variance measure. Users may
                           use 1 for absolute deviations.

        Returns
        -------
        the fitted estimator.

        Notes
        -----
        Technical details and derivations can be found in :cite:`OrdGetis2012`.

        Examples
        --------
        >>> import libpysal
        >>> w = libpysal.io.open(libpysal.examples.get_path("stl.gal")).read()
        >>> f = libpysal.io.open(libpysal.examples.get_path("stl_hom.txt"))
        >>> y = np.array(f.by_col['HR8893'])
        >>> from esda import losh
        >>> ls = losh(connectivity=w, inference="chi-square").fit(y)
        >>> np.round(ls.Hi[0], 3)
        >>> np.round(ls.pval[0], 3)

        Boston housing data replicating R spdep::LOSH()
        >>> import libpysal
        >>> import geopandas as gpd
        >>> boston = libpysal.examples.load_example('Bostonhsg')
        >>> boston_ds = gpd.read_file(boston.get_path('boston.shp'))
        >>> w = libpysal.weights.Queen.from_dataframe(boston_ds)
        >>> ls = losh(connectivity=w, inference="chi-square").fit(boston['NOX'])
        >>> np.round(ls.Hi[0], 3)
        >>> np.round(ls.VarHi[0], 3)
        """
        y = np.asarray(y).flatten()

        w = self.connectivity

        self.Hi, self.ylag, self.yresid, self.VarHi = self._statistic(y, w, a)

        if self.inference is None:
            return self
        elif self.inference == 'chi-square':
            if a != 2:
                warnings.warn(f'Chi-square inference assumes that a=2, but \
                a={a}. This means the inference will be invalid!')
            else:
                dof = 2/self.VarHi
                Zi = (2*self.Hi)/self.VarHi
                self.pval = 1 - stats.chi2.cdf(Zi, dof)
        else:
            raise NotImplementedError(f'The requested inference method \
            ({self.inference}) is not currently supported!')

        return self

    @staticmethod
    def _statistic(y, w, a):
        # Define what type of variance to use
        if a is None:
            a = 2
        else:
            a = a

        rowsum = np.array(w.sparse.sum(axis=1)).flatten()

        # Calculate spatial mean
        ylag = lp.weights.lag_spatial(w, y)/rowsum
        # Calculate and adjust residuals based on multiplier
        yresid = abs(y-ylag)**a
        # Calculate denominator of Hi equation
        denom = np.mean(yresid) * np.array(rowsum)
        # Carry out final Hi calculation
        Hi = lp.weights.lag_spatial(w, yresid) / denom
        # Calculate average of residuals
        yresid_mean = np.mean(yresid)
        # Calculate VarHi
        n = len(y)
        squared_rowsum = np.asarray(w.sparse.multiply(w.sparse).sum(axis=1)).flatten()

        VarHi = ((n-1)**-1) * \
                (denom**-2) * \
                ((np.sum(yresid**2)/n) - yresid_mean**2) * \
                ((n*squared_rowsum) - (rowsum**2))

        return (Hi, ylag, yresid, VarHi)

```
</details>
  
## Local Join Count (univariate)

Links to:
- [Documentation](https://github.com/jeffcsauer/GSOC2020/blob/master/docs/localjoincounts.ipynb)
- [Location in PySAL `esda`](https://github.com/pysal/esda/blob/master/esda/local_join_count.py)
- [Location in GSOC 2020 workbook](https://github.com/jeffcsauer/GSOC2020/blob/master/functions/local_join_count.py)

Stable release version:

<details>
  <summary>Click to expand code</summary>
  
``` python
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from libpysal import weights
from esda.crand import (
    crand as _crand_plus,
    njit as _njit,
    _prepare_univariate
)


PERMUTATIONS = 999


class Local_Join_Count(BaseEstimator):

    """Univariate Local Join Count Statistic"""

    def __init__(self, connectivity=None, permutations=PERMUTATIONS, n_jobs=1, 
                 keep_simulations=True, seed=None):
        """
        Initialize a Local_Join_Count estimator
        Arguments
        ---------
        connectivity     : scipy.sparse matrix object
                           the connectivity structure describing
                           the relationships between observed units.
                           Need not be row-standardized.
        permutations     : int
                           number of random permutations for calculation of pseudo
                           p_values
        n_jobs           : int
                           Number of cores to be used in the conditional randomisation. If -1,
                           all available cores are used.    
        keep_simulations : Boolean
                           (default=True)
                           If True, the entire matrix of replications under the null 
                           is stored in memory and accessible; otherwise, replications 
                           are not saved
        seed             : None/int
                           Seed to ensure reproducibility of conditional randomizations. 
                           Must be set here, and not outside of the function, since numba 
                           does not correctly interpret external seeds 
                           nor numpy.random.RandomState instances.              
                           
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
        self.n_jobs = n_jobs
        self.keep_simulations = keep_simulations
        self.seed = seed

    def fit(self, y, n_jobs=1, permutations=999):
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
        # Need to ensure that the np.array() are of
        # dtype='float' for numba
        y = np.array(y, dtype='float')

        w = self.connectivity
        # Fill the diagonal with 0s
        w = weights.util.fill_diagonal(w, val=0)
        w.transform = 'b'
        
        keep_simulations = self.keep_simulations
        n_jobs = self.n_jobs
        seed = self.seed
        
        self.y = y
        self.n = len(y)
        self.w = w

        self.LJC = self._statistic(y, w)
        
        if permutations:
            self.p_sim, self.rjoins = _crand_plus(
                z=self.y, 
                w=self.w, 
                observed=self.LJC,
                permutations=permutations, 
                keep=keep_simulations, 
                n_jobs=n_jobs,
                stat_func=_ljc_uni
            )
            # Set p-values for those with LJC of 0 to NaN
            self.p_sim[self.LJC == 0] = 'NaN'
        
        del (self.n, self.keep_simulations, self.n_jobs, 
             self.permutations, self.seed, self.w, self.y,
             self.connectivity, self.rjoins)
        
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
        LJC = np.array(adj_list_LJC.LJC.values, dtype='float')
        return (LJC)

# --------------------------------------------------------------
# Conditional Randomization Function Implementations
# --------------------------------------------------------------

# Note: scaling not used

@_njit(fastmath=True)
def _ljc_uni(i, z, permuted_ids, weights_i, scaling):
    zi, zrand = _prepare_univariate(i, z, permuted_ids, weights_i)
    return zi * (zrand @ weights_i)
```

</details>

## Local Join Count (bivariate)

Links to:
- [Documentation](https://github.com/jeffcsauer/GSOC2020/blob/master/docs/localjoincounts.ipynb)
- [Location in PySAL `esda`](https://github.com/pysal/esda/blob/master/esda/local_join_count_bv.py)
- [Location in GSOC 2020 workbook](https://github.com/jeffcsauer/GSOC2020/blob/master/functions/local_join_count_bv.py)

Stable release version:

<details>
  <summary>Click to expand code</summary>
  
``` python
import numpy as np
import pandas as pd
import warnings
from scipy import sparse
from sklearn.base import BaseEstimator
from libpysal import weights
from esda.crand import (
    crand as _crand_plus,
    njit as _njit,
    _prepare_univariate,
    _prepare_bivariate
)


PERMUTATIONS = 999


class Local_Join_Count_BV(BaseEstimator):

    """Univariate Local Join Count Statistic"""

    def __init__(self, connectivity=None, permutations=PERMUTATIONS, n_jobs=1, 
                 keep_simulations=True, seed=None):
        """
        Initialize a Local_Join_Count_BV estimator
        Arguments
        ---------
        connectivity     : scipy.sparse matrix object
                           the connectivity structure describing
                           the relationships between observed units.
                           Need not be row-standardized.
        permutations     : int
                           number of random permutations for calculation of pseudo
                           p_values
        n_jobs           : int
                           Number of cores to be used in the conditional randomisation. If -1,
                           all available cores are used.    
        keep_simulations : Boolean
                           (default=True)
                           If True, the entire matrix of replications under the null 
                           is stored in memory and accessible; otherwise, replications 
                           are not saved
        seed             : None/int
                           Seed to ensure reproducibility of conditional randomizations. 
                           Must be set here, and not outside of the function, since numba 
                           does not correctly interpret external seeds 
                           nor numpy.random.RandomState instances.              
                           
        """

        self.connectivity = connectivity
        self.permutations = permutations
        self.n_jobs = n_jobs
        self.keep_simulations = keep_simulations
        self.seed = seed

    def fit(self, x, z, case="CLC", n_jobs=1, permutations=999):
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
        >>> import libpysal
        >>> import geopandas as gpd
        >>> guerry = libpysal.examples.load_example('Guerry')
        >>> guerry_ds = gpd.read_file(guerry.get_path('Guerry.shp'))
        >>> guerry_ds['infq5'] = 0
        >>> guerry_ds['donq5'] = 0
        >>> guerry_ds.loc[(guerry_ds['Infants'] > 23574), 'infq5'] = 1
        >>> guerry_ds.loc[(guerry_ds['Donatns'] > 10973), 'donq5'] = 1
        >>> w = libpysal.weights.Queen.from_dataframe(guerry_ds)
        >>> LJC_BV_Case2 = Local_Join_Count_BV(connectivity=w).fit(guerry_ds['infq5'], guerry_ds['donq5'], case='CLC')
        >>> LJC_BV_Case2.LJC
        >>> LJC_BV_Case2.p_sim
        """
        # Need to ensure that the np.array() are of
        # dtype='float' for numba
        x = np.array(x, dtype='float')
        z = np.array(z, dtype='float')

        w = self.connectivity
        # Fill the diagonal with 0s
        w = weights.util.fill_diagonal(w, val=0)
        w.transform = 'b'

        self.x = x
        self.z = z
        self.n = len(x)
        self.w = w
        self.case = case
        
        keep_simulations = self.keep_simulations
        n_jobs = self.n_jobs
        seed = self.seed

        self.LJC = self._statistic(x, z, w, case=case)

        if permutations:
            if case == "BJC":
                self.p_sim, self.rjoins = _crand_plus(
                    z=np.column_stack((x, z)),
                    w=self.w, 
                    observed=self.LJC,
                    permutations=permutations, 
                    keep=True, 
                    n_jobs=n_jobs,
                    stat_func=_ljc_bv_case1
                )
                # Set p-values for those with LJC of 0 to NaN
                self.p_sim[self.LJC == 0] = 'NaN'
            elif case == "CLC":
                self.p_sim, self.rjoins = _crand_plus(
                    z=np.column_stack((x, z)),
                    w=self.w, 
                    observed=self.LJC,
                    permutations=permutations, 
                    keep=True, 
                    n_jobs=n_jobs,
                    stat_func=_ljc_bv_case2
                )
                # Set p-values for those with LJC of 0 to NaN
                self.p_sim[self.LJC == 0] = 'NaN'
            else:
                raise NotImplementedError(f'The requested LJC method ({case}) \
                is not currently supported!')

        del (self.n, self.keep_simulations, self.n_jobs, 
             self.permutations, self.seed, self.w, self.x,
             self.z, self.connectivity, self.rjoins)
                
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
            return (np.array(adj_list_BJC.BJC.values, dtype='float'))
        elif case == "CLC":
            CLC = (focal_x == 1) & (focal_z == 1) & \
                  (neighbor_x == 1) & (neighbor_z == 1)
            adj_list_CLC = pd.DataFrame(adj_list.focal.values,
                                        CLC.astype('uint8')).reset_index()
            adj_list_CLC.columns = ['CLC', 'ID']
            adj_list_CLC = adj_list_CLC.groupby(by='ID').sum()
            return (np.array(adj_list_CLC.CLC.values, dtype='float'))
        else:
            raise NotImplementedError(f'The requested LJC method ({case}) \
            is not currently supported!')

# --------------------------------------------------------------
# Conditional Randomization Function Implementations
# --------------------------------------------------------------

# Note: scaling not used

@_njit(fastmath=True)
def _ljc_bv_case1(i, z, permuted_ids, weights_i, scaling):
    zx = z[:, 0]
    zy = z[:, 1]
    zyi, zyrand = _prepare_univariate(i, zy, permuted_ids, weights_i)
    return zx[i] * (zyrand @ weights_i)

@_njit(fastmath=True)
def _ljc_bv_case2(i, z, permuted_ids, weights_i, scaling):
    zx = z[:, 0]
    zy = z[:, 1]
    zxi, zxrand, zyi, zyrand = _prepare_bivariate(i, z, permuted_ids, weights_i)
    zf = zxrand * zyrand
    return zy[i] * (zf @ weights_i)
```

</details>

## Local Join Count (multivariate)

Links to:
- [Documentation](https://github.com/jeffcsauer/GSOC2020/blob/master/docs/localjoincounts.ipynb)
- [Location in PySAL `esda`](https://github.com/pysal/esda/blob/master/esda/local_join_count_mv.py)
- [Location in GSOC 2020 workbook](https://github.com/jeffcsauer/GSOC2020/blob/master/functions/local_join_count_mv.py)

Stable release version:

<details>
  <summary>Click to expand code</summary>
  
``` python

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.base import BaseEstimator
from libpysal import weights
from esda.crand import (
    crand as _crand_plus,
    njit as _njit,
    _prepare_univariate
)


PERMUTATIONS = 999


class Local_Join_Count_MV(BaseEstimator):

    """Multivariate Local Join Count Statistic"""

    def __init__(self, connectivity=None, permutations=PERMUTATIONS, n_jobs=1, 
                 keep_simulations=True, seed=None):
        """
        Initialize a Local_Join_Count_MV estimator
        Arguments
        ---------
        connectivity     : scipy.sparse matrix object
                           the connectivity structure describing
                           the relationships between observed units.
                           Need not be row-standardized.
        permutations     : int
                           number of random permutations for calculation of pseudo
                           p_values
        n_jobs           : int
                           Number of cores to be used in the conditional randomisation. If -1,
                           all available cores are used.    
        keep_simulations : Boolean
                           (default=True)
                           If True, the entire matrix of replications under the null 
                           is stored in memory and accessible; otherwise, replications 
                           are not saved
        seed             : None/int
                           Seed to ensure reproducibility of conditional randomizations. 
                           Must be set here, and not outside of the function, since numba 
                           does not correctly interpret external seeds 
                           nor numpy.random.RandomState instances.              
                           
        """

        self.connectivity = connectivity
        self.permutations = permutations
        self.n_jobs = n_jobs
        self.keep_simulations = keep_simulations
        self.seed = seed

    def fit(self, variables, n_jobs=1, permutations=999):
        """
        Arguments
        ---------
        variables     : numpy.ndarray
                        array(s) containing binary (0/1) data
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
        >>> y = [0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,1]
        >>> LJC_MV = Local_Join_Count_MV(connectivity=w).fit([x, y, z])
        >>> LJC_MV.LJC
        >>> LJC_MV.p_sim

        Guerry data extending GeoDa tutorial
        >>> import libpysal
        >>> import geopandas as gpd
        >>> guerry = libpysal.examples.load_example('Guerry')
        >>> guerry_ds = gpd.read_file(guerry.get_path('Guerry.shp'))
        >>> guerry_ds['infq5'] = 0
        >>> guerry_ds['donq5'] = 0
        >>> guerry_ds['suic5'] = 0
        >>> guerry_ds.loc[(guerry_ds['Infants'] > 23574), 'infq5'] = 1
        >>> guerry_ds.loc[(guerry_ds['Donatns'] > 10973), 'donq5'] = 1
        >>> guerry_ds.loc[(guerry_ds['Suicids'] > 55564), 'suic5'] = 1
        >>> w = libpysal.weights.Queen.from_dataframe(guerry_ds)
        >>> LJC_MV = Local_Join_Count_MV(connectivity=w).fit([guerry_ds['infq5'], guerry_ds['donq5'], guerry_ds['suic5']])
        >>> LJC_MV.LJC
        >>> LJC_MV.p_sim
        """

        w = self.connectivity
        # Fill the diagonal with 0s
        w = weights.util.fill_diagonal(w, val=0)
        w.transform = 'b'

        self.n = len(variables[0])
        self.w = w

        self.variables = np.array(variables, dtype='float')
        
        keep_simulations = self.keep_simulations
        n_jobs = self.n_jobs
        seed = self.seed

        # Need to ensure that the product is an 
        # np.array() of dtype='float' for numba
        self.ext = np.array(np.prod(np.vstack(variables), axis=0), 
                            dtype='float')

        self.LJC = self._statistic(variables, w)

        if permutations:
            self.p_sim, self.rjoins = _crand_plus(
                z=self.ext, 
                w=self.w, 
                observed=self.LJC,
                permutations=permutations, 
                keep=True, 
                n_jobs=n_jobs,
                stat_func=_ljc_mv
            )
            # Set p-values for those with LJC of 0 to NaN
            self.p_sim[self.LJC == 0] = 'NaN'
        
        del (self.n, self.keep_simulations, self.n_jobs, 
             self.permutations, self.seed, self.w, self.ext,
             self.variables, self.connectivity, self.rjoins)

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
        focal_all = np.array(np.all(np.dstack(focal) == 1,
                                    axis=2))
        neighbor_all = np.array(np.all(np.dstack(neighbor) == 1,
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

        return (np.array(adj_list_MCLC.MCLC.values, dtype='float'))

# --------------------------------------------------------------
# Conditional Randomization Function Implementations
# --------------------------------------------------------------

# Note: scaling not used

@_njit(fastmath=True)
def _ljc_mv(i, z, permuted_ids, weights_i, scaling):
    zi, zrand = _prepare_univariate(i, z, permuted_ids, weights_i)
    return zi * (zrand @ weights_i)

```

</details>

## Local Geary (univariate)

Links to:
- [Documentation](https://github.com/jeffcsauer/GSOC2020/blob/master/docs/localgeary.ipynb)
- [Location in PySAL `esda`](https://github.com/pysal/esda/blob/master/esda/local_geary.py)
- [Location in GSOC 2020 workbook](https://github.com/jeffcsauer/GSOC2020/blob/master/functions/local_geary.py)

Stable release version:

<details>
  <summary>Click to expand code</summary>
  
``` python
    
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

    def __init__(self, connectivity=None, labels=False, sig=SIG,
                 permutations=PERMUTATIONS, n_jobs=1, keep_simulations=True,
                 seed=None):
        """
        connectivity     : scipy.sparse matrix object
                           the connectivity structure describing
                           the relationships between observed units.
                           Need not be row-standardized.
        labels           : boolean
                           (default=False)
                           If True use, label if an observation
                           belongs to an outlier, cluster, other,
                           or non-significant group. 1 = outlier,
                           2 = cluster, 3 = other, 4 = non-significant.
                           Note that this is not the exact same as the
                           cluster map produced by GeoDa.
        sig              : float
                           (default=0.05)
                           Default significance threshold used for
                           creation of labels groups.
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
        labs            : numpy array
                          array containing the labels for if each observation.
        """

        self.connectivity = connectivity
        self.labels = labels
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

        if self.labels:
            Eij_mean = np.mean(self.localG)
            x_mean = np.mean(x)
            # Create empty vector to fill
            self.labs = np.empty(len(x)) * np.nan
            # Outliers
            self.labs[(self.localG < Eij_mean) &
                      (y > y_mean) &
                      (self.p_sim <= self.sig)] = 1
            # Clusters
            self.labs[(self.localG < Eij_mean) &
                      (y < y_mean) &
                      (self.p_sim <= self.sig)] = 2
            # Other
            self.labs[(self.localG > Eij_mean) &
                      (self.p_sim <= self.sig)] = 3
            # Non-significant
            self.labs[self.p_sim > self.sig] = 4

        del (self.keep_simulations, self.n_jobs,
             self.permutations, self.seed, self.rlocalG,
             self.connectivity, self.labels)

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

```

</details>

## Local Geary (multivariate) - might need to remove?

Links to:
- [Documentation](https://github.com/jeffcsauer/GSOC2020/blob/master/docs/localgeary.ipynb)
- [Location in PySAL `esda`](https://github.com/pysal/esda/blob/master/esda/local_geary_mv.py)
- [Location in GSOC 2020 workbook](https://github.com/jeffcsauer/GSOC2020/blob/master/functions/local_geary_mv.py)

Stable release version:

<details>
  <summary>Click to expand code</summary>
  
``` python

```

</details>