import numpy
from scipy import sparse
from scipy import stats
from sklearn.base import BaseEstimator
import pysal.lib as lp


class losh(BaseEstimator):
    """Local spatial heteroscedasticity (LOSH)"""

    def __init__(self, connectivity=None, inference=None):
        """
        Initialize a losh estimator

        Arguments
        ---------
        connectivity: scipy.sparse matrix object
                      the connectivity structure describing the relationships
                      between observed units.
        inference: str
                   describes type of inference to be used. options are
                   "chi-square", "permutation", or "simulation".

        Attributes
        ----------
        Hi: numpy array
            Array of LOSH values for each spatial unit.    
        ylag: array
              Spatially lagged y values.
        yresid: array
                Spatially lagged residual values.
        VarHi: array
               Variance of Hi.
        pval: numpy array
              P-values for inference based on either "chi-square", 
              "permutation", or "simulation" approaches.
        """
        
        self.connectivity = connectivity
        self.inference = inference
        
    def fit(self, y, a=None):
        """
        Arguments
        ---------
        y       :   numpy.ndarray
                    array containing continuous data
        a       :   int
                    residual multiplier. Default is 2 in order to generate a
                    variance measure. Users may use 1 for absolute deviations.

        Returns
        -------
        the fitted estimator.

        Notes
        -----
        Technical details and derivations can be found in :cite:`OrdGetis2012`.
        """
        y = np.asarray(y).flatten()
        
        w = self.connectivity

        self.Hi, self.ylag, self.yresid, self.VarHi = self._statistic(y, w, a)
        
        if self.inference is None:
            return self
        elif self.inference == "chi-square":
            print("Note: chi-square inference selected. This assumes a=2.")
            dof = 2/self.VarHi
            Zi = (2*self.Hi)/self.VarHi
            self.pval = 1 - stats.chi2.cdf(Zi, dof)

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
        # Calculate denominator of Hi calculation 
        # as mean of residuals - does np.mean need to be an int?
        denom =  np.mean(yresid) * np.array(rowsum)
        # Carry out final $H_{i}$ calculation by dividing
        # spatial average of residuals by denom
        Hi = lp.weights.lag_spatial(w, yresid) / denom
        
        # Calculate variance
        n = len(y)
        # Calculate average of residuals
        yresid_mean = np.mean(yresid)
        # Calculate VarHi
        squared_rowsum = (n*np.array([np.sum(np.array(list(w[i].values()))**2) for i in w.id_order]))
        VarHi = ((n-1)**-1) * \
                 (denom**-2) * \
                 ((np.sum(yresid**2)/n) - yresid_mean**2) * \
                 (squared_rowsum - (rowsum**2))

        return (Hi, ylag, yresid, VarHi)