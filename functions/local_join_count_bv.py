import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.base import BaseEstimator

class Local_Join_Count_BV(BaseEstimator):

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
        BJC_:  numpy.ndarray (1,)
               array containing the estimated Bivariate Local Join Counts ...
        CLC_:  numpy.ndarray (1,)
               array containing the estimated Bivariate Local Join Counts ...
        """

        self.connectivity = connectivity

    def fit(self, x, z, case=None):
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

        self.LJC_ = self._statistic(x, z, w, case=case)

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
            return adj_list_BJC.BJC.values
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