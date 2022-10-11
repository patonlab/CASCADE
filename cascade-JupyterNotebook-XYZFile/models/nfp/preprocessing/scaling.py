import numpy as np

from sklearn.preprocessing import RobustScaler
from scipy import sparse
from sklearn.utils import check_array
from sklearn.utils.validation import FLOAT_DTYPES
from sklearn.preprocessing.data import _handle_zeros_in_scale

class RobustNanScaler(RobustScaler):
    
    def _check_array(self, X, copy):
        """Makes sure centering is not enabled for sparse matrices."""
        X = check_array(X, accept_sparse=('csr', 'csc'), copy=self.copy,
                        estimator=self, dtype=FLOAT_DTYPES, force_all_finite=False)

        if sparse.issparse(X):
            if self.with_centering:
                raise ValueError(
                    "Cannot center sparse matrices: use `with_centering=False`"
                    " instead. See docstring for motivation and alternatives.")
        return X

    
    def fit(self, X, y=None):
        if sparse.issparse(X):
            raise TypeError("RobustScaler cannot be fitted on sparse inputs")
        X = self._check_array(X, self.copy)
        if self.with_centering:
            self.center_ = np.nanmedian(X, axis=0)

        if self.with_scaling:
            q_min, q_max = self.quantile_range
            if not 0 <= q_min <= q_max <= 100:
                raise ValueError("Invalid quantile range: %s" %
                                 str(self.quantile_range))

            q = np.nanpercentile(X, self.quantile_range, axis=0)
            self.scale_ = (q[1] - q[0])
            self.scale_ = _handle_zeros_in_scale(self.scale_, copy=False)
        return self
