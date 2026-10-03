import numpy as np

class PCA:
    """
    Principal Component Analysis via SVD.

    Attributes after fit:
      components_               (n_features, n_components) principal axes (columns)
      scores_                   (n_samples,  n_components)  projected data
      explained_variance_       (n_components,)             eigenvalues
      explained_variance_ratio_ (n_components,)             variance ratios
      mean_                     (n_features,)               feature means
    """
    def __init__(self, data=None, n_components=None):
        """
        Parameters
        ----------
        n_components : int or None
            Number of principal components to keep. If None, keep all.
        """
        self.n_components = n_components
        if data is not None:
             self.fit(data)
    def fit(self, X):
        """
        Fit the model with X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        """
        X = np.asarray(X, dtype=float)
        # Compute mean
        self.mean_ = X.mean(axis=0)
        # Center data
        X_centered = X - self.mean_

        # Compute SVD
        U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)

        # Truncate if necessary
        if self.n_components is not None:
            Vt = Vt[:self.n_components, :]
            S = S[:self.n_components]
            U = U[:, :self.n_components]

        # components: columns are principal components
        self.components_ = Vt.T  # shape (n_features, n_components)

        # Scores: projection of X onto principal axes
        self.scores_ = U * S  # shape (n_samples, n_components)

        n_samples = X.shape[0]
        # explained variance
        self.explained_variance_ = (S**2) / (n_samples - 1)
        total_var = self.explained_variance_.sum()
        self.explained_variance_ratio_ = self.explained_variance_ / total_var

        return self

    def transform(self, X, subtract_mean=True):
        """
        Project X onto the principal axes.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        subtract_mean : bool, default True
            If True, center data by subtracting the training mean before projection.

        Returns
        -------
        scores : array, shape (n_samples, n_components)
        """
        X = np.asarray(X, dtype=float)
        if subtract_mean:
            X_centered = X - self.mean_
        else:
            X_centered = X
        # Multiply to get scores
        return X_centered.dot(self.components_)

    def inverse_transform(self, scores):
        """
        Reconstruct data from its principal-component scores.

        Parameters
        ----------
        scores : array, shape (n_samples, n_components)

        Returns
        -------
        X_reconstructed : array, shape (n_samples, n_features)
        """
        # Multiply by transpose of components and add mean
        return scores.dot(self.components_.T) + self.mean_
