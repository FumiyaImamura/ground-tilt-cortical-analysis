# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 23:26:16 2025

@author: Fumiya
"""
import numpy as np

class gcPCA:
    def __init__(self):
        pass
    def CovarianceAndNormalize(self, X):
        # X is trial x time x cells
        if X.ndim != 3:
            raise ValueError("X must be 3D: (trials, time, cells)")
        X = X - X.mean(axis=1, keepdims=True)
        sumXTX = np.einsum('rtc,rtd->cd', X, X)
        normXTX = sumXTX / np.trace(sumXTX)
        return normXTX
    def CovarianceAndNormalize_intraTrial(self,X):
        # X is trial x time x cells
        return self.CovarianceAndNormalize(X)
    def CovarianceAndNormalize_interTrial(self,X):
        # X is trial x time x cells
        X=X.transpose([1,0,2])
        return self.CovarianceAndNormalize(X)
    def fit(self,Xcov):
        Xcov = np.asarray(Xcov)
        Xcov = np.ascontiguousarray(Xcov, dtype=np.float64)
        D,V=np.linalg.eigh(Xcov)
        ind=np.argsort(D)[::-1]
        D=D[ind]
        V=V[:,ind]
        self.explained_variance=D
        self.components_=V
    def transform(self,X):
        X=np.asanyarray(X)
        return X@self.components_

        