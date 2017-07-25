Tensor Network Randomized SVD for Matlab&copy;/Octave&copy;
-------------------------------------------------------------------------------------------------------

Computes an SVD low-rank approximation of a given matrix in the MPO (Tensor Train Matrix) format using a randomized matrix algorithm. Also includes an algorithm for the fast conversion of a sparse matrix into an MPO-form. 

1. Functions
------------

* [UTN,S,VTN]=TNrSVD(ATN,k,q,tol)

Computes a *k*-rank approximation of given matrix *A* in MPO-form *ATN* using a randomized SVD algorithm. Both orthogonal matrices *U* and *V* are computed directly in MPO-form, *q* represents the exponent in the subspace iteration (*A*^T**A*)^*q* and *tol* is a relative tolerance used for the SVD-based MPO-rounding.

* TN=matrix2mpo(A,n)

Converts a given (sparse) matrix *A* into the Tensor Network MPO-form *TN* with given MPO-tensor dimensions as specified by each row of the matrix *n*. 


2. Reference
------------

Computing low-rank approximations of large-scale matrices with the Tensor Network Randomized SVD

Authors: Kim Batselier, Wenjian Yu, Luca Daniel, Ngai Wong
