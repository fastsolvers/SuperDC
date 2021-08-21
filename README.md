# SuperDC (Version 1.0.0)
Stable Superfast Divide-and-Conquer Eigensolver for rank-structured Hermitian matrices -- Hermitian matrices with small off-diagonal ranks or numerical ranks. Examples:
- Banded matrices
- Toeplitz matrices (in Fourier space)
- Some discretized integral equations
- Schur complements in direct factorizations of some sparse matrices
- Other dense or sparse matrices with small off-diagonal (numerical) ranks<br>
(The matrices needed to be represented or approximated by a hierarchical semiseparable matrix (HSS) form.)

*Provided as is. No warranty whatsoever. No liability whatsoever.*

# Authors
Any questions or feedbacks are welcome!
- Xiaofeng Ou: ou17  -at-  purdue.edu
- Jianlin Xia: xiaj  -at-  purdue.edu

# References
## To cite the work (the package is used for)
- X. Ou and J. Xia, SuperDC: Stable superfast divide-and-conquer eigenvalue decomposition, submitted (2021), arXiv:2108.04209.
- X. Ou, J. Vogel, J. Xia, and Z. Xin, Efficient numerical computations via superfast eigenvalue decompositions, preprint, 2021.
## Earlier related work
- J. Vogel, J. Xia, S. Cauley, and V. Balakrishnan, Superfast divide-and-conquer method and perturbation analysis for structured eigenvalue solutions, SIAM J. Sci. Comput., 38 (2016), pp. A1358-A1382. 
- D. Cai and J. Xia, A stable matrix version of the fast multipole method: stabilization strategies and examples, submitted.
## HSS related
- J. Xia, S. Chandrasekaran, M. Gu, X. S. Li, Fast algorithms for hierarchically semiseparable matrices, Numer. Linear Algebra Appl., 17 (2010), pp. 953-976.
- J. Xia, Y. Xi, and M. Gu, A superfast structured solver for Toeplitz linear systems via randomized sampling, SIAM J. Matrix Anal. Appl., 33 (2012), pp. 837-858.
- J. Xia, S. Chandrasekaran, M. Gu, X. S. Li, Superfast multifrontal method for large structured linear systems of equations, SIAM J. Matrix Anal. Appl., 31 (2009), pp. 1382-1411.
- X. Liu, J. Xia, and M. V. de Hoop, Parallel randomized and matrix-free direct solvers for large structured dense linear systems, SIAM J. Sci. Comput., 38 (2016), pp. S508-S538. (PDF. Journal article link)

# Usage

Other package needed:
- FMM1D (https://github.com/fastsolvers/FMM1D)
- HSS (https://github.com/fastsolvers/HSS)

Place the following folders under \<superdc\>:<br>
  \<superdc_1.0.0\>, \<fmm1d\>, \<hss\>

## superdc_1.0.0/tests/
- test_band.m: A banded matrix test example
- test_dense.m A dense matrix test example

## superdc_1.0.0/src/superdc.m
- Main eigenvalue decomposition routine
- Accepted inputs: Hermitian HSS matrices

## superdc_1.0.0/src/superdcmv.m
- Application of the eigenvector matrix Q or Q* to a vector



