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

# Reference
- Work the package is used for:\
  X. Ou and J. Xia, *SuperDC: Stable superfast divide-and-conquer eigenvalue decomposition*, submitted (2021), arXiv:2108.04209.
- Earlier preliminary work:\
  J. Vogel, J. Xia, S. Cauley, and V. Balakrishnan, *Superfast divide-and-conquer method and perturbation analysis for structured eigenvalue solutions*,\
  SIAM J. Sci. Comput., 38 (2016), pp. A1358-A1382. (PDF--erratum. Journal article link)
- Work in progress:\
  J. Vogel, J. Xia, Z. Xin, and X. Ou, *Structured numerical computations via superfast eigenvalue decompositions*, under preparation, 2021.

# Usage

Other package needed:
- FMM1D (https://github.com/fastsolvers/FMM1D)
- HSS (https://github.com/fastsolvers/HSS)

Place the following folders under \<superdc\>:<br>
  \<superdc_1.0.0\>, \<fmm1d\>, \<hss\>

## tests/
- test_band.m: A banded matrix test example
- test_dense.m A dense matrix test example

## src/superdc.m
- Main eigenvalue decomposition routine
- Accepted inputs: Hermitian HSS matrices

## src/superdcmv.m
- Application of the eigenvector matrix Q or Q* to a vector



