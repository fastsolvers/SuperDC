# SuperDC (Version 1.0.0)
Stable Superfast Divide and Conquer Eigensolver for
- Hierarchical semiseparable matrix (HSS)
- Banded matrix with small bandwidth
- Block tridiagonal matirx with small subdiagonal block
- Toeplitz matrix


# Authors
Any questions or feedbacks are welcomed!
- Jianlin Xia: xiaj@purdue.edu
- Xiaofeng Ou: ou17@purdue.edu


# Reference
- Work the package is used for:\
  X. Ou and J. Xia, *SuperDC: Stable superfast divide-and-conquer eigenvalue decomposition*, submitted, 2021.
- Earlier preliminary work:\
  J. Vogel, J. Xia, S. Cauley, and V. Balakrishnan, *Superfast divide-and-conquer method and perturbation analysis for structured eigenvalue solutions*,\
  SIAM J. Sci. Comput., 38 (2016), pp. A1358-A1382. (PDF--erratum. Journal article link)
- Work in progress:\
  J. Vogel, J. Xia, Z. Xin, and X. Ou, *Structured numerical computations via superfast eigenvalue decompositions*, under preparation, 2021.

# Usage

Other package needed (included): FMM1D (https://github.com/fastsolvers/FMM1D)\
Suggested package: HSS

## tests/test_band.m
- A banded matrix test example

## superdc.m
- Main eigenvalue decomposition routine
- Accepted inputs: Hermitian HSS matrices

## superdcmv.m
- Application of the eigenvector matrix Q or Q* to a vector



