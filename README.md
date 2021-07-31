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
X.Ou, J.Xia, *SuperDC: Stable superfast divide-and-conquer eigenvalue decomposition*, submitted

# Usage

## superdc_eigensolver.m
compute the SuperDC eigenvalue decomposition
- accept different parameters, see input/output comments inside

## superdc_matvec.m
apply the eigenmatrix or its transpose to vectors

## superdc_accuracy.m
check the accuracy of the decompostion and the orthogonality of eigenvectors


# License
GNU General Public License v3.0


