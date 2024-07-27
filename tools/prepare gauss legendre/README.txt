Small project that is ment to generate roots and weights for Gauss-Legendre quadrature to be used in NCS.

Currently it relies on Eigen for finding roots of Legendre polynomials through matrix diagonalization.

It provide accurate enough for calculations in NCS results only for small degrees (up to 10) and I itent on improving upon it in the future.

Usage:
1. build through "prepare gauss legendre.sln"
2. run and redirect output to a header file "prepare gauss legendre.exe" > gl.h
3. put gl.h into NCS directory 