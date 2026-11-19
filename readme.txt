
VSPLINE git repo

https://github.com/JeffFessler/vspline

This code dates from 1990:
https://netlib.org/gcv/vspline

Posting on git in 2026 for legacy reasons.

J. A. Fessler,
"Nonparametric fixed-interval smoothing with vector splines,"
IEEE Trans. Sig. Proc., 39(24):852-9, Apr. 1991

https://doi.org/10.1109/78.80907


Legacy netlib installation instructions:

(1) Edit Defs_make in accordance with the instructions therein.

(2) Change `defs*.h` if necessary for your compiler / environment.
e.g., on a sun 4, add `#include <alloca.h>`

(3) Type make.

(4) make the programs in the `Test/` directory (`cd Test; make`)
and verify them using the examples in `Data/`.

(5) To use, simply include `$(LIB_VSPLINE)` or `$(VLIBS)` (see `List_vspline`)
when linking; see examples in `Test/Makefile`.

(6) For nonlinear problems, you must implement a subroutine
that provides the nonlinear function and its derivatives.
For documentation, see Latex file in `Doc/`.  This is tricky.

This software can also (probably) be used with MATLAB `mex` files
by adding the flag `-DMATLAB` to `CFLAGS` in `matrix/Makefile`.
