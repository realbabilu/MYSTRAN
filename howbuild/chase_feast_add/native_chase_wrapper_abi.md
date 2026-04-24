# Native CHASE Wrapper ABI (Current Expectation)

When `MYSTRAN_USE_EXTERNAL_CHASE=ON`, MYSTRAN expects this external routine symbol:

```fortran
subroutine MYSTRAN_CHASE_DSYGV_CRS( &
    N,NNZA,IA,JA,A,NNZB,IB,JB,B,NEV,MAXSUB,TOL,MAXIT,EMIN,EMAX,EVAL,EVEC,RES,NFOUND,INFO)
```

Argument meaning:
- `N`: system size.
- `A`: stiffness matrix `KLL` in CSR (`IA,JA,A`), symmetric real.
- `B`: mass matrix `MLL` in CSR (`IB,JB,B`), symmetric real.
- `NEV`: requested eigenpairs.
- `MAXSUB`: backend subspace size suggestion.
- `TOL`, `MAXIT`: convergence controls.
- `EMIN`, `EMAX`: search interval in eigenvalue domain (`lambda = (2*pi*f)^2`).
- `EVAL`: output eigenvalues (ascending expected).
- `EVEC(N,*)`: output eigenvectors.
- `RES`: residual array (optional diagnostic from backend).
- `NFOUND`: number of converged eigenpairs returned.
- `INFO`: `0` success, nonzero means backend failure/fallback.

Runtime behavior:
- Success (`INFO=0` and `NFOUND>0`) -> MYSTRAN fills `EIGEN_VAL/EIGEN_VEC`.
- Otherwise -> warning `4917` and fallback to ARPACK Lanczos.
