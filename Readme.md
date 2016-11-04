Genetic Star Discrepancy
===============================

Testing repository for the Genetic Star Discrepancy algorithm.

Note: Since the test data is relatively small and self-contained (~250MB), I've chosen to include it in the repository proper.

In order to clone this repo, you'll need to get the BLIS submodule as well. You
can do this in one of two ways: the simplest is to do a recursive clone, i.e.

```
git clone --recursive https://github.com/chipbuster/genetic-discrepancy
```

If you've already cloned the project, you can get the submodules with two commands:

```
git submodule init
git submodule update
```

Hey you dopes: if we ever plan to release this to public, we'll need to add GPL headers to all the source code files (see: https://www.gnu.org/licenses/gpl-howto.html)


Develop path plan:

1. Use BLIS src to benchmark raw matrix-matrix multiplication.
   Results: about 2 TFLOPS peak for MKL, 1.5 TFLOPS peak for BLIS

2. Use BLISLab to instatiate a derp kernel.
   Ongoing in tests/fitness_kern

3. Transplant BLISLab kernel into actual KNL BLIS kernel.
