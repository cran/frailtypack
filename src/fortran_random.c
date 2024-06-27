#include <R.h>
#include <Rinternals.h>

void F77_SUB(updaterandomseed)(int *seed) {
  SEXP setSeedFunc, arg;
  PROTECT(setSeedFunc = findFun(install("set.seed"), R_GlobalEnv));
  PROTECT(arg = ScalarInteger(*seed));
  eval(lang2(setSeedFunc, arg), R_GlobalEnv);
  UNPROTECT(2);
}

void F77_SUB(rndstart)(void) {
  GetRNGstate();
}

void F77_SUB(rndend)(void) {
  PutRNGstate();
}

double F77_SUB(unifrand)(void) {
  return unif_rand();
}

