/**
 * @file test_suzuki.c
 * @brief Energy ranges of the Suzuki beam-stopping fit.
 *
 * Below 9 keV/amu the fit must fail rather than use the high-energy table.
 * The reference at 9.58 keV/amu is the low-energy table of Suzuki et al.,
 * Plasma Phys. Control. Fusion 40 (1998) 2097, Table 3, extrapolated below
 * its 10 keV/amu limit.
 */
#include <stdio.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../suzuki.h"

static int evaluate(real ekev, real* sigma) {
    real ni[1] = {4e19};
    int anum[1] = {2}, znum[1] = {1};
    return (int) suzuki_sigmav(sigma, ekev * 1e3 * CONST_E, 1.0, 4e19,
                               1500 * CONST_E, 1, ni, anum, znum);
}

int main(void) {
    real sigma = 0;
    int fail = 0;
    if(evaluate(9.58, &sigma) || fabs(sigma / 9.8326043781e-20 - 1) > 1e-9) {
        printf("9.58 keV/amu: %.10e, expected 9.8326043781e-20\n", sigma);
        fail = 1;
    }
    if(!evaluate(5.0, &sigma)) {
        printf("5 keV/amu evaluated without error: %.10e\n", sigma);
        fail = 1;
    }
    if(!evaluate(2e4, &sigma)) {
        printf("20000 keV/amu evaluated without error\n");
        fail = 1;
    }
    real ni[2] = {4e19, 2e18};
    int anum[2] = {2, 20}, znum[2] = {1, 10};
    sigma = 1;
    if(!suzuki_sigmav(&sigma, 30e3 * CONST_E, 1.0, 6e19, 1500 * CONST_E,
                     2, ni, anum, znum) || sigma != 0) {
        printf("Unsupported neon impurity did not return an error and zero rate\n");
        fail = 1;
    }
    printf(fail ? "FAIL\n" : "PASS\n");
    return fail;
}
