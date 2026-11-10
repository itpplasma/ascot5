#include <stdio.h>
#include <math.h>
#include "../simulate/mccc/collision_noise.h"

int main(void) {
    /* Fixture shared with SIMPLE test/test_collision_noise.f90. */
    const double expected[4] = {1.0115071454929696, -0.70803017369489951,
                                1.1794644411780411, -0.65541438399747554};
    const uint64_t key[4][4] = {{0, 1, 0, 3}, {20261110, 1, 0, 4},
                                {20261110, 16384, 1000000, 3}, {7, 42, 123456789, 0}};
    for(int i = 0; i < 4; i++) {
        double x = collision_noise_normal(key[i][0], key[i][1], key[i][2], key[i][3]);
        if(fabs(x - expected[i]) > 4e-15) {
            fprintf(stderr, "fixture %d: %.17g, expected %.17g\n", i, x, expected[i]);
            return 1;
        }
    }
    /* Moments over steps of one marker, and correlation of the two
     * components a collision step uses. */
    const int n = 2000000;
    double s1 = 0, s2 = 0, s4 = 0, c = 0;
    for(int k = 0; k < n; k++) {
        double a = collision_noise_normal(99, 5, k, 3);
        double b = collision_noise_normal(99, 5, k, 4);
        s1 += a; s2 += a*a; s4 += a*a*a*a; c += a*b;
    }
    double mean = s1/n, var = s2/n - mean*mean, kurt = s4/n, corr = c/n;
    printf("mean %.2e var-1 %.2e kurtosis-3 %.2e corr %.2e\n", mean, var - 1, kurt - 3, corr);
    if(fabs(mean) > 4e-3 || fabs(var - 1) > 6e-3 || fabs(kurt - 3) > 3e-2
       || fabs(corr) > 4e-3) {
        return 2;
    }
    puts("counter-based collision noise: PASS");
    return 0;
}
