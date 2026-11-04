#include <math.h>
#include <stdio.h>
#include "../physlib.h"

int main(void) {
    const real mass = 1.6726e-27;
    const real field = 2.0;
    const real velocities[] = {1.0, 1e4, 1e6, 1e7};
    const real pitches[] = {-1.0, -0.7, 0.0, 0.7, 1.0};
    for(int i = 0; i < 4; i++) {
        real momentum = physlib_pnorm_vnorm(mass, velocities[i]);
        for(int j = 0; j < 5; j++) {
            real parallel = physlib_gc_ppar(momentum, pitches[j]);
            real mu = physlib_gc_mu(mass, momentum, pitches[j], field);
            real recovered = physlib_gc_p(mass, mu, parallel, field);
            real pitch = physlib_gc_xi(mass, mu, parallel, field);
            if(!isfinite(recovered) || !isfinite(pitch)
               || fabs(recovered / momentum - 1) > 1e-14
               || fabs(pitch) > 1) {
                fprintf(stderr, "GC momentum invariant failed at v=%g xi=%g\n",
                        velocities[i], pitches[j]);
                return 1;
            }
        }
    }
    puts("GC momentum and pitch invariants: PASS");
    return 0;
}
