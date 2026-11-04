/**
 * @file test_particle_copy.c
 * @brief Copied markers must keep their error flag and queue index
 */
#include <stdio.h>
#include "../ascot5.h"
#include "../particle.h"

int main(void) {
    particle_simd_fo fo1, fo2;
    particle_simd_gc gc1, gc2;
    particle_allocate_fo(&fo1, 2);
    particle_allocate_fo(&fo2, 2);
    particle_allocate_gc(&gc1, 2);
    particle_allocate_gc(&gc2, 2);
    fo1.err[1] = 7;  fo1.index[1] = 11;
    gc1.err[1] = 13; gc1.index[1] = 17;
    fo2.err[0] = 0;  fo2.index[0] = 0;
    gc2.err[0] = 0;  gc2.index[0] = 0;
    particle_copy_fo(&fo1, 1, &fo2, 0);
    particle_copy_gc(&gc1, 1, &gc2, 0);
    if(fo2.err[0] != 7 || fo2.index[0] != 11
       || gc2.err[0] != 13 || gc2.index[0] != 17) {
        fprintf(stderr, "Particle copy lost the error flag or index\n");
        return 1;
    }
    puts("Particle copy keeps error flag and index: PASS");
    return 0;
}
