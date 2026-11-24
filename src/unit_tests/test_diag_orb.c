/*
 * test_diag_orb.c
 *
 *  Created on: Oct 12, 2020
 *      Author: sjjamsa
 */

#include "../diag/diag_orb.h"
#include "../consts.h"
#include "../simulate.h"

#include <stdio.h>

#define TOLERANCE 1.0e-12


int test_radial_crossings(){

    int nFails;

    nFails = 0;


    real correct,fr,ir,r0,got;

#define RADIAL_CROSSING(FR,IR,R0,CORRECT,GOT) \
    got = diag_orb_check_radial_crossing( FR, IR, R0);\
    if ( fabs( CORRECT - GOT)  > TOLERANCE ) {\
        nFails += 1;\
        printf(" Radial: ini=%g r0=%g final=%g FAIL: got=%g, expected=%g\n    ---> FAIL!\n",IR,R0,FR,GOT,CORRECT);\
    }\
    else {\
    printf(" Radial: ini=%g r0=%g final=%g got=%g, expected=%g ---> OK \n",IR,R0,FR,GOT,CORRECT);\
    }

    printf("Testing radial crossings.\n");

    ir = 0.1; fr = 0.3; r0 = 0.2; correct = 0.5; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.3; fr = 0.1; r0 = 0.2; correct = 0.5; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.4; fr = 0.3; r0 = 0.2; correct = -1.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.9; fr = 1.1; r0 = 1.0; correct = 0.5; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.0; fr = 0.2; r0 = 0.0; correct = 0.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.0; fr = 0.2; r0 = 0.2; correct = -1.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.0; fr = 0.200000000000001; r0 = 0.2; correct = 1.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.0; fr = 0.0; r0 = 0.2; correct = -1.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.4; fr = 0.4; r0 = 0.4; correct = -1.0; RADIAL_CROSSING(fr,ir,r0,correct,got)

    printf("Testing finished with radial crossings.\n");


    return nFails;



}

int test_plane_crossings(){

    // Remember, we work in radians here!

    int nFails = 0;
    real correct,fa,ia,a0,got;

#define PLANE_CROSSING(FA,IA,A0,CORRECT,GOT) \
    got = diag_orb_check_plane_crossing( FA, IA, A0);\
    if ( fabs( CORRECT - GOT)  > TOLERANCE ) {\
        nFails += 1;\
        printf(" Plane: ini=%g r0=%g final=%g FAIL: got=%g, expected=%g\n    ---> FAIL!\n",IA,A0,FA,GOT,CORRECT);\
    }\
    else {\
    printf(" Plane: ini=%g r0=%g final=%g got=%g, expected=%g ---> OK \n",IA,A0,FA,GOT,CORRECT);\
    }

    printf("Testing plane crossings.\n");


    ia = 0.1; fa = 0.3; a0 = 0.2; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1; fa = 0.1; a0 = 0.2; correct = 0.0; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = CONST_PI-0.1; fa = CONST_PI+0.1; a0 = CONST_PI; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = CONST_2PI-0.1; fa = CONST_2PI+0.1; a0 = CONST_2PI; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = -0.1; fa =0.1; a0 = 0.0; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1+3*CONST_2PI; fa = 0.3+3*CONST_2PI; a0 = 0.2; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1+3*CONST_2PI; fa = 0.1+3*CONST_2PI; a0 = 0.2; correct = 0.0; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1-3*CONST_2PI; fa = 0.3-3*CONST_2PI; a0 = 0.2; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1-3*CONST_2PI; fa = 0.1-3*CONST_2PI; a0 = 0.2; correct = 0.0; PLANE_CROSSING(fa,ia,a0,correct,got)



    printf("Testing finished with plane crossings.\n");

    return nFails;



}


int test_update_gc(){

    int nFail = 0;

    int iSIMD;

//diag_orb_update_gc(diag_orb_data* data, particle_simd_gc* p_f,
//                        particle_simd_gc* p_i)

#define RADIALSTEP 0.025

    diag_orb_data data = {
        .mode = DIAG_ORB_POINCARE, .Nfld = DIAG_ORB_GCFIELDS, .Nmrk = NSIMD,
        .Npnt = 1, .record_mode = simulate_mode_gc, .nradialplots = 2,
        .radialdistances = {0.11, 0.21},
    };
    real correct, got;
    real zero[NSIMD] = {0}, rho_i[NSIMD], rho_f[NSIMD];
    integer id[NSIMD], index[NSIMD];
    particle_simd_gc p_i = {
        .r = zero, .phi = zero, .z = zero, .ppar = zero, .mu = zero,
        .zeta = zero, .weight = zero, .charge = zero, .theta = zero,
        .B_r = zero, .B_phi = zero, .B_z = zero, .id = id, .index = index,
        .rho = rho_i, .mileage = rho_i, .n_mrk = NSIMD,
    };
    particle_simd_gc p_f = p_i;
    p_f.rho = rho_f;
    p_f.mileage = rho_f;

    printf("Testing update_gc in poincare mode for radial steps.\n");
    diag_orb_init(&data);

    for ( iSIMD=0; iSIMD<NSIMD; iSIMD++ ){
        rho_i[iSIMD] = RADIALSTEP*iSIMD;
        rho_f[iSIMD] = RADIALSTEP*(iSIMD+1);
        id[iSIMD] = iSIMD + 20;
        index[iSIMD] = iSIMD;
    }


    diag_orb_update_gc(&data, &p_f, &p_i);


    real rhoi,rhof;
    integer imrk,ipoint,idx;

    for ( iSIMD=0; iSIMD<NSIMD; iSIMD++ ){
        imrk   = p_f.index[iSIMD];
        ipoint = data.mrk_pnt[imrk];
        idx    = imrk * data.Npnt + ipoint;
        rhoi = p_i.rho[iSIMD];
        rhof = p_f.rho[iSIMD];
        if(       rhoi < data.radialdistances[0] &&  rhof > data.radialdistances[0]    ) {
            correct = data.radialdistances[0];
        }else if (rhoi < data.radialdistances[1] &&  rhof > data.radialdistances[1]   ) {
            correct = data.radialdistances[1];
        }else{
            correct = 0.0;
        }

        got = data.rho[idx];
        if ( fabs( got-correct ) > TOLERANCE ){
            nFail++;
            printf(" Rho %6g->%6g expected=%g but got=%g.\n FAIL\n",rhoi,rhof,correct,got);
            printf("   time=%g\n",data.mileage[idx]);
            printf(" Index %ld\n",idx);
        }else{
            printf(" Rho %-6g->%-6g expected=%g and got=%g. OK\n",rhoi,rhof,correct,got);
        }




        if (nFail > 0) {
                /* We don't need to print errors again and again...*/
                //break;
        }
    }


    diag_orb_free(&data);


    printf("   The test does not cover the whole modality. Just a single quantity in one modality.\n");

    printf("Finished with update_gc.\n");

    return nFail;
}



static int test_update_fo(void) {
    const real rho[] = {0.9, 1.3, 0.9, 1.0, 1.3, 1.0};
    const real mileage[] = {0.0, 4.0, 8.0, 9.0, 12.0, 15.0};
    const real p_phi[] = {0.0, 8.0, 0.0, 2.0, 8.0, 2.0};
    const real expected_mileage[] = {1.0, 7.0, 9.0, 15.0};
    const real expected_direction[] = {1.0, -1.0, 1.0, -1.0};
    diag_orb_data data = {
        .mode = DIAG_ORB_POINCARE, .Nfld = DIAG_ORB_FOFIELDS, .Nmrk = 1,
        .Npnt = 5, .record_mode = simulate_mode_fo, .nradialplots = 1,
        .radialdistances = {1.0},
    };
    real zero[NSIMD] = {0}, rho_i[NSIMD], rho_f[NSIMD];
    real mileage_i[NSIMD], mileage_f[NSIMD], p_phi_i[NSIMD], p_phi_f[NSIMD];
    integer id[NSIMD] = {1}, index[NSIMD] = {0};
    particle_simd_fo p_i = {
        .r = zero, .phi = zero, .z = zero, .p_r = zero, .p_z = zero,
        .p_phi = p_phi_i, .weight = zero, .charge = zero, .rho = rho_i,
        .theta = zero, .B_r = zero, .B_phi = zero, .B_z = zero, .id = id,
        .mileage = mileage_i, .index = index, .n_mrk = NSIMD,
    };
    particle_simd_fo p_f = p_i;
    p_f.rho = rho_f;
    p_f.mileage = mileage_f;
    p_f.p_phi = p_phi_f;

    diag_orb_init(&data);
    for(int step = 0; step < 5; step++) {
        rho_i[0] = rho[step];
        rho_f[0] = rho[step + 1];
        mileage_i[0] = mileage[step];
        mileage_f[0] = mileage[step + 1];
        p_phi_i[0] = p_phi[step];
        p_phi_f[0] = p_phi[step + 1];
        diag_orb_update_fo(&data, &p_f, &p_i);
    }

    int fails = (data.mrk_pnt[0] != 4);
    for(int i = 0; i < 4; i++) {
        if(fabs(data.rho[i] - 1.0) > TOLERANCE
           || fabs(data.mileage[i] - expected_mileage[i]) > TOLERANCE
           || fabs(data.p_phi[i] - 2.0) > TOLERANCE
           || data.pncrdi[i] != expected_direction[i]) {
            printf(" Full orbit radial crossing %d: rho=%g mileage=%g "
                   "p_phi=%g direction=%g ---> FAIL\n", i, data.rho[i],
                   data.mileage[i], data.p_phi[i], data.pncrdi[i]);
            fails++;
        }
    }
    diag_orb_free(&data);
    return fails;
}

int main(int argc, char **argv) {

    printf("\n\nTesting diag_orb.c\n\n");


    int fails=0;

    fails += test_radial_crossings();

    fails += test_plane_crossings();

    fails += test_update_gc();
    fails += test_update_fo();

    printf("\nFinished testing diag_orb.c\n\n");


    return fails;

}
