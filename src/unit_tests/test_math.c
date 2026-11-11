/**
 * @file test_math.c
 * @brief Test program for math functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "../math.h"
#include "../gctransform.h"

/**
 * Evaluate some mathematical function to be used for testing purposes
 */
double myfunc(double x){
    return sin(10*x)*exp(-x);
}

/**
 * Main function for the test program.
 */
int main(void) {

    int i;

    /* A uniform Cartesian field has a zero Cartesian Jacobian at every angle. */
    for(i = 0; i < 12; i++) {
        real phi = -3.0 + 0.5*i;
        real c = cos(phi), s = sin(phi);
        real br = 2.4*c + 0.7*s;
        real bp = -2.4*s + 0.7*c;
        real cylindrical[12] = {br, 0, bp, 0, bp, 0, -br, 0,
                                0.4, 0, 0, 0};
        real cartesian[12];
        math_jac_rpz2xyz(cylindrical, cartesian, 0.7 + i, phi);
        real expected[12] = {2.4, 0, 0, 0, 0.7, 0, 0, 0, 0.4, 0, 0, 0};
        for(int j = 0; j < 12; j++) {
            if(fabs(cartesian[j] - expected[j]) > 2e-14) {
                printf("Uniform-field Jacobian failed: angle %g, component %d: %g\n",
                       phi, j, cartesian[j]);
                return 1;
            }
        }
    }

    /* A purely toroidal field of constant magnitude has curvature -Rhat/R
     * and neither grad|B| nor torsion, so the first-order guiding-centre
     * momentum and moment follow from the curvature alone:
     * ppar = ppar0 (1 + rho rhohat_R / R),
     * mu = mu0 - rho ppar0^2 rhohat_R / (m B R). */
    for(i = 0; i < 6; i++) {
        real mass = 1.67262192369e-27, charge = 1.602176634e-19, b0 = 2.5;
        real r = 5.0 + 0.3*i, phi = -3.0 + 1.1*i;
        real pr = mass*(1.0e6 - 3.0e5*i), pphi = mass*(2.0e6 - 4.5e5*i),
             pz = mass*(0.4e6 + 2.0e5*i);
        real toroidal[12] = {0, 0, 0, 0, b0, 0, 0, 0, 0, 0, 0, 0};
        real R, Phi, Z, ppar, mu, zeta;
        gctransform_particle2guidingcenter(mass, charge, toroidal, r, phi,
                                           0.2, pr, pphi, pz, &R, &Phi, &Z,
                                           &ppar, &mu, &zeta);
        real pperp = sqrt(pr*pr + pz*pz);
        real mu0 = pperp*pperp/(2*mass*b0);
        real rho = pperp/(charge*b0), rhohat_r = pz/pperp;
        real ppar_exact = pphi*(1 + rho*rhohat_r/r);
        real mu_exact = mu0 - rho*pphi*pphi*rhohat_r/(mass*b0*r);
        if(fabs(ppar/ppar_exact - 1) > 1e-13 || fabs(mu/mu_exact - 1) > 1e-13) {
            printf("Curvature of a toroidal field failed: case %d, "
                   "ppar/exact-1 %g, mu/exact-1 %g\n", i,
                   ppar/ppar_exact - 1, mu/mu_exact - 1);
            return 1;
        }
    }

    double accurate=(10.0-(sin(100.0)+10.0*cos(100.0))*exp(-10.0))/101.0;
    double approx=math_simpson(myfunc,0.0,10.0,1.e-10);


    printf("Checking simpson rule.\n");
    printf("accurate result: %20.18g\n",accurate);
    printf("numerical value: %20.18g\n",approx);
    printf("relative error: %20.18g\n",fabs((accurate-approx)/accurate));
    if( fabs((accurate-approx)/accurate) > 1.0e-14){
        printf("Fail!\n");
        return 1;
    }
    else {
        printf("ok!\n");
    }

    /*
      srand48(0);
      for(i=0; i < 1000000; i++) {
      real r = math_normal_rand();
      printf("%lf\n", r);
      }
    */

    printf("Testing matrix operations\n");
    real matA[6] = {1, 2,
                    3, 4,
                    5, 6};/* 3 x 2 matrix*/
    real matB[8] = { 7, 8, 9, 10,
                     11, 12, 13, 14};/* 2 x 4 matrix*/
    real matC[12];/* 3 x 4 matrix*/
    math_matmul(matA, matB, 3, 2, 4, matC);
    printf("Matrix multiplication\n");
    printf("Calculated    Correct\n");
    printf(" %g %g %g %g   %g %g %g %g\n",matC[0],matC[1], matC[2], matC[3], 39.0, 49.0,  59.0,  69.0);
    printf(" %g %g %g %g   %g %g %g %g\n",matC[4],matC[5], matC[6], matC[7], 54.0, 68.0,  82.0,  96.0);
    printf(" %g %g %g %g   %g %g %g %g\n",matC[8],matC[9],matC[10],matC[11], 69.0, 87.0, 105.0, 123.0);

    if( matC[ 0] != 39.0 ||  matC[ 1] != 49.0 || matC[ 2] != 59.0 || matC[ 3] != 69.0 ||
        matC[ 4] != 54.0 ||  matC[ 5] != 68.0 || matC[ 6] != 82.0 || matC[ 7] != 96.0 ||
        matC[ 8] != 69.0 ||  matC[ 9] != 87.0 || matC[10] !=105.0 || matC[11] != 123.0 ){
        printf("Fail!\n");
        return 1;
    }
    else {
        printf("ok!\n");
    }

    printf("Testing 3x3 matrix determinant\n");
    real matdeta[9],det;
    for(i=0;i<9;i++){
        matdeta[i] = (float) i+1;
    }
    matdeta[8] = 10.0;

    det = math_determinant3x3(matdeta[0],matdeta[1],matdeta[2],matdeta[3],matdeta[4],matdeta[5],matdeta[6],matdeta[7],matdeta[8] );
    printf("det( [1 2 3; 4 5 6; 7 8 10] ) == %f  (correct: -3)\n",det);
    if( fabs(det+3) > 1.0e-14){
        return 1;
    }
    else {
        printf("ok!\n");
    }

    printf("Testing triple product\n");
    real *a,*b,*c;
    a=&matdeta[0];
    b=&matdeta[3];
    c=&matdeta[6];

    real tp = math_scalar_triple_product( a, b, c);
    printf("dot( cross([1 2 3],[4 5 6]),[7 8 10] ) == %f  (correct: -3)\n",tp);
    if( fabs(tp+3) > 1.0e-14){
        return 1;
    }
    else {
        printf("ok!\n");
    }

    return 0;
}
