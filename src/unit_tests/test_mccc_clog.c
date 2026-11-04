#include <stdio.h>
#include <math.h>
#define MCCC_CLOG_NRL_REFERENCE_EV 55000.0
#include "../simulate/mccc/mccc_coefs.h"

int main(void) {
    const real proton_mass = 1.6726e-27;
    real mass[3] = {CONST_M_E, proton_mass, 12*proton_mass};
    real charge[3] = {-CONST_E, CONST_E, 6*CONST_E};
    real density[3] = {9.2e19, 8e19, 2e18};
    real temp[3] = {2000*CONST_E, 1000*CONST_E, 1000*CONST_E};
    real reference[3], changed[3];
    /* NRL 2019 thermal-log values, densities in cm^-3 and temperatures in eV. */
    const real expected_electron = 15.524497613053288;
    mccc_coefs_clog(reference, proton_mass, CONST_E, 3e6, 3,
                     mass, charge, density, temp);
    if(fabs(reference[0] - expected_electron) > 1e-10) {
        fprintf(stderr, "NRL electron logarithm: %.15g, expected %.15g\n",
                 reference[0], expected_electron);
        return 1;
    }
    mccc_coefs_clog(changed, proton_mass, CONST_E, 1e6, 3,
                     mass, charge, density, temp);
    for(int i = 0; i < 3; i++) {
        if(fabs(changed[i] - reference[i]) > 1e-12) return 2;
        density[i] *= 4;
    }
    mccc_coefs_clog(changed, proton_mass, CONST_E, 3e6, 3,
                     mass, charge, density, temp);
    for(int i = 0; i < 3; i++) {
        if(fabs(changed[i] - reference[i] + log(2)) > 1e-12) return 3;
    }
    puts("NRL reference-energy Coulomb logarithm: PASS");
    return 0;
}
