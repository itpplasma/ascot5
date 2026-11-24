#include <math.h>
#include <stdio.h>
#include "../consts.h"
#include "../endcond.h"
#include "../error.h"
#include "../physlib.h"
#include "../simulate/simulate_gc_adaptive.h"

static int test_timestep(const char* label, real initial_step, real time_limit,
                         real mileage_limit, int reverse, real electric_z,
                         int hit_wall, int expect_error) {
    sim_data sim = {0};
    real b[3] = {0, 0, 2}, dB[9] = {0}, e[3] = {0, 0, electric_z};
    sim.B_data.type = B_field_type_TC;
    sim.E_data.type = E_field_type_TC;
    B_TC_init(&sim.B_data.BTC, 1, 0, 0.25, 0.5, b, dB);
    E_TC_init(&sim.E_data.ETC, e);
    sim.enable_ada = 1;
    sim.enable_orbfol = 1;
    sim.fix_usrdef_use = 1;
    sim.fix_usrdef_val = initial_step;
    sim.ada_tol_orbfol = 1e-8;
    sim.ada_max_dphi = 0.1;
    sim.ada_max_drho = 0.1;
    sim.endcond_active = endcond_tlim;
    sim.endcond_lim_simtime = time_limit;
    sim.endcond_max_mileage = mileage_limit;
    sim.reverse_time = reverse;

    particle_gc marker = {
        .r = 2, .energy = 1e4 * CONST_E, .pitch = 0.7, .mass = CONST_U,
        .charge = CONST_E, .anum = 1, .znum = 1, .weight = 1, .id = 1,
    };
    real gamma = physlib_gamma_Ekin(marker.mass, marker.energy);
    real velocity = marker.pitch * physlib_vnorm_gamma(gamma);
    real wall_zmax = velocity * time_limit / 2;
    real wall_r[] = {1, 3, 3, 1};
    real wall_z[] = {-1, -1, wall_zmax, wall_zmax};
    if(hit_wall) {
        sim.endcond_active |= endcond_wall;
        sim.wall_data.type = wall_type_2D;
        sim.wall_data.w2d = (wall_2d_data){
            .n = 4, .wall_r = wall_r, .wall_z = wall_z,
        };
    }
    particle_state state;
    a5err err = particle_input_gc_to_state(&marker, &state, &sim.B_data);
    if(err) {
        fprintf(stderr, "%s: marker initialization failed: %lu\n", label,
                (unsigned long)err);
        return 1;
    }
    particle_state* markers[] = {&state};
    particle_queue queue = {.n = 1, .p = markers};
    simulate_gc_adaptive(&queue, &sim, 1);

    if(expect_error) {
        int err_type, err_line, err_file;
        error_parse(state.err, &err_type, &err_line, &err_file);
        if(err_type == ERR_INVALID_TIMESTEP && state.endcond == 0
           && state.mileage == 0 && state.time == 0) {
            return 0;
        }
    }
    else if(hit_wall) {
        if(!state.err && state.endcond == endcond_wall
           && state.walltile > 0 && state.time == state.mileage
           && fabs(state.time / time_limit - 0.5) < 5e-11
           && fabs(state.z / wall_zmax - 1) < 5e-11) {
            return 0;
        }
    }
    else {
        real mileage = fmin(fabs(time_limit), mileage_limit);
        real time = reverse ? -mileage : mileage;
        real z = velocity * time;
        if(!state.err && state.endcond == endcond_tlim && state.time == time
           && state.mileage == mileage && fabs(state.z - z) < fabs(z)*5e-11) {
            return 0;
        }
    }
    fprintf(stderr, "%s: time=%.17g, mileage=%.17g, z=%.17g, end=%lu, error=%lu\n",
            label, state.time, state.mileage, state.z,
            (unsigned long)state.endcond, (unsigned long)state.err);
    return 1;
}

int main(void) {
    int failed = 0;
    failed += test_timestep("forward time limit", 1e-8, 4e-13, 1, 0, 0, 0, 0);
    failed += test_timestep("reverse time limit", 1e-8, -4e-13, 1, 1, 0, 0, 0);
    failed += test_timestep("mileage limit", 1e-8, 1, 4e-13, 0, 0, 0, 0);
    failed += test_timestep("wall before time limit", 1e-8, 4e-13, 1, 0, 0, 1, 0);
    failed += test_timestep("non-final underflow", 4e-13, 1, 1, 0, 0, 0, 1);
    failed += test_timestep("rejected final step", 1e-8, 4e-13, 1, 0, 1e12, 0, 1);
    if(!failed) {
        puts("Adaptive guiding-centre time limits: PASS");
    }
    return failed != 0;
}
