/**
 * @file endcond.c
 * @brief Marker simulation end conditions
 *
 * In the absence of errors, marker simulation is ended when marker meets even
 * one of the active end conditions. User can choose which end conditions are
 * active.
 *
 * The end conditions are:
 * - tlim: Marker time has passed the simulation time limit
 *
 * - emin: Marker energy is below minimum value
 *
 * - therm: Energy is below value derived from local thermal electron energy
 *
 * - wall: Marker has intersected wall
 *
 * - rhomin: Marker has reached minimum rho (normalized poloidal flux) value
 *
 * - rhomax: Marker has reached maximum rho value
 *
 * - polmax: The total cumulative distance marker has travelled poloidally
 *   exceeds maximum value
 *
 * - tormax: The total cumulative distance marker has travelled toroidally
 *   exceeds maximum value
 *
 * - cpumax: Marker simulation has exceeded maximum wall time
 *
 * - neutr: Marker has been neutralized by an atomic reaction
 *
 * - ioniz: Marker has been ionized by an atomic reaction
 *
 * - hybrid: Not an end condition per se but used to notate that the guiding
 *   center simulation will be resumed as a gyro-orbit simulation
 *
 * As magnetic field lines have no energy, emin and therm are never checked for
 * them. Guiding centers are the only markers for which hybrid is checked.
 *
 * In the code, the end conditions are represented as bit arrays with each bit
 * corresponding to a specific end condition. Each marker has a field "endcond",
 * and when marker meets an end condition, the corresponding bit is flagged.
 * This way if marker simultaneously meets several end conditions, all can be
 * flagged.
 *
 * Additionally, when marker meets an end condition, its running state is set to
 * False which notates its simulation should be discontinued. If the end
 * condition is wall collision, the ID of the wall element the marker collided
 * with is stored in the marker fields.
 *
 * @todo Error checking would be a good idea
 */
#include <math.h>
#include "endcond.h"
#include "particle.h"
#include "simulate.h"
#include "physlib.h"
#include "consts.h"
#include "math.h"
#include "plasma.h"

/**
 * @brief Check end conditions for FO markers
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 * @param sim pointer to simulation data struct
 */
void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i,
                      sim_data* sim) {

    /* Note which end conditions are set as active.
       Only these ones are checked */
    int active_tlim      = sim->endcond_active & endcond_tlim;
    int active_wall      = sim->endcond_active & endcond_wall;
    int active_emin      = sim->endcond_active & endcond_emin;
    int active_therm     = sim->endcond_active & endcond_therm;
    int active_rhomax    = sim->endcond_active & endcond_rhomax;
    int active_rhomin    = sim->endcond_active & endcond_rhomin;
    int active_polmax    = sim->endcond_active & endcond_polmax;
    int active_tormax    = sim->endcond_active & endcond_tormax;
    int active_cpumax    = sim->endcond_active & endcond_cpumax;
    int active_neutr     = sim->endcond_active & endcond_neutr;
    int active_ioniz     = sim->endcond_active & endcond_ioniz;

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < p_f->n_mrk; i++) {
        if(p_f->running[i]) {

            /* Update bounces if pitch changed sign */
            if( (p_i->p_r[i]*p_i->B_r[i] + p_i->p_phi[i]*p_i->B_phi[i]
                 + p_i->p_z[i]*p_i->B_z[i])
                * (p_f->p_r[i]*p_f->B_r[i] + p_f->p_phi[i]*p_f->B_phi[i]
                   + p_f->p_z[i]*p_f->B_z[i]) < 0 ) {
                if(p_f->bounces[i] > 0) {
                    /* Half bounce */
                    p_f->bounces[i] *= -1;
                }
                else if(p_f->bounces[i] < 0) {
                    /* Bounce complete */
                    p_f->bounces[i] *= -1;
                    p_f->bounces[i] +=  1;
                }
                else {
                    /* Initial bounce */
                    p_f->bounces[i] +=  1;
                }
            }

            /* Check if the marker time exceeds simulation time */
            if(active_tlim) {
                if(!sim->reverse_time && p_f->time[i] > sim->endcond_lim_simtime) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if(sim->reverse_time && p_f->time[i] < sim->endcond_lim_simtime) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if(p_f->mileage[i] > sim->endcond_max_mileage) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if(active_wall) {
                real w_coll = 0;
                int tile = wall_hit_wall(
                    p_i->r[i], p_i->phi[i], p_i->z[i],
                    p_f->r[i], p_f->phi[i], p_f->z[i], &sim->wall_data, &w_coll);
                if(tile > 0) {
                    real w = w_coll;
                    real x = (1-w)*p_i->r[i]*cos(p_i->phi[i])
                             + w*p_f->r[i]*cos(p_f->phi[i]);
                    real y = (1-w)*p_i->r[i]*sin(p_i->phi[i])
                             + w*p_f->r[i]*sin(p_f->phi[i]);
                    real pi[3], pf[3], pxyz_i[3], pxyz_f[3];
                    pi[0] = p_i->p_r[i]; pi[1] = p_i->p_phi[i]; pi[2] = p_i->p_z[i];
                    pf[0] = p_f->p_r[i]; pf[1] = p_f->p_phi[i]; pf[2] = p_f->p_z[i];
                    math_vec_rpz2xyz(pi, pxyz_i, p_i->phi[i]);
                    math_vec_rpz2xyz(pf, pxyz_f, p_f->phi[i]);
                    for(int j = 0; j < 3; j++)
                        pxyz_f[j] = pxyz_i[j] + w*(pxyz_f[j] - pxyz_i[j]);

                    real phi = atan2(y, x) - p_i->phi[i];
                    p_f->time[i] = p_i->time[i] + w*(p_f->time[i] - p_i->time[i]);
                    p_f->mileage[i] = p_i->mileage[i]
                                     + w*(p_f->mileage[i] - p_i->mileage[i]);
                    if(sim->wall_data.type == wall_type_3D) {
                        p_f->r[i] = hypot(x, y);
                        p_f->phi[i] = p_i->phi[i] + atan2(sin(phi), cos(phi));
                    }
                    else {
                        p_f->r[i] = p_i->r[i] + w*(p_f->r[i] - p_i->r[i]);
                        p_f->phi[i] = p_i->phi[i] + w*(p_f->phi[i] - p_i->phi[i]);
                    }
                    p_f->z[i]    = p_i->z[i]    + w*(p_f->z[i] - p_i->z[i]);
                    /* Use the recorded angle: the 2D wall branch interpolates phi. */
                    math_vec_xyz2rpz(pxyz_f, pf, p_f->phi[i]);
                    p_f->p_r[i] = pf[0]; p_f->p_phi[i] = pf[1]; p_f->p_z[i] = pf[2];

                    real B_dB[15], psi[1], rho[2];
                    a5err wallerr = B_field_eval_B_dB(B_dB, p_f->r[i], p_f->phi[i],
                        p_f->z[i], p_f->time[i], &sim->B_data);
                    if(!wallerr) wallerr = B_field_eval_psi(psi, p_f->r[i], p_f->phi[i],
                        p_f->z[i], p_f->time[i], &sim->B_data);
                    if(!wallerr) wallerr = B_field_eval_rho(rho, psi[0], &sim->B_data);
                    if(wallerr) p_f->err[i] = wallerr;
                    else {
                        p_f->rho[i] = rho[0];
                        p_f->B_r[i] = B_dB[0]; p_f->B_r_dr[i] = B_dB[1];
                        p_f->B_r_dphi[i] = B_dB[2]; p_f->B_r_dz[i] = B_dB[3];
                        p_f->B_phi[i] = B_dB[4]; p_f->B_phi_dr[i] = B_dB[5];
                        p_f->B_phi_dphi[i] = B_dB[6]; p_f->B_phi_dz[i] = B_dB[7];
                        p_f->B_z[i] = B_dB[8]; p_f->B_z_dr[i] = B_dB[9];
                        p_f->B_z_dphi[i] = B_dB[10]; p_f->B_z_dz[i] = B_dB[11];
                    }
                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= endcond_wall;
                    p_f->running[i] = 0;
                }
            }

            /* Evaluate marker energy, and check if it is below the minimum
             * energy limit or local thermal energy limit */
            if(active_emin || active_therm) {

                real Ti, vflow;
                a5err errflag =
                    plasma_eval_temp(&Ti, p_f->rho[i], p_f->r[i], p_f->phi[i],
                                     p_f->z[i], p_f->time[i], 1,
                                     &sim->plasma_data);
                plasma_eval_flow(&vflow, p_f->rho[i], p_f->r[i], p_f->phi[i],
                                 p_f->z[i], p_f->time[i], &sim->plasma_data);

                /* Error handling */
                if(errflag) {
                    p_f->err[i]     = errflag;
                    p_f->running[i] = 0;
                    Ti = 0;
                    vflow = 0;
                }
                real pnorm = math_normc(
                    p_f->p_r[i], p_f->p_phi[i], p_f->p_z[i]);
                real gamma = physlib_gamma_pnorm(p_f->mass[i], pnorm);
                real vplasma[3];
                real bnorm = math_normc(p_f->B_r[i], p_f->B_phi[i], p_f->B_z[i]);
                vplasma[0] = ( p_f->p_r[i] / ( gamma * p_f->mass[i] )
                    - vflow * p_f->B_r[i] / bnorm );
                vplasma[1] = ( p_f->p_phi[i] / ( gamma * p_f->mass[i] )
                    - vflow * p_f->B_phi[i] / bnorm );
                vplasma[2] = p_f->p_z[i] / ( gamma * p_f->mass[i] )
                    - vflow * p_f->B_z[i] / bnorm;

                real vnorm = math_norm(vplasma);
                pnorm = physlib_pnorm_vnorm(p_f->mass[i], vnorm);
                real ekin  = physlib_Ekin_pnorm(p_f->mass[i], pnorm);

                if( active_emin && (ekin < sim->endcond_min_ekin) ) {
                    p_f->endcond[i] |= endcond_emin;
                    p_f->running[i] = 0;
                }
                if( active_therm && (ekin < (sim->endcond_min_thermal * Ti)) ) {
                    p_f->endcond[i] |= endcond_therm;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if(active_rhomax) {
                if(p_f->rho[i] > sim->endcond_max_rho) {
                    p_f->endcond[i] |= endcond_rhomax;
                    p_f->running[i] = 0;
                }
            }
            if(active_rhomin) {
                if(p_f->rho[i] < sim->endcond_min_rho) {
                    p_f->endcond[i] |= endcond_rhomin;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if(active_tormax) {
                if(fabs(p_f->phi[i]) > sim->endcond_max_tororb) {
                    maxorb |= endcond_tormax;
                }
            }
            if(active_polmax) {
                if(fabs(p_f->theta[i]) > sim->endcond_max_polorb) {
                    maxorb |= endcond_polmax;
                }
                else if( p_f->bounces[i] - 1 >=
                         (int)(sim->endcond_max_polorb / CONST_2PI )) {
                    maxorb |= endcond_polmax;
                }
            }
            if( sim->endcond_torandpol &&
                maxorb & endcond_tormax && maxorb & endcond_polmax ) {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }
            else if(!sim->endcond_torandpol && maxorb) {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }
            /* Check if the time spent simulating this marker exceeds the
             * given limit*/
            if(active_cpumax) {
                if(p_f->cputime[i] > sim->endcond_max_cputime) {
                    p_f->endcond[i] |= endcond_cpumax;
                    p_f->running[i] = 0;
                }
            }
            /* Check if the particle has been neutralized */
            if(active_neutr) {
                if(p_i->charge[i] != 0.0 && p_f->charge[i] == 0.0) {
                    p_f->endcond[i] |= endcond_neutr;
                    p_f->running[i] = 0;
                }
            }

            /* Check if the particle has been ionized */
            if(active_ioniz) {
                if(p_i->charge[i] == 0.0 && p_f->charge[i] != 0.0) {
                    p_f->endcond[i] |= endcond_ioniz;
                    p_f->running[i] = 0;
                }
            }

            /* Zero end condition if error happened in this function */
            if(p_f->err[i]) {
                p_f->endcond[i] = 0;
            }
        }
    }
}

/**
 * @brief Check end conditions for GC markers
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 * @param sim pointer to simulation data struct
 */
void endcond_check_gc(particle_simd_gc* p_f, particle_simd_gc* p_i,
                      sim_data* sim) {
    int i;

    int active_tlim      = sim->endcond_active & endcond_tlim;
    int active_wall      = sim->endcond_active & endcond_wall;
    int active_emin      = sim->endcond_active & endcond_emin;
    int active_therm     = sim->endcond_active & endcond_therm;
    int active_rhomax    = sim->endcond_active & endcond_rhomax;
    int active_rhomin    = sim->endcond_active & endcond_rhomin;
    int active_polmax    = sim->endcond_active & endcond_polmax;
    int active_tormax    = sim->endcond_active & endcond_tormax;
    int active_cpumax    = sim->endcond_active & endcond_cpumax;

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(i = 0; i < p_f->n_mrk; i++) {
        if(p_f->running[i]) {
            /* Update bounces if pitch changed sign */
            if( p_i->ppar[i] * p_f->ppar[i] < 0 ) {
                if(p_f->bounces[i] > 0) {
                    /* Half bounce */
                    p_f->bounces[i] *= -1;
                }
                else if(p_f->bounces[i] < 0) {
                    /* Bounce complete */
                    p_f->bounces[i] *= -1;
                    p_f->bounces[i] +=  1;
                }
                else {
                    /* Initial bounce */
                    p_f->bounces[i] +=  1;
                }
            }

            /* Check if the marker time exceeds simulation time */
            if(active_tlim) {
                if(!sim->reverse_time && p_f->time[i] > sim->endcond_lim_simtime) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if(sim->reverse_time && p_f->time[i] < sim->endcond_lim_simtime) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if(p_f->mileage[i] > sim->endcond_max_mileage) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if(active_wall) {
                real w_coll = 0;
                int tile = wall_hit_wall(p_i->r[i], p_i->phi[i], p_i->z[i],
                                         p_f->r[i], p_f->phi[i], p_f->z[i],
                                         &sim->wall_data, &w_coll);
                if(tile > 0) {
                    real w = w_coll;
                    real x = (1-w)*p_i->r[i]*cos(p_i->phi[i])
                             + w*p_f->r[i]*cos(p_f->phi[i]);
                    real y = (1-w)*p_i->r[i]*sin(p_i->phi[i])
                             + w*p_f->r[i]*sin(p_f->phi[i]);
                    real phi = atan2(y, x) - p_i->phi[i];
                    p_f->time[i] = p_i->time[i] + w*(p_f->time[i] - p_i->time[i]);
                    p_f->mileage[i] = p_i->mileage[i]
                                     + w*(p_f->mileage[i] - p_i->mileage[i]);
                    if(sim->wall_data.type == wall_type_3D) {
                        p_f->r[i] = hypot(x, y);
                        p_f->phi[i] = p_i->phi[i] + atan2(sin(phi), cos(phi));
                    }
                    else {
                        p_f->r[i] = p_i->r[i] + w*(p_f->r[i] - p_i->r[i]);
                        p_f->phi[i] = p_i->phi[i] + w*(p_f->phi[i] - p_i->phi[i]);
                    }
                    p_f->z[i] = p_i->z[i] + w*(p_f->z[i] - p_i->z[i]);
                    p_f->ppar[i] = p_i->ppar[i] + w*(p_f->ppar[i] - p_i->ppar[i]);
                    p_f->mu[i] = p_i->mu[i] + w*(p_f->mu[i] - p_i->mu[i]);
                    p_f->zeta[i] = p_i->zeta[i] + w*(p_f->zeta[i] - p_i->zeta[i]);
                    real B_dB[15], psi[1], rho[2];
                    a5err wallerr = B_field_eval_B_dB(B_dB, p_f->r[i], p_f->phi[i],
                        p_f->z[i], p_f->time[i], &sim->B_data);
                    if(!wallerr) wallerr = B_field_eval_psi(psi, p_f->r[i], p_f->phi[i],
                        p_f->z[i], p_f->time[i], &sim->B_data);
                    if(!wallerr) wallerr = B_field_eval_rho(rho, psi[0], &sim->B_data);
                    if(wallerr) p_f->err[i] = wallerr;
                    else {
                        p_f->rho[i] = rho[0];
                        p_f->B_r[i] = B_dB[0]; p_f->B_r_dr[i] = B_dB[1];
                        p_f->B_r_dphi[i] = B_dB[2]; p_f->B_r_dz[i] = B_dB[3];
                        p_f->B_phi[i] = B_dB[4]; p_f->B_phi_dr[i] = B_dB[5];
                        p_f->B_phi_dphi[i] = B_dB[6]; p_f->B_phi_dz[i] = B_dB[7];
                        p_f->B_z[i] = B_dB[8]; p_f->B_z_dr[i] = B_dB[9];
                        p_f->B_z_dphi[i] = B_dB[10]; p_f->B_z_dz[i] = B_dB[11];
                    }
                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= endcond_wall;
                    p_f->running[i] = 0;
                }
            }

            /* Evaluate marker energy, and check if it is below the minimum
             * energy limit or local thermal energy limit */
            if(active_emin || active_therm) {
                real Bnorm = math_normc(
                    p_f->B_r[i], p_f->B_phi[i], p_f->B_z[i]);


                real Ti, vflow;
                a5err errflag =
                    plasma_eval_temp(&Ti, p_f->rho[i], p_f->r[i], p_f->phi[i],
                                     p_f->z[i], p_f->time[i], 1,
                                     &sim->plasma_data);
                plasma_eval_flow(&vflow, p_f->rho[i], p_f->r[i], p_f->phi[i],
                                 p_f->z[i], p_f->time[i], &sim->plasma_data);

                /* Error handling */
                if(errflag) {
                    p_f->err[i]     = errflag;
                    p_f->running[i] = 0;
                    Ti = 0;
                    vflow = 0;
                }
                real pnorm = physlib_gc_p(
                    p_f->mass[i], p_f->mu[i], p_f->ppar[i], Bnorm);
                real xi = physlib_gc_xi(
                    p_f->mass[i], p_f->mu[i], p_f->ppar[i], Bnorm);
                real vnorm = physlib_vnorm_pnorm(p_f->mass[i], pnorm);
                real vpar = xi * vnorm;
                real vperp2 = (1 - xi * xi) * vnorm * vnorm;
                vnorm = sqrt((vpar - vflow) * (vpar - vflow) + vperp2);
                real gamma = physlib_gamma_vnorm(vnorm);
                real ekin = physlib_Ekin_gamma(p_f->mass[i], gamma);

                if(active_emin && (ekin < sim->endcond_min_ekin) ) {
                    p_f->endcond[i] |= endcond_emin;
                    p_f->running[i] = 0;
                }
                if( active_therm && (ekin < (sim->endcond_min_thermal * Ti)) ) {
                    p_f->endcond[i] |= endcond_therm;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if(active_rhomax) {
                if(p_f->rho[i] > sim->endcond_max_rho) {
                    p_f->endcond[i] |= endcond_rhomax;
                    p_f->running[i] = 0;
                }
            }
            if(active_rhomin) {
                if(p_f->rho[i] < sim->endcond_min_rho) {
                    p_f->endcond[i] |= endcond_rhomin;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if(active_tormax) {
                if(fabs(p_f->phi[i]) > sim->endcond_max_tororb) {
                    maxorb |= endcond_tormax;
                }
            }
            if(active_polmax) {
                if(fabs(p_f->theta[i]) > sim->endcond_max_polorb) {
                    maxorb |= endcond_polmax;
                }
                else if(p_f->bounces[i] - 1 >=
                        (int)(sim->endcond_max_polorb / CONST_2PI)) {
                    maxorb |= endcond_polmax;
                }
            }
            if( sim->endcond_torandpol &&
                maxorb & endcond_tormax && maxorb & endcond_polmax ) {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }
            else if(!sim->endcond_torandpol && maxorb) {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }

            /* Check if the time spent simulating this marker exceeds the
             * given limit*/
            if(active_cpumax) {
                if(p_f->cputime[i] > sim->endcond_max_cputime) {
                    p_f->endcond[i] |= endcond_cpumax;
                    p_f->running[i] = 0;
                }
            }

            /* If hybrid mode is used, check whether this marker meets the hybrid
             * condition. */
            if(sim->sim_mode == 3) {
                if(p_f->rho[i] > sim->endcond_max_rho) {
                    p_f->endcond[i] |= endcond_hybrid;
                    p_f->running[i] = 0;
                }
            }

            /* Zero end condition if error happened in this function */
            if(p_f->err[i]) {
                p_f->endcond[i] = 0;
            }
        }
    }
}

/**
 * @brief Check end conditions for ML markers
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 * @param sim pointer to simulation data struct
 */
void endcond_check_ml(particle_simd_ml* p_f, particle_simd_ml* p_i,
                      sim_data* sim) {
    int i;

    int active_tlim      = sim->endcond_active & endcond_tlim;
    int active_wall      = sim->endcond_active & endcond_wall;
    int active_rhomax    = sim->endcond_active & endcond_rhomax;
    int active_rhomin    = sim->endcond_active & endcond_rhomin;
    int active_polmax    = sim->endcond_active & endcond_polmax;
    int active_tormax    = sim->endcond_active & endcond_tormax;
    int active_cpumax    = sim->endcond_active & endcond_cpumax;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {
            /* Check if the marker time exceeds simulation time */
            if(active_tlim) {
                if(!sim->reverse_time && p_f->time[i] > sim->endcond_lim_simtime) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if(sim->reverse_time && p_f->time[i] < sim->endcond_lim_simtime) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if(p_f->mileage[i] > sim->endcond_max_mileage) {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if(active_wall) {
                real w_coll = 0;
                int tile = wall_hit_wall(p_i->r[i], p_i->phi[i], p_i->z[i],
                                         p_f->r[i], p_f->phi[i], p_f->z[i],
                                         &sim->wall_data, &w_coll);
                if(tile > 0) {
                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= endcond_wall;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if(active_rhomax) {
                if(p_f->rho[i] > sim->endcond_max_rho) {
                    p_f->endcond[i] |= endcond_rhomax;
                    p_f->running[i] = 0;
                }
            }
            if(active_rhomin) {
                if(p_f->rho[i] < sim->endcond_min_rho) {
                    p_f->endcond[i] |= endcond_rhomin;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if(active_tormax) {
                if(fabs(p_f->phi[i]) > sim->endcond_max_tororb) {
                    maxorb |= endcond_tormax;
                }
            }
            if(active_polmax) {
                if(fabs(p_f->theta[i]) > sim->endcond_max_polorb) {
                    maxorb |= endcond_polmax;
                }
            }
            if( sim->endcond_torandpol &&
                maxorb & endcond_tormax && maxorb & endcond_polmax ) {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }
            else if(!sim->endcond_torandpol && maxorb) {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }

            /* Check if the time spent simulating this marker exceeds the
             * given limit*/
            if(active_cpumax) {
                if(p_f->cputime[i] > sim->endcond_max_cputime) {
                    p_f->endcond[i] |= endcond_cpumax;
                    p_f->running[i] = 0;
                }
            }
        }
    }
}

/**
 * @brief Split endcond to an array of end conditions
 *
 * This function splits the bit array end condition to an array
 * where the active end conditions are presented by numbers.
 * Number for each end condition are defined in this function.
 *
 * @param endcond bit array representing marker end conditions
 * @param endconds integer array large enough to hold all end conditions
 */
void endcond_parse(int endcond, int* endconds) {
    int i = 0;

    if(endcond & endcond_tlim)   {endconds[i++] =  1;};
    if(endcond & endcond_emin)   {endconds[i++] =  2;};
    if(endcond & endcond_therm)  {endconds[i++] =  3;};
    if(endcond & endcond_wall)   {endconds[i++] =  4;};
    if(endcond & endcond_rhomin) {endconds[i++] =  5;};
    if(endcond & endcond_rhomax) {endconds[i++] =  6;};
    if(endcond & endcond_polmax) {endconds[i++] =  7;};
    if(endcond & endcond_tormax) {endconds[i++] =  8;};
    if(endcond & endcond_cpumax) {endconds[i++] =  9;};
    if(endcond & endcond_hybrid) {endconds[i++] = 10;};
    if(endcond & endcond_neutr)  {endconds[i++] = 11;};
    if(endcond & endcond_ioniz)  {endconds[i++] = 12;};
}

/**
 * @brief Represent end condition in human-readable format
 *
 * This function takes end condition represented as integer, as given by
 * endcond_parse().
 *
 * @param endcond end condition integer representation
 * @param str end condition as human-readable string
 */
void endcond_parse2str(int endcond, char* str) {
    int endconds[32];
    endcond_parse(endcond, endconds);

    switch(endcond) {
        case 1:
            sprintf(str, "Sim time limit");
            break;
        case 2:
            sprintf(str, "Min energy");
            break;
        case 3:
            sprintf(str, "Thermalization");
            break;
        case 4:
            sprintf(str, "Wall collision");
            break;
        case 5:
            sprintf(str, "Min rho");
            break;
        case 6:
            sprintf(str, "Max rho");
            break;
        case 7:
            sprintf(str, "Max poloidal orbits");
            break;
        case 8:
            sprintf(str, "Max toroidal orbits");
            break;
        case 9:
            sprintf(str, "CPU time exceeded");
            break;
        case 10:
            sprintf(str, "Hybrid condition");
            break;
        case 11:
            sprintf(str, "Neutralization");
            break;
        case 12:
            sprintf(str, "Ionization");
            break;
    }
}
