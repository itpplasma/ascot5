/**
 * @file mccc_coefs.h
 * @brief Routines to evaluate coefficients needed to evaluate collisions.
 */
#ifndef MCCC_COEFS_H
#define MCCC_COEFS_H

#include <math.h>
#include <float.h>
#include "../../ascot5.h"
#include "../../consts.h"
#include "mccc.h"

/**
 * @brief Evaluate collision parameter [kg^2 m^3 / s^4]
 *
 *\f$c_{ab} = \frac{n_b q_a^2q_b^2 \ln\Lambda_{ab}}{4\pi\epsilon_0^2}\f$
 *
 * where
 *
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_cab(qa, qb, nb, clogab) (                            \
        nb * qa*qa * qb*qb * clogab / ( 4 * CONST_PI * CONST_E0*CONST_E0 ) )

/**
 * @brief Evaluate non-relativistic drag coefficient [m/s^2]
 *
 *\f$Q =-c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0(v_a/v_b) / (m_a m_b v_b^2)\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$m_b\f$ is plasma species mass [kg]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_Q(ma, qa, mb, qb, nb, vb, clogab, mu0) (           \
        -mccc_coefs_cab(qa, qb, nb, clogab) * mu0 / ( ma * mb * vb*vb ) )

/**
 * @brief Evaluate derivative of non-relativistic drag coefficient [m/s^2]
 *
 *\f$Q'=-c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0'(v_a/v_b) / (m_a m_b v_b^2)\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$m_b\f$ is plasma species mass [kg]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_dQ(ma, qa, mb, qb, nb, vb, clogab, dmu0) (           \
        -mccc_coefs_cab(qa, qb, nb, clogab) * dmu0 / ( ma * mb * vb*vb*vb ) )

/**
 * @brief Evaluate non-relativistic friction coefficient [m/s^2]
 *
 *\f$F = -c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\left(m_a^{-1} + m_b^{-1}\right)
  \mu0(v_a/v_b) / (m_a v_b^2)\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$m_b\f$ is plasma species mass [kg]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_F(ma, qa, mb, qb, nb, vb, clogab, mu0) (             \
        -( 1/mb + 1/ma ) * mccc_coefs_cab(qa, qb, nb, clogab) * mu0     \
        / ( ma * vb*vb ) )

/**
 * @brief Evaluate non-relativistic parallel diffusion coefficient [m^2/s^3]
 *
 *\f$D_\parallel=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0(v_a/v_b)
  /(2m_a^2v_a)\f$
 *
 * or \f$D_\parallel = 4c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})/(6\sqrt{\pi}m_a^2vb)
      \;\mathrm{when}\; v_a=0\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_Dpara(ma, qa, va, qb, nb, vb, clogab, mu0) (         \
        ( va > 0 ) ?                                                    \
        mccc_coefs_cab(qa, qb, nb, clogab) * mu0 / ( 2 * ma*ma * va ) : \
        mccc_coefs_cab(qa, qb, nb, clogab) * 4                          \
        / ( 6 * CONST_SQRTPI * ma*ma * vb ) )

/**
 * @brief Evaluate derivative of non-relativistic parallel diffusion
 * coefficient [m/s^2]
 *
 * \f$D_\parallel'=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})(\mu_0'(v_a/v_b)/v_b
  - \mu(v_a/v_b)/v_a)/(2m_a^2v_a)\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_dDpara(ma, qa, va, qb, nb, vb, clogab, mu0, dmu0) (  \
        mccc_coefs_cab(qa, qb, nb, clogab) * ( dmu0/vb - mu0/va )       \
        / ( 2 * ma*ma * va ) )

/**
 * @brief Evaluate non-relativistic perpendicular diffusion
 * coefficient [m^2/s^3]
 *
 *\f$D_\perp=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_1(v_a/v_b)
  /(2m_a^2v_a)\f$
 *
 * or \f$D_\perp = 4c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})/(6\sqrt{\pi}m_a^2vb)
      \;\mathrm{when}\; v_a=0\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_Dperp(ma, qa, va, qb, nb, vb, clogab, mu1) (         \
        ( va > 0 ) ?                                                    \
        mccc_coefs_cab(qa, qb, nb, clogab) * mu1 / ( 2 * ma*ma * va ) : \
        mccc_coefs_cab(qa, qb, nb, clogab) * 4                          \
        / ( 6 * CONST_SQRTPI * ma*ma * vb ) )

/**
 * @brief Evaluate guiding center drag coefficient [m/s^2]
 *
 *\f$K = Q + D_\parallel' + 2D_\parallel/va\f$
 *
 * where
 *
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$D_parallel\f$ is  [1/s]
 * - \f$Q\f$ is  [m/s^2]
 */
#define mccc_coefs_K(va, Dpara, dDpara, Q) (    \
        Q + dDpara + 2*Dpara / va )

/**
 * @brief Evaluate pitch collision frequency [1/s]
 *
 *\f$\nu = 2D_\perp/v_a^2\f$
 *
 * where
 *
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$D_perp\f$ is  []
 */
 #define mccc_coefs_nu(va, Dperp) (             \
        2 * Dperp / ( va * va ) )

/**
 * @brief Evaluate spatial diffusion coefficient [m^2/s]
 *
 *\f$D_X = (\frac{1}{2}(D_\parallel-D_\perp)(1-\xi^2) + D_\perp)/\omega_g^2\f$
 *
 * where
 *
 * - \f$xi\f$ is test particle pitch
 * - \f$D_\parallel\f$ is parallel diffusion coefficient [m^2/s^3]
 * - \f$D_perp\f$ is perpendicular diffusion coefficient [m^2/s^3]
 * - \f$\omega_g\f$ is gyrofrequency [1/s]
 */
#define mccc_coefs_DX(xi, Dpara, Dperp, gyrofreq) (              \
        ( 0.5 * ( Dpara - Dperp ) * ( 1 - xi*xi ) + Dperp )      \
        / (gyrofreq*gyrofreq) )

GPU_DECLARE_TARGET_SIMD_UNIFORM(mdata)
static void mccc_coefs_mufun(real mufun[3], real x, mccc_data* mdata);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(nspec, mb, qb, nb, Tb)
inline static void mccc_coefs_clog(
        real* clogab, real ma, real qa, real va, int nspec,
        const real* mb, const real* qb, const real* nb, const real* Tb);
DECLARE_TARGET_END

/**
 * @brief Evaluate Coulomb logarithm.
 *
 * Coulomb logarithm is evaluated separately with respect to each plasma
 * species. It is calculated as a logarithm of the ratio of maximum and
 * minimum impact parameters. Maximum impact parameter is the Debye length
 * and minimum impact parameter is either classical particle radius or
 * inverse of De Broglie wavelength.
 *
 * @param clogab array where evaluated values for Coulomb logarithm are stored.
 * @param ma test particle mass [kg]
 * @param qa test particle charge [C]
 * @param va test particle velocity [m/s]
 * @param nspec number of plasma species
 * @param mb plasma species masses [kg]
 * @param qb plasma species charges [C]
 * @param nb plasma species densities [m^-3]
 * @param Tb plasma species temperatures [J]
 */
inline static void mccc_coefs_clog(
        real* clogab, real ma, real qa, real va, int nspec, const real* mb,
        const real* qb, const real* nb, const real* Tb) {

#ifdef MCCC_CLOG_NRL_REFERENCE_EV
    /* Diagnostic: NRL thermal logarithms at a fixed test-ion reference energy.
     * This isolates the logarithm prescription from the stochastic integrator.
     * Species zero is the electron, as required by the plasma interface. */
    const real proton_mass = 1.6726e-27;
    real a = ma / proton_mass;
    real z = fabs(qa / CONST_E);
    real te = Tb[0] / CONST_E;
    clogab[0] = 24 - log(sqrt(nb[0]*1e-6) / te);
    for(int i = 1; i < nspec; i++) {
        real ab = mb[i] / proton_mass;
        real zb = fabs(qb[i] / CONST_E);
        real ti = Tb[i] / CONST_E;
        real argument = sqrt(nb[i]*1e-6*zb*zb/ti) * z*zb*(a+ab)
            / (a*ti + ab*MCCC_CLOG_NRL_REFERENCE_EV);
        clogab[i] = 23 - log(fmax(DBL_EPSILON, argument));
    }
    return;
#endif

    /* Evaluate Debye length */
    real sum = 0;
    GPU_SEQUENTIAL_LOOP
    for(int i = 0; i < nspec; i++){
        sum += nb[i] * qb[i] * qb[i] / Tb[i];
    }
    real debyeLength = sqrt(CONST_E0/sum);

    /* Evaluate classical and quantum mechanical impact parameter. The one *
     * that is larger is used to evaluate Coulomb logarithm.               */
    GPU_SEQUENTIAL_LOOP
    for(int i=0; i < nspec; i++){
        real vbar = va * va + 2 * Tb[i] / mb[i];
        real mr   = ma * mb[i] / ( ma + mb[i] );
        real bcl  = fabs( qa * qb[i] / ( 4*CONST_PI*CONST_E0 * mr * vbar ) );
        real bqm  = fabs( CONST_HBAR / ( 2 * mr * sqrt( vbar ) ) );

        if(bcl > bqm){
            clogab[i] = log( debyeLength / bcl );
        }
        else{
            clogab[i] = log( debyeLength / bqm );
        }
    }
}

/** @brief Speed derivative of the Coulomb logarithms [s/m]. */
GPU_DECLARE_TARGET_SIMD_UNIFORM(nspec, mb, qb, Tb)
inline static void mccc_coefs_dclog(
        real* dclogab, real ma, real qa, real va, int nspec,
        const real* mb, const real* qb, const real* Tb) {
    GPU_SEQUENTIAL_LOOP
    for(int i = 0; i < nspec; i++) {
#ifdef MCCC_CLOG_NRL_REFERENCE_EV
        dclogab[i] = 0;
#else
        real vbar = va*va + 2*Tb[i]/mb[i];
        real mr = ma*mb[i]/(ma+mb[i]);
        real bcl = fabs(qa*qb[i]/(4*CONST_PI*CONST_E0*mr*vbar));
        real bqm = fabs(CONST_HBAR/(2*mr*sqrt(vbar)));
        dclogab[i] = (bcl > bqm ? 2 : 1)*va/vbar;
#endif
    }
}
DECLARE_TARGET_END

/**
 * @brief Evaluate special functions needed by collision coefficients
 *
 * This function either evaluates the special functions directly or interpolates
 * them from look-up table which should be initialized with mccc_init() before
 * calling this function.
 *
 * Special functions are
 *
 * - mufun[0] = \f$\mu_0(x) = (\mathrm{erf}(x) -2 x \pi^{-1/2} e^{-x^2})/x^2\f$
 * - mufun[1] = \f$\mu_1(x) = \mathrm{erf}(x) - \frac{1}{2}\mu_0(x)\f$
 * - mufun[2] = \f$\mu_0'(x)\f$
 *
 * @param mufun pointer to array where values are stored
 * @param x argument for the special functions
 * @param mdata pointer to mccc data
 */
inline static void mccc_coefs_mufun(real mufun[3], real x, mccc_data* mdata) {

    if(!mdata->usetabulated && fabs(x) < 0.1) {
        /* The direct difference loses precision before its derivative does. */
        real y = x*x;
        mufun[0] = 4*x/CONST_SQRTPI * (1.0/3 + y*(-1.0/5 + y*(1.0/14
            + y*(-1.0/54 + y*(1.0/264 + y*(-1.0/1560 + y/10800))))));
        mufun[1] = erf(x) - 0.5*mufun[0];
        mufun[2] = 4/CONST_SQRTPI * (1.0/3 + y*(-3.0/5 + y*(5.0/14
            + y*(-7.0/54 + y*(9.0/264 + y*(-11.0/1560 + 13*y/10800))))));
    }
    else if(!mdata->usetabulated && x!= 0) {
        real expm2x = exp(-x*x);
        real erfx   = erf(x);

        mufun[0] = ( erfx - 2 * x * expm2x / CONST_SQRTPI ) / (x*x);
        mufun[1] = erfx - 0.5 * mufun[0];
        mufun[2] = 4 * expm2x / CONST_SQRTPI - 2 * mufun[0] / x;
    }
    else if(mdata->usetabulated && x != 0) {
        // TODO implement me
    }
    else {
        mufun[0] = 0;
        mufun[1] = 0;
        mufun[2] = 4 / ( 3 * CONST_SQRTPI );
    }

}
#endif
