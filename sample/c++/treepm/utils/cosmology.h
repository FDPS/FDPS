#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef __COSMOLOGY__
#define __COSMOLOGY__

typedef float COSM_REAL;

struct cosmology {
  COSM_REAL omega_m, omega_v, omega_b, omega_nu, hubble;
};

COSM_REAL fomega(COSM_REAL anow, struct cosmology cosm);
COSM_REAL dladt(COSM_REAL anow, struct cosmology cosm);
COSM_REAL ztotime(COSM_REAL znow, struct cosmology cosm);
COSM_REAL atotime(COSM_REAL znow, struct cosmology cosm);
COSM_REAL timetoz(COSM_REAL tnow, struct cosmology cosm);
COSM_REAL timetoa(COSM_REAL tnow, struct cosmology cosm);
#endif /* __COSMOLOGY__ */

#ifdef __cplusplus
}
#endif
