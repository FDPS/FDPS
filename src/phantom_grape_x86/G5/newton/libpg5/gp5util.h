#ifndef GP5UTIL_H
#define GP5UTIL_H

#if defined(__cplusplus)
extern "C" {
#endif

#define JMEMSIZE (1<<16) /* 64K */

void g5_set_xmjMC(int devid, int adr, int nj, double (*xj)[3], double *mj);
void g5_set_xmjMC0(int devid, int adr, int nj, double (*xj)[3], double *mj, double *epsj);
void g5_set_nMC(int devid, int n);
void g5_set_xiMC(int devid, int ni, double (*xi)[3]);
void g5_set_xiMC0(int devid, int ni, double (*xi)[3], double *eps);
void g5_runMC(int devid);
void g5_runMC0(int devid);
void g5_get_forceMC(int devid, int ni, double (*a)[3], double *p);
void g5_calculate_force_on_xMC(int devid, double (*x)[3], double (*a)[3], double *p, int ni);
void g5_calculate_force_on_xMC0(int devid, double (*x)[3], double (*a)[3], double *p, int ni, double *eps);

void g5_open(void);
void g5_close(void);
int  g5_get_number_of_pipelines(void);
int  g5_get_jmemsize(void);
void g5_set_range(double xmin, double xmax, double mmin);
void g5_set_n(int n);
void g5_set_mj(int adr, int nj, double *mj);
void g5_set_xj(int adr, int nj, double (*xj)[3]);
void g5_set_xmj(int adr, int nj, double (*xj)[3], double *mj);
void g5_set_xmj0(int adr, int nj, double (*xj)[3], double *mj, double *epsj);
void g5_set_eps_to_all(double eps);
void g5_set_eps(int ni, double *eps);
void g5_set_xi(int ni, double (*xi)[3]);
void g5_run(void);
void g5_get_force(int ni, double (*a)[3], double *p);
void g5_calculate_force_on_x(double (*x)[3], double (*a)[3], double *p, int ni);
void g5_calculate_force_on_x0(double (*x)[3], double (*a)[3], double *p, int ni, double *eps2);

#if defined(__cplusplus)
}
#endif

#endif /* GP5UTIL_H */
