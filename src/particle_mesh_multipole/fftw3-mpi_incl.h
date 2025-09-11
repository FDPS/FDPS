#pragma once
/*
 * This file was created for the purpose of wrapping the FFTW3 library.
 * The implementation below is based on fftw3-mpi.h in the library of
 * version 3.3.8 installed on Fugaku as fftw/qowphew via Spack.
 * FDPS DEVELOPERS MUST KEEP THE BELOW IMPLEMENTATION THE SAME AS THAT
 * IN THE ORIGINAL LIBRARY.
 *      Daisuke Namekata (16th Feb, 2021)
 */
#if !defined(PARTICLE_SIMULATOR_WORKAROUND_FFTW3_PROBLEM)
#include <fftw3-mpi.h>
#else // PARTICLE_SIMULATOR_WORKAROUND_FFTW3_PROBLEM

#include "fftw3_incl.h"
#include <mpi.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

struct fftw_mpi_ddim_do_not_use_me {
     ptrdiff_t n;                     /* dimension size */
     ptrdiff_t ib;                    /* input block */
     ptrdiff_t ob;                    /* output block */
};

/*
  huge second-order macro that defines prototypes for all API
  functions.  We expand this macro for each supported precision
 
  XM: name-mangling macro (MPI)
  X: name-mangling macro (serial)
  R: real data type
  C: complex data type
*/

#define FFTW_MPI_DEFINE_API(XM, X, R, C)			\
								\
typedef struct fftw_mpi_ddim_do_not_use_me XM(ddim);		\
								\
FFTW_EXTERN void XM(init)(void);				\
FFTW_EXTERN void XM(cleanup)(void);				\
								\
FFTW_EXTERN ptrdiff_t XM(local_size_many_transposed)		\
     (int rnk, const ptrdiff_t *n, ptrdiff_t howmany,		\
      ptrdiff_t block0, ptrdiff_t block1, MPI_Comm comm,	\
      ptrdiff_t *local_n0, ptrdiff_t *local_0_start,		\
      ptrdiff_t *local_n1, ptrdiff_t *local_1_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size_many)			\
     (int rnk, const ptrdiff_t *n, ptrdiff_t howmany,		\
      ptrdiff_t block0, MPI_Comm comm,				\
      ptrdiff_t *local_n0, ptrdiff_t *local_0_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size_transposed)			\
     (int rnk, const ptrdiff_t *n, MPI_Comm comm,		\
      ptrdiff_t *local_n0, ptrdiff_t *local_0_start,		\
      ptrdiff_t *local_n1, ptrdiff_t *local_1_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size)				\
     (int rnk, const ptrdiff_t *n, MPI_Comm comm,		\
      ptrdiff_t *local_n0, ptrdiff_t *local_0_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size_many_1d)(			\
     ptrdiff_t n0, ptrdiff_t howmany,				\
     MPI_Comm comm, int sign, unsigned flags,			\
     ptrdiff_t *local_ni, ptrdiff_t *local_i_start,		\
     ptrdiff_t *local_no, ptrdiff_t *local_o_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size_1d)(			\
     ptrdiff_t n0, MPI_Comm comm, int sign, unsigned flags,	\
     ptrdiff_t *local_ni, ptrdiff_t *local_i_start,		\
     ptrdiff_t *local_no, ptrdiff_t *local_o_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size_2d)(			\
     ptrdiff_t n0, ptrdiff_t n1, MPI_Comm comm,			\
     ptrdiff_t *local_n0, ptrdiff_t *local_0_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size_2d_transposed)(		\
     ptrdiff_t n0, ptrdiff_t n1, MPI_Comm comm,			\
     ptrdiff_t *local_n0, ptrdiff_t *local_0_start,		\
     ptrdiff_t *local_n1, ptrdiff_t *local_1_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size_3d)(			\
     ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2, MPI_Comm comm,	\
     ptrdiff_t *local_n0, ptrdiff_t *local_0_start);		\
FFTW_EXTERN ptrdiff_t XM(local_size_3d_transposed)(		\
     ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2, MPI_Comm comm,	\
     ptrdiff_t *local_n0, ptrdiff_t *local_0_start,		\
     ptrdiff_t *local_n1, ptrdiff_t *local_1_start);		\
								\
FFTW_EXTERN X(plan) XM(plan_many_transpose)			\
     (ptrdiff_t n0, ptrdiff_t n1,				\
      ptrdiff_t howmany, ptrdiff_t block0, ptrdiff_t block1,	\
      R *in, R *out, MPI_Comm comm, unsigned flags);		\
FFTW_EXTERN X(plan) XM(plan_transpose)				\
     (ptrdiff_t n0, ptrdiff_t n1,				\
      R *in, R *out, MPI_Comm comm, unsigned flags);		\
								\
FFTW_EXTERN X(plan) XM(plan_many_dft)				\
     (int rnk, const ptrdiff_t *n, ptrdiff_t howmany,		\
      ptrdiff_t block, ptrdiff_t tblock, C *in, C *out,		\
      MPI_Comm comm, int sign, unsigned flags);			\
FFTW_EXTERN X(plan) XM(plan_dft)				\
     (int rnk, const ptrdiff_t *n, C *in, C *out,		\
      MPI_Comm comm, int sign, unsigned flags);			\
FFTW_EXTERN X(plan) XM(plan_dft_1d)				\
     (ptrdiff_t n0, C *in, C *out,				\
      MPI_Comm comm, int sign, unsigned flags);			\
FFTW_EXTERN X(plan) XM(plan_dft_2d)				\
     (ptrdiff_t n0, ptrdiff_t n1, C *in, C *out,		\
      MPI_Comm comm, int sign, unsigned flags);			\
FFTW_EXTERN X(plan) XM(plan_dft_3d)				\
     (ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2, C *in, C *out,	\
      MPI_Comm comm, int sign, unsigned flags);			\
								\
FFTW_EXTERN X(plan) XM(plan_many_r2r)				\
     (int rnk, const ptrdiff_t *n, ptrdiff_t howmany,		\
      ptrdiff_t iblock, ptrdiff_t oblock, R *in, R *out,	\
      MPI_Comm comm, const X(r2r_kind) *kind, unsigned flags);	\
FFTW_EXTERN X(plan) XM(plan_r2r)				\
     (int rnk, const ptrdiff_t *n, R *in, R *out,		\
      MPI_Comm comm, const X(r2r_kind) *kind, unsigned flags);	\
FFTW_EXTERN X(plan) XM(plan_r2r_2d)				\
     (ptrdiff_t n0, ptrdiff_t n1, R *in, R *out, MPI_Comm comm,	\
      X(r2r_kind) kind0, X(r2r_kind) kind1, unsigned flags);	\
FFTW_EXTERN X(plan) XM(plan_r2r_3d)				\
     (ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2,			\
      R *in, R *out, MPI_Comm comm, X(r2r_kind) kind0,		\
      X(r2r_kind) kind1, X(r2r_kind) kind2, unsigned flags);	\
								\
FFTW_EXTERN X(plan) XM(plan_many_dft_r2c)			\
     (int rnk, const ptrdiff_t *n, ptrdiff_t howmany,		\
      ptrdiff_t iblock, ptrdiff_t oblock, R *in, C *out,	\
      MPI_Comm comm, unsigned flags);				\
FFTW_EXTERN X(plan) XM(plan_dft_r2c)				\
     (int rnk, const ptrdiff_t *n, R *in, C *out,		\
      MPI_Comm comm, unsigned flags);				\
FFTW_EXTERN X(plan) XM(plan_dft_r2c_2d)				\
     (ptrdiff_t n0, ptrdiff_t n1, R *in, C *out,		\
      MPI_Comm comm, unsigned flags);				\
FFTW_EXTERN X(plan) XM(plan_dft_r2c_3d)				\
     (ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2, R *in, C *out,	\
      MPI_Comm comm, unsigned flags);				\
								\
FFTW_EXTERN X(plan) XM(plan_many_dft_c2r)			\
     (int rnk, const ptrdiff_t *n, ptrdiff_t howmany,		\
      ptrdiff_t iblock, ptrdiff_t oblock, C *in, R *out,	\
      MPI_Comm comm, unsigned flags);				\
FFTW_EXTERN X(plan) XM(plan_dft_c2r)				\
     (int rnk, const ptrdiff_t *n, C *in, R *out,		\
      MPI_Comm comm, unsigned flags);				\
FFTW_EXTERN X(plan) XM(plan_dft_c2r_2d)				\
     (ptrdiff_t n0, ptrdiff_t n1, C *in, R *out,		\
      MPI_Comm comm, unsigned flags);				\
FFTW_EXTERN X(plan) XM(plan_dft_c2r_3d)				\
     (ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2, C *in, R *out,	\
      MPI_Comm comm, unsigned flags);				\
								\
FFTW_EXTERN void XM(gather_wisdom)(MPI_Comm comm_);		\
FFTW_EXTERN void XM(broadcast_wisdom)(MPI_Comm comm_);          \
								\
FFTW_EXTERN void XM(execute_dft)(X(plan) p, C *in, C *out);	\
FFTW_EXTERN void XM(execute_dft_r2c)(X(plan) p, R *in, C *out);	\
FFTW_EXTERN void XM(execute_dft_c2r)(X(plan) p, C *in, R *out);	\
FFTW_EXTERN void XM(execute_r2r)(X(plan) p, R *in, R *out); 



/* end of FFTW_MPI_DEFINE_API macro */

#define FFTW_MPI_MANGLE_DOUBLE(name) FFTW_MANGLE_DOUBLE(FFTW_CONCAT(mpi_,name))
#define FFTW_MPI_MANGLE_FLOAT(name) FFTW_MANGLE_FLOAT(FFTW_CONCAT(mpi_,name))
#define FFTW_MPI_MANGLE_LONG_DOUBLE(name) FFTW_MANGLE_LONG_DOUBLE(FFTW_CONCAT(mpi_,name))

FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_DOUBLE, FFTW_MANGLE_DOUBLE, double, fftw_complex)
FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_FLOAT, FFTW_MANGLE_FLOAT, float, fftwf_complex)
FFTW_MPI_DEFINE_API(FFTW_MPI_MANGLE_LONG_DOUBLE, FFTW_MANGLE_LONG_DOUBLE, long double, fftwl_complex)

#define FFTW_MPI_DEFAULT_BLOCK (0)

/* MPI-specific flags */
#define FFTW_MPI_SCRAMBLED_IN (1U << 27)
#define FFTW_MPI_SCRAMBLED_OUT (1U << 28)
#define FFTW_MPI_TRANSPOSED_IN (1U << 29)
#define FFTW_MPI_TRANSPOSED_OUT (1U << 30)

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif // PARTICLE_SIMULATOR_WORKAROUND_FFTW3_PROBLEM
