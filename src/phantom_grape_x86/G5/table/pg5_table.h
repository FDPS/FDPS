// default softning
//#define SFT_FOR_PP          (7.8125e-4) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
//#define SFT_FOR_PM          (2.34375e-2) /* 3.0/SIZE_OF_MESH */
#include "pppm.h"

// table size
#if 1
#define EXP_BIT 4
#define FRC_BIT 5
//#define EXP_BIT 0
//#define FRC_BIT 6
#else
#define EXP_BIT 4
#define FRC_BIT 2
#endif
#define TICK (1<<(23-FRC_BIT))
#define TBL_SIZE (1<<(EXP_BIT+FRC_BIT))

#ifdef GLOBAL_DEFINE
#define EXTERN 
#else
#define EXTERN extern
#endif

// #define ALLOCATE_TABLE

#ifdef ALLOCATE_TABLE
#include <malloc.h>
EXTERN float (*Force_table)[2];
#else
EXTERN float Force_table[TBL_SIZE][2];
#endif

void pg5_gen_force_table( double (*force_func)(double), double rcut);
void pg5_gen_s2_force_table(double sft_for_PP, double sft_for_PM);
void pg5_set_xscale(double);
