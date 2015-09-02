#ifndef __PPPM__
#define __PPPM__

#ifdef N32_GRAPE
#define NUMBER_OF_PART      (32768)
#define SIZE_OF_MESH        (32) /* 2*NUMBER_OF_PART**(1/3) */
#define SIZE_OF_MESH_P2     (34) /* SIZE_OF_MESH+2 */
#define SIZE_OF_GREEN       (17) /* SIZE_OF_MESH/2+1 */
#define SFT_FOR_PP          (3.125e-3) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
#define SFT_FOR_PM          (9.375e-2) /* 3.0/SIZE_OF_MESH */
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (8.7890625e-3)
#define SIZE_OF_SFT         (16384) /* number of bins + 2 */
#define BIN_OF_SFT          (5.36474547e-7) /* RADIUS_FOR_PP2/(SIZE_OF_SFT-1) */
#define SIZE_OF_HOC         (10) /* < 1.0/RADIUS_FOR_PP */
#define SIZE_OF_HOC3        (1000) /* SIZE_OF_HOC*SIZE_OF_HOC*SIZE_OF_HOC */
#define SIZE_OF_CACHE       (4096) /* NUMBER_OF_PART/8 */
#endif

#ifdef N32
#define NUMBER_OF_PART      (32768)
#define SIZE_OF_MESH        (64) /* 2*NUMBER_OF_PART**(1/3) */
#define SIZE_OF_MESH_P2     (66) /* SIZE_OF_MESH+2 */
#define SIZE_OF_GREEN       (33) /* SIZE_OF_MESH/2+1 */
#define SFT_FOR_PP          (3.125e-3) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
#define SFT_FOR_PM          (4.6875e-2) /* 3.0/SIZE_OF_MESH */
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (2.197265625e-3)
#define SIZE_OF_SFT         (16384) /* number of bins + 2 */
#define BIN_OF_SFT          (1.341186367e-7) /*RADIUS_FOR_PP2/(SIZE_OF_SFT-1)*/
#define SIZE_OF_HOC         (21) /* < 1.0/RADIUS_FOR_PP */
#define SIZE_OF_HOC3        (9261) /* SIZE_OF_HOC*SIZE_OF_HOC*SIZE_OF_HOC */
#define SIZE_OF_CACHE       (4096) /* NUMBER_OF_PART/8 */
#endif

#ifdef N64_GRAPE
#define NUMBER_OF_PART      (262144)
#define SIZE_OF_MESH        (64) /* 2*NUMBER_OF_PART**(1/3) */
#define SIZE_OF_MESH_P2     (66) /* SIZE_OF_MESH+2 */
#define SIZE_OF_GREEN       (33) /* SIZE_OF_MESH/2+1 */
#define SFT_FOR_PP          (1.5625e-3) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
#define SFT_FOR_PM          (4.6875e-2) /* 3.0/SIZE_OF_MESH */
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (2.197265625e-3)
#define SIZE_OF_SFT         (16384) /* number of bins + 2 */
#define BIN_OF_SFT          (1.341186367e-7) /*RADIUS_FOR_PP2/(SIZE_OF_SFT-1)*/
#define SIZE_OF_HOC         (15) /* < 1.0/RADIUS_FOR_PP */
#define SIZE_OF_HOC3        (3375) /* SIZE_OF_HOC*SIZE_OF_HOC*SIZE_OF_HOC */
#define SIZE_OF_CACHE       (32768) /* NUMBER_OF_PART/8 */
#endif

#ifdef N64
#define NUMBER_OF_PART      (262144)
#define SIZE_OF_MESH        (128) /* 2*NUMBER_OF_PART**(1/3) */
#define SIZE_OF_MESH_P2     (130) /* SIZE_OF_MESH+2 */
#define SIZE_OF_GREEN       (65) /* SIZE_OF_MESH/2+1 */
#define SFT_FOR_PP          (1.5625e-3) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
#define SFT_FOR_PM          (2.34375e-2) /* 3.0/SIZE_OF_MESH */
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (5.4931640625e-4)
#define SIZE_OF_SFT         (16384) /* number of bins + 2 */
#define BIN_OF_SFT          (3.3529659174e-8) /*RADIUS_FOR_PP2/(SIZE_OF_SFT-1)*/
#define SIZE_OF_HOC         (42) /* < 1.0/RADIUS_FOR_PP */
#define SIZE_OF_HOC3        (74088) /* SIZE_OF_HOC*SIZE_OF_HOC*SIZE_OF_HOC */
#define SIZE_OF_CACHE       (32768) /* NUMBER_OF_PART/8 */
#endif

#ifdef N128_GRAPE
#define NUMBER_OF_PART      (2097152)
#define SIZE_OF_MESH        (128) /* 2*NUMBER_OF_PART**(1/3) */
#define SIZE_OF_MESH_P2     (130) /* SIZE_OF_MESH+2 */
#define SIZE_OF_GREEN       (65) /* SIZE_OF_MESH/2+1 */
#define SFT_FOR_PP          (7.8125e-4) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
#define SFT_FOR_PM          (2.34375e-2) /* 3.0/SIZE_OF_MESH */
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (5.4931640625e-4)
#define SIZE_OF_SFT         (16384) /* number of bins + 2 */
#define BIN_OF_SFT          (3.3529659174e-8) /*RADIUS_FOR_PP2/(SIZE_OF_SFT-1)*/
#define SIZE_OF_HOC         (40) /* < 1.0/RADIUS_FOR_PP */
#define SIZE_OF_HOC3        (64000) /* SIZE_OF_HOC*SIZE_OF_HOC*SIZE_OF_HOC */
#define SIZE_OF_CACHE       (262144) /* NUMBER_OF_PART/8 */
#endif


#ifdef N128
#define NUMBER_OF_PART      (2097152)
#define SIZE_OF_MESH        (256) /* 2*NUMBER_OF_PART**(1/3) */
#define SIZE_OF_MESH_P2     (258) /* SIZE_OF_MESH+2 */
#define SIZE_OF_GREEN       (129) /* SIZE_OF_MESH/2+1 */
#define SFT_FOR_PP          (7.8125e-4) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
#define SFT_FOR_PM          (1.171875e-2) /* 3.0/SIZE_OF_MESH */
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (1.373291015625e-4)
#define SIZE_OF_SFT         (16384) /* number of bins + 2 */
#define BIN_OF_SFT          (8.382414793536e-9) /* RADIUS_FOR_PP2/(SIZE_OF_SFT-2) */
#define SIZE_OF_HOC         (85) /* < 1.0/RADIUS_FOR_PP */
#define SIZE_OF_HOC3        (614125) /* SIZE_OF_HOC*SIZE_OF_HOC*SIZE_OF_HOC */
#define SIZE_OF_CACHE       (262144) /* NUMBER_OF_PART/8 */
#endif

#ifdef N256_GRAPE
#define NUMBER_OF_PART      (16777216)
#define SIZE_OF_MESH        (256) /* 2*NUMBER_OF_PART**(1/3) */
#define SIZE_OF_MESH_P2     (258) /* SIZE_OF_MESH+2 */
#define SIZE_OF_GREEN       (129) /* SIZE_OF_MESH/2+1 */
#define SFT_FOR_PP          (3.90625e-4) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
#define SFT_FOR_PM          (1.171875e-2) /* 3.0/SIZE_OF_MESH */
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (1.373291016e-4)
#define SIZE_OF_SFT         (16384) /* number of bins + 2 */
#define BIN_OF_SFT          (8.382414793536e-9) /* RADIUS_FOR_PP2/(SIZE_OF_SFT-2) */
#define SIZE_OF_HOC         (85) /* < 1.0/RADIUS_FOR_PP */
#define SIZE_OF_HOC3        (614125) /* SIZE_OF_HOC*SIZE_OF_HOC*SIZE_OF_HOC */
#define SIZE_OF_CACHE       (2097152) /* NUMBER_OF_PART/8 */
#endif

#ifdef N256
#define NUMBER_OF_PART      (16777216)
#define SIZE_OF_MESH        (512) /* 2*NUMBER_OF_PART**(1/3) */
#define SIZE_OF_MESH_P2     (514) /* SIZE_OF_MESH+2 */
#define SIZE_OF_GREEN       (257) /* SIZE_OF_MESH/2+1 */
#define SFT_FOR_PP          (3.90625e-4) /* 1.0/(10*NUMBER_OF_PART**(1/3))*/
#define SFT_FOR_PM          (5.859375e-3) /* 3.0/SIZE_OF_MESH */
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (3.4332275390625e-5)
#define SIZE_OF_SFT         (16384) /* number of bins + 2 */
#define BIN_OF_SFT          (2.0956036983840e-9) /*RADIUS_FOR_PP2/(SIZE_OF_SFT-2)*/
#define SIZE_OF_HOC         (170) /* < 1.0/RADIUS_FOR_PP */
#define SIZE_OF_HOC3        (4913000) /* SIZE_OF_HOC*SIZE_OF_HOC*SIZE_OF_HOC */
#define SIZE_OF_CACHE       (2097152) /* NUMBER_OF_PART/8 */
#endif


#endif /* __PPPM__ */
