#define ALIGN16 __attribute__ ((aligned(16)))
#define ALIGN64 __attribute__ ((aligned(64)))

typedef double v2df __attribute__ ((vector_size(16)));
typedef float  v4sf __attribute__ ((vector_size(16)));
typedef int    v4si __attribute__ ((vector_size(16)));
typedef short  v8hi __attribute__ ((vector_size(16)));

typedef struct iptdata{
	float x[4];
	float y[4];
	float z[4];
	float eps2[4]; // not used in this implementation
} Ipdata, *pIpdata;

typedef struct fodata{
	float ax[4];
	float ay[4];
	float az[4];
	float phi[4];
} Fodata, *pFodata;

typedef struct jpdata{
	float x, y, z, m;
} Jpdata, *pJpdata;

