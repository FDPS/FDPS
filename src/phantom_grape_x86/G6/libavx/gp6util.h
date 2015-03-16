void g6_open(int);
void g6_open_(int *);

void g6_set_tunit(int);
void g6_set_tunit_(int *);

void g6_set_xunit(int);
void g6_set_xunit_(int *);

void g6_reset(int);
void g6_reset_(int *);

void g6_reset_fofpga(int);
void g6_reset_fofpga_(int *);

void g6_close(int);
void g6_close_(int *);

int g6_npipes(void);
int g6_npipes_(void);

void g6_set_j_particle(int, int, int,
		       double, double, double,
		       double *, double *, double *,
		       double *, double *);
void g6_set_j_particle_(int *, int *, int *,
		       double *, double *, double *,
		       double *, double *, double *,
		       double *, double *);

void g6_set_ti(int, double);
void g6_set_ti_(int *, double *);

void g6calc_firsthalf(int, int, int, int *,
		      double (*)[3], double (*)[3],
		      double (*)[3], double (*)[3], double *,
		      double, double *);
void g6calc_firsthalf_(int *, int *, int *, int *,
		      double (*)[3], double (*)[3],
		      double (*)[3], double (*)[3], double *,
		      double *, double *);

int g6calc_lasthalf(int, int, int, int *,
		    double (*)[3], double (*)[3], double, double *,
		    double (*)[3], double (*)[3], double *);
int g6calc_lasthalf_(int *, int *, int *, int *,
		    double (*)[3], double (*)[3], double *, double *,
		    double (*)[3], double (*)[3], double *);

int g6calc_lasthalf2(int, int, int, int *,
		     double (*)[3], double (*)[3], double, double *,
		     double (*)[3], double (*)[3], double *, int *);
int g6calc_lasthalf2_(int *, int *, int *, int *,
		      double (*)[3], double (*)[3], double *, double *,
		      double (*)[3], double (*)[3], double *, int *);

int g6calc_lasthalf2p(int, int, int, int *,
		      double (*)[3], double (*)[3], double, double *,
		      double (*)[3], double (*)[3], double *, int *, float *);
int g6calc_lasthalf2p_(int *, int *, int *, int *,
		       double (*)[3], double (*)[3], double *, double *,
		       double (*)[3], double (*)[3], double *, int *, float *);

int g6calc_lasthalfn(int, int, int, int *,
		     double (*)[3], double (*)[3], double, double *,
		     double (*)[3], double (*)[3], double *);
int g6calc_lasthalfn_(int *, int *, int *, int *,
		      double (*)[3], double (*)[3], double *, double *,
		      double (*)[3], double (*)[3], double *);

int g6calc_lasthalf2n(int, int, int, int *,
		     double (*)[3], double (*)[3], double, double *,
		     double (*)[3], double (*)[3], double *, int *);
int g6calc_lasthalf2n_(int *, int *, int *, int *,
		      double (*)[3], double (*)[3], double *, double *,
		      double (*)[3], double (*)[3], double *, int *);

int g6calc_lasthalf2np(int, int, int, int *,
		     double (*)[3], double (*)[3], double, double *,
		       double (*)[3], double (*)[3], double *, int *, float *);
int g6calc_lasthalf2np_(int *, int *, int *, int *,
		      double (*)[3], double (*)[3], double *, double *,
			double (*)[3], double (*)[3], double *, int *, float *);

void g6_setup_njdata(int, int);
void g6_setup_njdata_(int *, int *);

int g6_read_neighbour_list(int);
int g6_read_neighbour_list_(int *);

int g6_get_neighbour_list(int, int, int, int *, int *);
int g6_get_neighbour_list_(int *, int *, int *, int *, int *);

void g6_debugfunc(int);
void g6_debugfunc_(int *);

void g6_debugfunc_double(double);
void g6_debugfunc_double_(double *);

void g6_dump(double, int, double *, double (*)[3], double (*)[3]);
void g6_dump_(double *,int *, double *, double (*)[3], double (*)[3]);
