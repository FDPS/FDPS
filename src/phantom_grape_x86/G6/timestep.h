void initialize_timestep(REAL);
REAL dt_startup(REAL, REAL, struct Particle);
REAL dt_criterion(REAL, struct Particle *, struct Particle *, REAL *, REAL *);
void global_time(REAL *, int, int *, int *, struct Particle *);
void distribute_index(struct Actives *, struct Particle *, struct Address *);
void change_actives(int, struct Actives *);
void reset_time(int, struct Particle *, struct Globaltime *);
