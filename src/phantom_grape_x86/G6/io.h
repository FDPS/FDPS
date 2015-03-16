void input_startup(FILE *, REAL, int *, struct Particle *, struct Address *);
void input_snapshot(FILE *, struct Globaltime *, int *, struct Particle *, struct Energy *);
void output_snapshot(FILE *, struct Globaltime *, int, struct Particle *, struct Energy *);
void output_parameter(FILE *, char *, REAL, REAL);
void output_data(FILE *, int, REAL, struct Energy *, struct Step *);
