void initialize_accuracy(REAL, REAL);
void startup_sgmotion(int, int *, struct Particle *);
void correct_sgmotion(int, int *, struct Particle *);
void integrate_exmotion_g6avx(int, double, int, int *, REAL, struct Particle *, struct Address *, void (*corrector)(int, int *, struct Particle *));
void integrate_exmotion_g6avx2(int, double, int, int *, REAL, struct Particle *, struct Address *, void (*corrector)(int, int *, struct Particle *));
void integrate_exmotion_g6avxn(int, double, int, int *, REAL, struct Particle *, struct Address *, void (*corrector)(int, int *, struct Particle *));
void integrate_exmotion_g6avx2n(int, double, int, int *, REAL, struct Particle *, struct Address *, void (*corrector)(int, int *, struct Particle *));
