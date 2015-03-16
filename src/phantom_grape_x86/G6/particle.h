// Index is address.
struct Particle{
  REAL  t;    // this particle's t, dt, mass, pot
  REAL  dt;   //
  float mass; //
  float pot;  //
  // ------------------------
  REAL  xpos; // if twobody,
  REAL  ypos; // the time of x etc. is rl particles.
  REAL  zpos; // t, dt are cm's.
  REAL  xvel; // It is confused.
  REAL  yvel;
  REAL  zvel;
  REAL  xacc;
  REAL  yacc;
  REAL  zacc;
  float xjrk;
  float yjrk;
  float zjrk;
  // -------------------------
  unsigned int id; // ID of this particle
  unsigned int in; // ID of the nearest particle (twobody: -1)
  int nn;
  int *nptr;
};

// index is ID
// ID -> address
struct Address{
  unsigned int jgp; // Which node ?
  unsigned int adr; // Which address in the node ?
};

// Index is integration NO.
// integration NO -> address
struct Actives{
  int ni;
  int index[NMAX-1];
  int nisg;
  int indexsg[NMAX-1];
  int nicm;
  int indexcm[NMAX-1];
  int nirl;
  int indexrl[NMAX-1];
};

struct Predposvel{
  REAL  xpos, ypos, zpos;
  REAL  xvel, yvel, zvel;
  float id;
  float mass; // only in rl twobody force
  float pad[2];
};

struct NewAccJerk{
  REAL  xacc;
  REAL  yacc;
  REAL  zacc;
  float pot;
  float xjrk;
  float yjrk;
  float zjrk;
  unsigned int in;
  int nn;
  int *nptr;
};
