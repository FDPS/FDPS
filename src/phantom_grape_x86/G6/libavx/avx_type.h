#define NPIPES 48
#define MAXLEN 1024

typedef struct prdposvel{
  double xpos, ypos, zpos;
  float  xvel, yvel, zvel;
  float  id;
  float  eps2;
  float  h2;
  float  pad[5];
} PrdPosVel;
typedef PrdPosVel * pPrdPosVel;

typedef struct newaccjrk{
  double xacc;
  double yacc;
  double zacc;
  float  pot;
  float  xjrk;
  float  yjrk;
  float  zjrk;
  int    nnb;
  float  rnnb;
  float  pad[4];
} NewAccJrk;
typedef NewAccJrk * pNewAccJrk;
