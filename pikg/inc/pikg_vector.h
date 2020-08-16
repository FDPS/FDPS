#ifndef H_PIKG_VECTOR
#define H_PIKG_VECTOR
typedef struct  {
   float x,y,z;
} pikg_f32vec;

typedef struct  {
   double x,y,z;
} pikg_f64vec;

typedef struct  {
  float x,y;
} pikg_f32vec2;

typedef struct  {
  double x,y;
} pikg_f64vec2;

typedef struct  {
  float x,y,z;
} pikg_f32vec3;

typedef struct  {
  double x,y,z;
} pikg_f64vec3;

typedef struct  {
  float x,y,z,w;
} pikg_f32vec4;

typedef struct  {
  double x,y,z,w;
} pikg_f64vec4;

#endif
