#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdint>
#include <array>
#include "tipsy_file_reader.h"

/* Structure definitions */
// The following two classes are used in function readTipsyFile,
// which is used to read particle data created by MAGI.
typedef struct magi_tipsy_header {
   double time;
   int nbodies;
   int ndim;
   int nsph;
   int ndark;
   int nstar;
} Magi_tipsy_header;

typedef struct magi_tipsy_particle {
   float mass;
   float pos[3];
   float vel[3];
   float eps;
   int idx;
} Magi_tipsy_particle;


template <class T>
T byteswap(const T val) {
    constexpr int size = sizeof(T);
    constexpr int block_size = sizeof(uint16_t);
    if (size > block_size) {
        assert(size % block_size == 0);
        constexpr int n_block = size / block_size;
        std::array<uint16_t *, n_block> block;
        T T_tmp = val;
        uint16_t * head = reinterpret_cast<uint16_t *>(&T_tmp);
        for (int i = 0; i < n_block; i++) block[i] = head + i;
        for (int i = 0; i < n_block/2; i++) {
            uint16_t high = *block[i];
            uint16_t low  = *block[n_block - 1 - i];
            *block[n_block - 1 - i] = (high >> 8) | (high << 8);
            *block[i] = (low >> 8) | (low << 8);
        }
        return T_tmp;
    } else {
        return val;
    }
}

extern "C" {

int read_ptcl_number(char file_name[]) {
   // This function is used to read the header information
   // from a file of the TIPSY format  created by
   // MAGI (https://bitbucket.org/ymiki/magi).
   FILE *fp;
   if ((fp = fopen(file_name,"rb")) == NULL) {
      fprintf(stderr,"Cannot open file %s\n.",file_name);
      exit(EXIT_FAILURE);
   }
   Magi_tipsy_header header;
   fread(&header, sizeof(Magi_tipsy_header), 1, fp);
#ifdef READ_DATA_WITH_BYTESWAP
   header.nbodies = byteswap(header.nbodies);
#endif
   fclose(fp);
   return header.nbodies;
}

void read_ptcl_data(char file_name[],
                    int n_ptcl,
                    float *mass,
                    float *pos,
                    float *vel) {
   // This function is used to read particle data from a file
   // of the TIPSY format created by
   // MAGI (https://bitbucket.org/ymiki/magi).
   FILE *fp;
   if ((fp = fopen(file_name,"rb")) == NULL) {
      fprintf(stderr,"Cannot open file %s\n.",file_name);
      exit(EXIT_FAILURE);
   }
   Magi_tipsy_header header;
   fread(&header, sizeof(Magi_tipsy_header), 1, fp);
   int i;
   for (i = 0; i < header.nbodies; i++) {
      Magi_tipsy_particle ptcl_tipsy;
      fread(&ptcl_tipsy,sizeof(Magi_tipsy_particle),1,fp);
#ifdef READ_DATA_WITH_BYTESWAP
      mass[i]      = byteswap(ptcl_tipsy.mass);
      pos[3*i]     = byteswap(ptcl_tipsy.pos[0]);
      pos[3*i + 1] = byteswap(ptcl_tipsy.pos[1]);
      pos[3*i + 2] = byteswap(ptcl_tipsy.pos[2]);
      vel[3*i]     = byteswap(ptcl_tipsy.vel[0]);
      vel[3*i + 1] = byteswap(ptcl_tipsy.vel[1]);
      vel[3*i + 2] = byteswap(ptcl_tipsy.vel[2]);
#else
      mass[i]      = ptcl_tipsy.mass;
      pos[3*i    ] = ptcl_tipsy.pos[0];
      pos[3*i + 1] = ptcl_tipsy.pos[1];
      pos[3*i + 2] = ptcl_tipsy.pos[2];
      vel[3*i    ] = ptcl_tipsy.vel[0];
      vel[3*i + 1] = ptcl_tipsy.vel[1];
      vel[3*i + 2] = ptcl_tipsy.vel[2];
#endif
   }
   fclose(fp);
}

} // END of extern "C"
