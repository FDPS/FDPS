#pragma once
/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

/*
   I modified the original version of http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c as
   - delete line 47 "#include<stdio.h>"
   - delete line 174 int main(void){...}
   - change N -> MT_N
   - change N -> MT_N
   - change the file name "mt19937ar.c" -> "MT.h"
*/

/* Period parameters */
namespace ParticleSimulator{
    class MT{
    private:
        enum{
            MT_N = 624,
            MT_M = 397,
        };
        enum{ 
            MATRIX_A = 0x9908b0dfUL,   /* constant vector a */
            UPPER_MASK = 0x80000000UL, /* most significant w-r bits */
            LOWER_MASK = 0x7fffffffUL, /* least significant r bits */
        };

        int mti; /* mti==MT_N+1 means mt[MT_N] is not initialized */
        unsigned long mt[MT_N]; /* the array for the state vector  */
        MT() : mti(MT_N+1){}

        ~MT(){}
        MT(const MT &){};
        MT & operator=(const MT &);
        static MT & getInstance(){
            static MT inst;
            return inst;
        }

/* initializes mt[MT_N] with a seed */
#if 0
        void init_genrand(unsigned long s)
        {
            mt[0]= s & 0xffffffffUL;
            for (mti=1; mti<MT_N; mti++) {
                mt[mti] = 
                    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
                /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
                /* In the previous versions, MSBs of the seed affect   */
                /* only MSBs of the array mt[].                        */
                /* 2002/01/09 modified by Makoto Matsumoto             */
                mt[mti] &= 0xffffffffUL;
                /* for >32 bit machines */
            }
        }
#endif

        /* initialize by an array with array-length */
        /* init_key is the array for initializing keys */
        /* key_length is its length */
        /* slight change for C++, 2004/2/26 */
        void init_by_array(unsigned long init_key[],
                           int key_length)
        {
            int i, j, k;
            init_genrand(19650218UL);
            i=1; j=0;
            k = (MT_N>key_length ? MT_N : key_length);
            for (; k; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
                    + init_key[j] + j; /* non linear */
                mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
                i++; j++;
                if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
                if (j>=key_length) j=0;
            }
            for (k=MT_N-1; k; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
                    - i; /* non linear */
                mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
                i++;
                if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
            }
            mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
        }

        /* generates a random number on [0,0xffffffff]-interval */
        unsigned long genrand_int32(void)
        {
            unsigned long y;
            static unsigned long mag01[2]={0x0UL, MATRIX_A};
            /* mag01[x] = x * MATRIX_A  for x=0,1 */

            if (mti >= MT_N) { /* generate N words at one time */
                int kk;

                if (mti == MT_N+1)   /* if init_genrand() has not been called, */
                    init_genrand(5489UL); /* a default initial seed is used */
                
                for (kk=0;kk<MT_N-MT_M;kk++) {
                    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
                    mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
                }
                for (;kk<MT_N-1;kk++) {
                    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
                    mt[kk] = mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
                }
                y = (mt[MT_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
                mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
                
                mti = 0;
            }
  
            y = mt[mti++];

            /* Tempering */
            y ^= (y >> 11);
            y ^= (y << 7) & 0x9d2c5680UL;
            y ^= (y << 15) & 0xefc60000UL;
            y ^= (y >> 18);
            
            return y;
        }

    public:
/* initializes mt[MT_N] with a seed */
        static void init_genrand(unsigned long s)
        {
            unsigned long *tmp_mt = getInstance().mt;
            int & tmp_mti         = getInstance().mti;

            tmp_mt[0]= s & 0xffffffffUL;
            for (tmp_mti=1; tmp_mti<MT_N; tmp_mti++) {
                tmp_mt[tmp_mti] = 
                    (1812433253UL * (tmp_mt[tmp_mti-1] ^ (tmp_mt[tmp_mti-1] >> 30)) + tmp_mti); 
                tmp_mt[tmp_mti] &= 0xffffffffUL;
            }
        }

        /* generates a random number on [0,0x7fffffff]-interval */
        static long genrand_int31(void)
        {
            return (long)(getInstance().genrand_int32()>>1);
        }

        /* generates a random number on [0,1]-real-interval */
        static double genrand_real1(void)
        {
            return getInstance().genrand_int32()*(1.0/4294967295.0); 
            /* divided by 2^32-1 */ 
        }

        /* generates a random number on [0,1)-real-interval */
        static double genrand_real2(void)
        {
            return getInstance().genrand_int32()*(1.0/4294967296.0); 
            /* divided by 2^32 */
        }

        /* generates a random number on (0,1)-real-interval */
        static double genrand_real3(void)
        {
            return (((double)getInstance().genrand_int32()) + 0.5)*(1.0/4294967296.0); 
            /* divided by 2^32 */
        }

        /* generates a random number on [0,1) with 53-bit resolution*/
        static double genrand_res53(void) 
        { 
            unsigned long a=getInstance().genrand_int32()>>5, b=getInstance().genrand_int32()>>6; 
            return(a*67108864.0+b)*(1.0/9007199254740992.0); 
        } 
        /* These real versions are due to Isaku Wada, 2002/01/09 added */

        static void initialize(){
            getInstance().mti = MT_N + 1;
        }
    };




    class MTTS{
    private:
        enum{
            MT_N = 624,
            MT_M = 397,
        };
        enum{ 
            MATRIX_A = 0x9908b0dfUL,   /* constant vector a */
            UPPER_MASK = 0x80000000UL, /* most significant w-r bits */
            LOWER_MASK = 0x7fffffffUL, /* least significant r bits */
        };

        int mti; /* mti==MT_N+1 means mt[MT_N] is not initialized */
        unsigned long mt[MT_N]; /* the array for the state vector  */

#if 0
        MTTS(const MTTS &){};
        MTTS & operator=(const MTTS &);
#endif

/* initializes mt[MT_N] with a seed */
#if 0
        void init_genrand(unsigned long s)
        {
            mt[0]= s & 0xffffffffUL;
            for (mti=1; mti<MT_N; mti++) {
                mt[mti] = 
                    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
                /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
                /* In the previous versions, MSBs of the seed affect   */
                /* only MSBs of the array mt[].                        */
                /* 2002/01/09 modified by Makoto Matsumoto             */
                mt[mti] &= 0xffffffffUL;
                /* for >32 bit machines */
            }
        }
#endif

        /* initialize by an array with array-length */
        /* init_key is the array for initializing keys */
        /* key_length is its length */
        /* slight change for C++, 2004/2/26 */
        void init_by_array(unsigned long init_key[],
                           int key_length)
        {
            int i, j, k;
            init_genrand(19650218UL);
            i=1; j=0;
            k = (MT_N>key_length ? MT_N : key_length);
            for (; k; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
                    + init_key[j] + j; /* non linear */
                mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
                i++; j++;
                if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
                if (j>=key_length) j=0;
            }
            for (k=MT_N-1; k; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
                    - i; /* non linear */
                mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
                i++;
                if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
            }
            mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
        }

        /* generates a random number on [0,0xffffffff]-interval */
        unsigned long genrand_int32(void)
        {
            unsigned long y;
            static unsigned long mag01[2]={0x0UL, MATRIX_A};
            /* mag01[x] = x * MATRIX_A  for x=0,1 */

            if (mti >= MT_N) { /* generate N words at one time */
                int kk;

                if (mti == MT_N+1)   /* if init_genrand() has not been called, */
                    init_genrand(5489UL); /* a default initial seed is used */
                
                for (kk=0;kk<MT_N-MT_M;kk++) {
                    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
                    mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
                }
                for (;kk<MT_N-1;kk++) {
                    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
                    mt[kk] = mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
                }
                y = (mt[MT_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
                mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
                
                mti = 0;
            }
  
            y = mt[mti++];

            /* Tempering */
            y ^= (y >> 11);
            y ^= (y << 7) & 0x9d2c5680UL;
            y ^= (y << 15) & 0xefc60000UL;
            y ^= (y >> 18);
            
            return y;
        }

    public:
        MTTS() : mti(MT_N+1){}
        ~MTTS(){}
/* initializes mt[MT_N] with a seed */
        void init_genrand(unsigned long s)
        {
            unsigned long *tmp_mt = mt;
            int & tmp_mti         = mti;

            tmp_mt[0]= s & 0xffffffffUL;
            for (tmp_mti=1; tmp_mti<MT_N; tmp_mti++) {
                tmp_mt[tmp_mti] = 
                    (1812433253UL * (tmp_mt[tmp_mti-1] ^ (tmp_mt[tmp_mti-1] >> 30)) + tmp_mti); 
                tmp_mt[tmp_mti] &= 0xffffffffUL;
            }
        }

        /* generates a random number on [0,0x7fffffff]-interval */
        long genrand_int31(void)
        {
            return (long)(genrand_int32()>>1);
        }

        /* generates a random number on [0,1]-real-interval */
        double genrand_real1(void)
        {
            return genrand_int32()*(1.0/4294967295.0); 
            /* divided by 2^32-1 */ 
        }

        /* generates a random number on [0,1)-real-interval */
        double genrand_real2(void)
        {
            return genrand_int32()*(1.0/4294967296.0); 
            /* divided by 2^32 */
        }

        /* generates a random number on (0,1)-real-interval */
        double genrand_real3(void)
        {
            return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
            /* divided by 2^32 */
        }

        /* generates a random number on [0,1) with 53-bit resolution*/
        double genrand_res53(void) 
        { 
            unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
            return(a*67108864.0+b)*(1.0/9007199254740992.0); 
        } 
        /* These real versions are due to Isaku Wada, 2002/01/09 added */
    };

}
