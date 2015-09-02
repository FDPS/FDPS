#ifndef __AVX2_MACRO__
#define __AVX2_MACRO__

#define VFMADDPS(dst, reg1, reg2) asm("vfmadd231ps " reg1 "," reg2 "," dst);
// dst = dst + reg1 * reg2

#define VFNMADDPS(dst, reg1, reg2) asm("vfnmadd231ps " reg1 "," reg2 "," dst);
// dst = dst - reg1 * reg2

#endif /* __AVX2_MACRO__ */

