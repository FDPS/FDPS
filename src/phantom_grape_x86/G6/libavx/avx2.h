#ifndef __AVX2_MACRO__
#define __AVX2_MACRO__

#define VFMADDPS(dst, reg1, reg2) asm("vfmadd231ps " reg1 "," reg2 "," dst);
// dst = dst + reg1 * reg2

#define VFNMADDPS(dst, reg1, reg2) asm("vfnmadd231ps " reg1 "," reg2 "," dst);
// dst = dst - reg1 * reg2

#define VFMSUB213PS_M(mem, reg, dst) \
    asm("vfmsub213ps %0, %"reg ", %"dst::"m"(mem));
// dst = reg * dst - mem                                                                             
#endif /* __AVX2_MACRO__ */
