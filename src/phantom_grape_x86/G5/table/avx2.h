#ifndef __AVX2_MACRO__
#define __AVX2_MACRO__

#define VLOADDQA(mem, dst) asm("vmovdqa %0, %"dst::"m"(mem));
#define VSTORDQA(src, mem) asm("vmovdqa %"src ", %0"::"m"(mem));

#define VFMADDPS(dst, reg1, reg2) asm("vfmadd231ps " reg1 "," reg2 "," dst);
// dst = dst + reg1 * reg2

#define VFNMADDPS(dst, reg1, reg2) asm("vfnmadd231ps " reg1 "," reg2 "," dst);
// dst = dst - reg1 * reg2

#define VPGATHERF32ARRAY2D(mem, idx, msk, dst) \
    asm("vpgatherdd %" msk", (%[ptr], %" idx", 8), %"dst:: [ptr] "r" (mem) : "memory");

#endif /* __AVX2_MACRO__ */

