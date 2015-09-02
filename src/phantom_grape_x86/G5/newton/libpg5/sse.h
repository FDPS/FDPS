#ifndef __SSE_MACRO__
#define __SSE_MACRO__

#define XMM00 "%xmm0"
#define XMM01 "%xmm1"
#define XMM02 "%xmm2"
#define XMM03 "%xmm3"
#define XMM04 "%xmm4"
#define XMM05 "%xmm5"
#define XMM06 "%xmm6"
#define XMM07 "%xmm7"
#define XMM08 "%xmm8"
#define XMM09 "%xmm9"
#define XMM10 "%xmm10"
#define XMM11 "%xmm11"
#define XMM12 "%xmm12"
#define XMM13 "%xmm13"
#define XMM14 "%xmm14"
#define XMM15 "%xmm15"

#define XORPS(a, b)         asm("xorps "  a  ","  b );
#define LOADPS(mem, reg)    asm("movaps %0, %"reg::"m"(mem));
#define LOADLPS(mem, reg)   asm("movlps %0, %"reg::"m"(mem));
#define LOADHPS(mem, reg)   asm("movhps %0, %"reg::"m"(mem));
#define LDDQU(mem, reg)     asm("lddqu %0, %"reg::"m"(mem));
#define MOVAPSX(mem, reg)   asm("movaps %0, %"reg::"x"(mem));
#define STORPS(reg, mem)    asm("movaps %"reg " , %0"::"m"(mem));
#define MOVAPS(src, dst)    asm("movaps " src "," dst);
#define MOVQ(src, dst)      asm("movq " src "," dst);
#define BCAST0(reg)         asm("shufps $0x00, " reg ","  reg);
#define BCAST1(reg)         asm("shufps $0x55, " reg ","  reg);
#define BCAST2(reg)         asm("shufps $0xaa, " reg ","  reg);
#define BCAST3(reg)         asm("shufps $0xff, " reg ","  reg);
#define MULPS(src, dst)     asm("mulps " src "," dst);
#define MULPS_M(mem, reg)   asm("mulps %0, %"reg::"m"(mem));
#define ADDPS(src, dst)     asm("addps " src ","  dst);
#define ADDPS_M(mem, reg)   asm("addps %0, %"reg::"m"(mem));
#define SUBPS(src, dst)     asm("subps "  src "," dst);
#define SUBPS_M(mem, reg)   asm("subps %0, %"reg::"m"(mem));
#define MINPS(src,dst)      asm("minps  " src "," dst);
#define MINPS_M(mem, reg)   asm("minps %0, %"reg ::"m"(mem));
#define RSQRTPS(src, dst)   asm("rsqrtps " src "," dst);
#define RSQRTPS_M(mem, reg) asm("rsqrtps %0, %"reg ::"m"(mem));
#define SQRTPS(src, dst)    asm("sqrtps " src "," dst);
#define SQRTPS_M(mem, reg)  asm("sqrtps %0, %"reg ::"m"(mem));
#define RCPPS(src, dst)     asm("rcpps " src "," dst);
#define RCPPS_M(mem, reg)   asm("rcpps %0, %"reg ::"m"(mem));
#define CVTTPS2PI(src, dst) asm("cvttps2pi " src "," dst);
#define CVTTPS2DQ(src, dst) asm("cvttps2dq " src "," dst);
#define CVTDQ2PS(src, dst)  asm("cvtdq2ps " src "," dst);
#define UNPCKLPS(src, dst)  asm("unpcklps "  src "," dst);
#define UNPCKHPS(src, dst)  asm("unpckhps "  src "," dst);
#define CMPLTPS(src, dst)   asm("cmpltps "  src "," dst);
#define CMPLTPS_M(mem, reg) asm("cmpltps %0, %"reg ::"m"(mem));
#define ANDPS(src, dst)     asm("andps " src "," dst);
#define ANDNPS(src, dst)    asm("andnps " src "," dst);
#define EXTINT0(reg, ireg)  asm volatile("pextrw $0, %"reg " , %0":"=r"(ireg));
#define EXTINT1(reg, ireg)  asm volatile("pextrw $2, %"reg " , %0":"=r"(ireg));
#define EXTINT2(reg, ireg)  asm volatile("pextrw $4, %"reg " , %0":"=r"(ireg));
#define EXTINT3(reg, ireg)  asm volatile("pextrw $6, %"reg " , %0":"=r"(ireg));
#define PREFETCH(mem)       asm ("prefetcht0 %0"::"m"(mem))
#define NOP                 asm("nop");

#endif /* __SSE_MACRO__ */
