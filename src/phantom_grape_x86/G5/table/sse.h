#define X2_X    "%xmm0"
#define Y2_X    "%xmm1"
#define Z2_X    "%xmm2"
#define R2CUT_X "%xmm3"
#define BUF0_X  "%xmm4"
#define BUF1_X  "%xmm5"
#define BUF2_X  "%xmm6"
#define BUF3_X  "%xmm7"

#define TWO_X   "%xmm8"
#define MJ_X    "%xmm9"
#define DX_X    "%xmm10"
#define DY_X    "%xmm11"
#define DZ_X    "%xmm12"
#define AX_X    "%xmm13"
#define AY_X    "%xmm14"
#define AZ_X    "%xmm15"

// instructions
#define XORPS(a, b)       asm("xorps "  a  ","  b );
#define CLEARPS(a) XORPS(a,a)
#define LOADPS(mem, reg)  asm("movaps %0, %"reg::"m"(mem));
#define LOADLPS(mem, reg)  asm("movlps %0, %"reg::"m"(mem));
#define MOVAPSX(mem, reg)  asm("movaps %0, %"reg::"x"(mem));
#define STORPS(reg, mem)  asm("movaps %"reg " , %0"::"m"(mem));
#define MOVAPS(src, dst)  asm("movaps " src ", " dst);
#define MOVQ(src, dst)    asm("movq " src ", " dst);
#define BCAST0(reg)       asm("shufps $0x00, " reg ","  reg);
#define BCAST1(reg)       asm("shufps $0x55, " reg ","  reg);
#define BCAST2(reg)       asm("shufps $0xaa, " reg ","  reg);
#define BCAST3(reg)       asm("shufps $0xff, " reg ","  reg);
#define MULPS(src, dst)   asm("mulps " src ", " dst);
#define MULPS_M(mem, reg) asm("mulps %0, %"reg::"m"(mem));
#define ADDPS(src, dst)   asm("addps " src ", "  dst);
#define SUBPS(src, dst)   asm("subps "  src ", " dst);
#define SUBPS_M(mem, reg) asm("subps %0, %"reg::"m"(mem));
#define MINPS(src, dst)   asm("minps "  src ", " dst);
#define UNPCKLPS(src, dst)   asm("unpcklps "  src ", " dst);
#define UNPCKHPS(src, dst)   asm("unpckhps "  src ", " dst);
#define PSRLD(imm, reg)   asm("psrld %0, %" reg::"I"(imm));
#define PSLLD(imm, reg)   asm("pslld %0, %" reg::"I"(imm));
// #define CVTTPS2DQ(src, dst)   asm("cvttps2dq "  src "," dst);
// #define CVTDQ2PS(src, dst)   asm("cvtdq2ps "  src "," dst);
#define NOP               asm("nop");
#define EXTINT0(reg, ireg) asm volatile("pextrw $0, %"reg ", %0":"=r"(ireg));
#define EXTINT1(reg, ireg) asm volatile("pextrw $2, %"reg ", %0":"=r"(ireg));
#define EXTINT2(reg, ireg) asm volatile("pextrw $4, %"reg ", %0":"=r"(ireg));
#define EXTINT3(reg, ireg) asm volatile("pextrw $6, %"reg ", %0":"=r"(ireg));
// #define RSQRTPS(src, dst) asm("rsqrtps " src "," dst);
// #define MOVHLPS(src, dst) asm("movhlps " src "," dst);
#define PREFETCH(mem) asm ("prefetcht0 %0"::"m"(mem))

#define LOADDQA(mem, reg) asm("movdqa %0, %"reg::"m"(mem));
#define STORDQA(src, mem) asm("movdqa %"src " , %0"::"m"(mem));
