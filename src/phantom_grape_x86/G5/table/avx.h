#ifndef __AVX_MACRO__
#define __AVX_MACRO__

#define X2    "%ymm0"
#define Y2    "%ymm1"
#define Z2    "%ymm2"
#define R2CUT "%ymm3"  // fixed
#define BUF0  "%ymm4"
#define BUF1  "%ymm5"
#define BUF2  "%ymm6"
#define BUF3  "%ymm7"  // fixed
#define TWO   "%ymm8"  // fixed
#define MJ    "%ymm9"
#define DX    "%ymm10"
#define DY    "%ymm11"
#define DZ    "%ymm12"
#define AX    "%ymm13" // fixed
#define AY    "%ymm14" // fixed
#define AZ    "%ymm15" // fixed

#define VXORPS(reg1, reg2, dst) asm("vxorps " reg1 "," reg2 "," dst);
#define VXORPD(reg1, reg2, dst) asm("vxorpd " reg1 "," reg2 "," dst);
#define VMOVAPS(src, dst)    asm("vmovaps " src "," dst);

#define VLOADPS(mem, reg) asm("vmovaps %0, %"reg::"m"(mem));
#define VSTORPS(reg, mem) asm("vmovaps %"reg ", %0" ::"m"(mem));
#define VLOADPD(mem, reg) asm("vmovapd %0, %"reg::"m"(mem));
#define VSTORPD(reg, mem) asm("vmovapd %"reg ", %0" ::"m"(mem));
#define VLOADDQA(mem, dst) asm("vmovdqa %0, %"dst::"m"(mem));
#define VSTORDQA(src, mem) asm("vmovdqa %"src ", %0"::"m"(mem));
#define VLOADDQU(mem, dst) asm("vmovdqu %0, %"dst::"m"(mem));
#define VSTORDQU(src, mem) asm("vmovdqu %"src ", %0"::"m"(mem));

#define VLOADLPS(mem, src1, src2) asm("vmovlps %0, %"src1 ", %"src2::"m"(mem));
#define VLOADHPS(mem, src1, src2) asm("vmovhps %0, %"src1 ", %"src2::"m"(mem));

#define VBROADCASTF128(mem, reg) asm("vbroadcastf128 %0, %"reg::"m"(mem));

#define VADDPS(reg1, reg2, dst) asm("vaddps " reg1 "," reg2 "," dst);
#define VADDPS_M(mem, reg, dst) asm("vaddps %0, %"reg ", %"dst " "::"m"(mem));
#define VADDPD(reg1, reg2, dst) asm("vaddpd " reg1 "," reg2 "," dst);
#define VADDPD_M(mem, reg, dst) asm("vaddpd %0, %"reg ", %"dst " "::"m"(mem));

#define VSUBPS(reg1, reg2, dst) asm("vsubps " reg1 "," reg2 "," dst);
#define VSUBPS_M(mem, reg, dst) asm("vsubps %0, %"reg ", %"dst " "::"m"(mem));
#define VSUBPD(reg1, reg2, dst) asm("vsubpd " reg1 "," reg2 "," dst);
#define VSUBPD_M(mem, reg, dst) asm("vsubpd %0, %"reg ", %"dst " "::"m"(mem));

#define VMULPS(reg1, reg2, dst) asm("vmulps " reg1 "," reg2 "," dst);
#define VMULPS_M(mem, reg, dst) asm("vmulps %0, %"reg ", %"dst " "::"m"(mem));
#define VMULPD(reg1, reg2, dst) asm("vmulpd " reg1 "," reg2 "," dst);
#define VMULPD_M(mem, reg, dst) asm("vmulpd %0, %"reg ", %"dst " "::"m"(mem));

#define VDIVPS(reg1, reg2, dst) asm("vdivps " reg1 "," reg2 "," dst);
#define VDIVPS_M(mem, reg, dst) asm("vdivps %0, %"reg ", %"dst " "::"m"(mem));
#define VDIVPD(reg1, reg2, dst) asm("vdivpd " reg1 "," reg2 "," dst);
#define VDIVPD_M(mem, reg, dst) asm("vdivpd %0, %"reg ", %"dst " "::"m"(mem));

#define VHADDPD(src1, src2, dst) asm("vhaddpd " src1 "," src2 "," dst);
// src2_1 + src2_2 -> dst_1
// src1_1 + src1_2 -> dst_2
// src2_3 + src2_4 -> dst_3
// src1_3 + src1_4 -> dst_4

#define VRSQRTPS(reg, dst) asm("vrsqrtps " reg "," dst);

#define VZEROALL asm("vzeroall");
#define VZEROUPPER asm("vzeroupper");

#define VCVTPD2PS(reg, dst) asm("vcvtpd2ps " reg "," dst);
#define VCVTPS2PD(reg, dst) asm("vcvtps2pd " reg "," dst);
#define VCVTDQ2PS(src, dst) asm("vcvtdq2ps " src "," dst);
#define VCVTPS2DQ(src, dst) asm("vcvtps2dq " src "," dst);

#define VPERM2F128(src1, src2, dest, imm) asm("vperm2f128 %0, %"src2 ", %"src1 ", %"dest " "::"g"(imm));
#define VPERM2F128_M02(src, mem, dst) asm("vperm2f128 $0x02, %0, %"src ", %"dst::"m"(mem));
// "mem" is second operand, and "src" is first operand!

#define VEXTRACTF128(src, dest, imm)      asm("vextractf128 %0, %"src ", %"dest " "::"g"(imm));
#define VSHUFPS(src1, src2, dest, imm) \
  asm("vshufps %0, %"src2 ", %"src1 ", %"dest " "::"g"(imm));

#define VANDPS(reg1, reg2, dst)      asm("vandps " reg1 "," reg2 "," dst);
#define VANDPS_M(mem, src, dst)      asm("vandps %0, %"src ", %"dst::"m"(mem));
#define VPAND(src1, src2, dst) asm("vpand " src1 "," src2 "," dst);
#define VPAND_M(mem, src, dst) asm("vpand %0, %"src ", %"dst::"m"(mem));

#define VCMPPS(reg1, reg2, dst, imm) asm("vcmpps %0, %"reg1 ", %"reg2 ", %"dst " "::"g"(imm));
#define VMINPS(src1, src2, dst)    asm("vminps "  src1 ", " src2 "," dst);
#define VMINPS_M(mem, src, dst)    asm("vminps %0, %"src ", %"dst::"m"(mem));

#define VUNPCKLPS(src1, src2, dst) asm("vunpcklps " src1 "," src2 "," dst);
#define VUNPCKHPS(src1, src2, dst) asm("vunpckhps " src1 "," src2 "," dst);

#define VPSRLD(imm, src1, src2) asm("vpsrld %0, %"src1 ", %"src2::"I"(imm));
#define VPSLLD(imm, src1, src2) asm("vpslld %0, %"src1 ", %"src2::"I"(imm));

#define VEXTINT0(reg, ireg) asm volatile("vpextrw $0, %"reg ", %0":"=r"(ireg));
#define VEXTINT1(reg, ireg) asm volatile("vpextrw $2, %"reg ", %0":"=r"(ireg));
#define VEXTINT2(reg, ireg) asm volatile("vpextrw $4, %"reg ", %0":"=r"(ireg));
#define VEXTINT3(reg, ireg) asm volatile("vpextrw $6, %"reg ", %0":"=r"(ireg));

#define MEMORYFLUSH asm("":::"memory");

#endif /* __AVX_MACRO__ */

