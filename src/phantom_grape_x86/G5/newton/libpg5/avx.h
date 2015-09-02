#ifndef __AVX_MACRO__
#define __AVX_MACRO__

#define YMM00 "%ymm0"
#define YMM01 "%ymm1"
#define YMM02 "%ymm2"
#define YMM03 "%ymm3"
#define YMM04 "%ymm4"
#define YMM05 "%ymm5"
#define YMM06 "%ymm6"
#define YMM07 "%ymm7"
#define YMM08 "%ymm8"
#define YMM09 "%ymm9"
#define YMM10 "%ymm10"
#define YMM11 "%ymm11"
#define YMM12 "%ymm12"
#define YMM13 "%ymm13"
#define YMM14 "%ymm14"
#define YMM15 "%ymm15"

#define VXORPS(reg1, reg2, dst) asm("vxorps " reg1 "," reg2 "," dst);
#define VXORPD(reg1, reg2, dst) asm("vxorpd " reg1 "," reg2 "," dst);
#define VMOVAPS(src, dst)    asm("vmovaps " src "," dst);

#define VLOADPS(mem, reg) asm("vmovaps %0, %"reg::"m"(mem));
#define VSTORPS(reg, mem) asm("vmovaps %"reg ", %0" ::"m"(mem));
#define VLOADPD(mem, reg) asm("vmovapd %0, %"reg::"m"(mem));
#define VSTORPD(reg, mem) asm("vmovapd %"reg ", %0" ::"m"(mem));

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

#define VSQRTPS(reg, dst) asm("vsqrtps " reg "," dst);

#define VZEROALL asm("vzeroall");

#define VCVTPD2PS(reg, dst) asm("vcvtpd2ps " reg "," dst);
#define VCVTPS2PD(reg, dst) asm("vcvtps2pd " reg "," dst);

#define VPERM2F128(src1, src2, dest, imm) asm("vperm2f128 %0, %"src2 ", %"src1 ", %"dest " "::"g"(imm));
#define VEXTRACTF128(src, dest, imm)      asm("vextractf128 %0, %"src ", %"dest " "::"g"(imm));
#define VSHUFPS(src1, src2, dest, imm) \
  asm("vshufps %0, %"src2 ", %"src1 ", %"dest " "::"g"(imm));

#define VANDPS(reg1, reg2, dst)      asm("vandps " reg1 "," reg2 "," dst);
#define VCMPPS(reg1, reg2, dst, imm) asm("vcmpps %0, %"reg1 ", %"reg2 ", %"dst " "::"g"(imm));

#define PREFETCH(mem)        asm ("prefetcht0 %0"::"m"(mem))

#endif /* __AVX_MACRO__ */
