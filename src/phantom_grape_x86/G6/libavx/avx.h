#define NVECS 8
#define NVECD 4

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

#define EQ   0
#define LT   1
#define LE   2
#define NEQ  4
#define NLT  5
#define NLE  6
#define NGE  9
#define NGT 10
#define GE  13
#define GT  14

#define VLOADPS(mem, reg) asm("vmovaps %0, %"reg::"m"(mem));
#define VSTORPS(reg, mem) asm("vmovaps %"reg ", %0" ::"m"(mem));
#define VLOADPD(mem, reg) asm("vmovapd %0, %"reg::"m"(mem));
#define VSTORPD(reg, mem) asm("vmovapd %"reg ", %0" ::"m"(mem));

#define VBROADCASTF128(mem, reg) asm("vbroadcastf128 %0, %"reg::"m"(mem));
#define VBROADCASTSS(mem, reg)   asm("vbroadcastss %0, %"reg::"m"(mem));
#define VBROADCASTSD(mem, reg)   asm("vbroadcastsd %0, %"reg::"m"(mem));

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

#define VPERM2F128(src1, src2, dest, imm) asm("vperm2f128 %0, %"src1 ", %"src2 ", %"dest " "::"g"(imm));
#define VEXTRACTF128(src, dest, imm)      asm("vextractf128 %0, %"src ", %"dest " "::"g"(imm));

#define VMERGE(src1, src2, dst)           asm("vperm2f128 %0, %"src1 ", %"src2 ", %"dst " "::"g"(2));
#define VUP2LOW(src, dst)                 asm("vextractf128 %0, %"src ", %"dst " "::"g"(1));

#define VANDPS(reg1, reg2, dst)      asm("vandps " reg1 "," reg2 "," dst);
#define VANDPS_M(mem, src, dst)      asm("vandps %0, %"src ", %"dst " "::"m"(mem));
#define VCMPPS(reg1, reg2, dst, imm) asm("vcmpps %0, %"reg1 ", %"reg2 ", %"dst " "::"g"(imm));
#define VMINPS(src1, src2, dst)      asm("vminps " src1 "," src2 "," dst);

#define VCMPNEQPS(reg1, reg2, dst) asm("vcmpps %0, %"reg1 ", %"reg2 ", %"dst " "::"g"(4));

#define PREFETCH(mem)        asm ("prefetcht0 %0"::"m"(mem))
