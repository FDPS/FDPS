require_relative "common.rb"

def get_declare_type_avx2(type)
  case type
  when "F64" then
    decl = "__m256d"
  when "F32" then
    decl = "__m256"
  when "F64vec" then
    decl = "__m256dx3"
  when "F32vec" then
    decl = "__m256x3"
  when "F64vec3" then
    decl = "__m256dx3"
  when "F32vec3" then
    decl = "__m256x3"
  when "F64vec4" then
    decl = "__m256dx4"
  when "F32vec4" then
    decl = "__m256x4"
  when "S64" then
    decl = "__m256i"
  when "S32" then
    decl = "__m256i"
  when "S16" then
    decl = "__m256i"
  when "U64" then
    decl = "__m256i"
  when "U32" then
    decl = "__m256i"
  when "U16" then
    decl = "__m256i"
  when "B64" then
    decl = "__m256"
  when "B32" then
    decl = "__m256"
  when "B16" then
    decl = "__m256"
  else
    abort "error: unsupported declaration of type #{type} for AVX2"
  end
  decl
end

def get_type_suffix_avx2(type)
  case type
  when "F64" then
    suffix = "pd"
  when "F32" then
    suffix = "ps"
  when "F64vec" then
    suffix = "pdx3"
  when "F32vec" then
    suffix = "psx3"
  when "F64vec3" then
    suffix = "pdx3"
  when "F32vec3" then
    suffix = "psx3"
  when "F64vec4" then
    suffix = "pdx4"
  when "F32vec4" then
    suffix = "psx4"
  when "S64" then
    suffix = "epi64"
  when "S32" then
    suffix = "epi32"
  when "S16" then
    suffix = "epi16"
  when "U64" then
    suffix = "epu64"
  when "U32" then
    suffix = "epu32"
  when "U16" then
    suffix = "epu16"
  when "B64" then
    suffix = "mask"
  when "B32" then
    suffix = "mask"
  when "B16" then
    suffix = "mask"
  else
    abort "error: unsupported suffix of type #{type} for AVX2"
  end
  suffix
end

class Kernelprogram
  def reserved_func_def_avx2(conversion_type)
    code = ""
    code += "__m256 rsqrt(__m256 op){\n"
    # code += "  __m256 rinv = _mm256_rsqrt14_ps(op);\n"
    # code += "  __m256 h = _mm256_mul_ps(op,rinv);\n"
    # code += "  h = _mm256_fnmadd_ps(h,rinv,_mm256_set1_ps(1.f));\n"
    # code += "  __m256 poly = _mm256_fmadd_ps(h,_mm256_set1_ps(0.375f),_mm256_set1_ps(0.5f));\n"
    # code += "  poly = _mm256_mul_ps(poly,h);\n"
    # code += "  rinv = _mm256_fmadd_ps(rinv,poly,rinv);\n"
    # code += "  return rinv;\n"
    code += "  return _mm256_rsqrt_ps(op);\n"
    code += "}\n"

    code += "__m256 sqrt(__m256 op){ return _mm256_sqrt_ps(op); }\n"

    # http://math-koshimizu.hatenablog.jp/entry/2017/07/28/083000
    code += "__m256 inv(__m256 op){\n"
    # code += "  __m256 x1 = _mm256_rcp14_ps(op);\n"
    # code += "  __m256 x2 = _mm256_fnmadd_ps(op,x1,_mm256_set1_ps(2.f));\n"
    # code += "  x2 = _mm256_mul_ps(x2,x1);\n"
    # code += "  __m256 ret = _mm256_fnmadd_ps(op,x2,_mm256_set1_ps(2.f));\n"
    # code += "  ret = _mm256_mul_ps(ret,x2);\n"
    # code += "  return ret;\n"
    code += "return _mm256_rcp_ps(op);\n"
    code += "}\n"

    code += "__m256d rsqrt(__m256d op){\n"

    # impl using float rsqrt
    #code += "  __m256 x = _mm256_castps128_ps256(_mm256_cvtpd_ps(op));\n"
    #code += "  x = rsqrt(x);\n"
    #code += "  return _mm256_cvtps_pd(_mm256_castps256_ps128(x));\n"

    # https://en.wikipedia.org/wiki/Fast_inverse_square_root
    code += "  __m256d y = _mm256_castsi256_pd(_mm256_sub_epi64(_mm256_set1_epi64x(0x5fe6eb50c7b537a9LL),_mm256_srlv_epi64(_mm256_castpd_si256(op),_mm256_set1_epi64x(1))));\n"
    code += "  __m256d h = _mm256_mul_pd(op,y);\n"
    code += "  h = _mm256_fnmadd_pd(h,y,_mm256_set1_pd(1.0));\n"
    code += "  __m256d poly = _mm256_fmadd_pd(h,_mm256_set1_pd(0.375),_mm256_set1_pd(0.5));\n"
    code += "  poly = _mm256_mul_pd(poly,h);\n"
    code += "  y = _mm256_fmadd_pd(y,poly,y);\n"
    code += "  return y;\n"

    code += "}"


    code += "__m256d sqrt(__m256d op){\n"
    code += "  return _mm256_sqrt_pd(op);\n"
    code += "}"
    code += "__m256d inv(__m256d op){\n"
    code += "  __m256 x = _mm256_castps128_ps256(_mm256_cvtpd_ps(op));\n"
    code += "  x = inv(x);\n"
    code += "  return _mm256_cvtps_pd(_mm256_castps256_ps128(x));\n"
    code += "}"

    code += "__m256d max(__m256d a,__m256d b){ return _mm256_max_pd(a,b);}\n"
    code += "__m256d min(__m256d a,__m256d b){ return _mm256_min_pd(a,b);}\n"
    code += "__m256  max(__m256  a,__m256  b){ return _mm256_max_ps(a,b);}\n"
    code += "__m256  min(__m256  a,__m256  b){ return _mm256_min_ps(a,b);}\n"
    code += "__m256i max(__m256i a,__m256i b){ return _mm256_max_epi32(a,b);}\n"
    code += "__m256i min(__m256i a,__m256i b){ return _mm256_min_epi32(a,b);}\n"

    code += "__m256d table(__m256d tab,__m256i index){ return _mm256_permutexvar_pd(index,tab);}\n"
    code += "__m256  table(__m256  tab,__m256i index){ return _mm256_permutexvar_ps(index,tab);}\n"

    #code += "__m256d to_double(__m256i op){ return _mm256_cvt_roundepi64_pd(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    code += "__m256  to_float(__m256i op){ return _mm256_cvtepi32_ps(op);}\n"
    #code += "__m256d to_float(__m256i op){ return _mm256_cvt_roundepu64_pd(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    #code += "__m256  to_float(__m256i op){ return _mm256_cvt_roundepu32_ps(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    #code += "__m256i  to_long(__m256d op){ return _mm256_cvt_roundpd_epi64(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    code += "__m256i  to_int(__m256  op){ return _mm256_cvtps_epi32(op);}\n"
    #code += "__m256i  to_ulong(__m256d op){ return _mm256_cvt_roundpd_epu64(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    #code += "__m256i  to_uint(__m256  op){ return _mm256_cvt_roundps_epu32(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"

    code += "void transpose2x2_pd(__m256d& a, __m256d& b){\n"
    code += "  __m256d tmp = _mm256_unpacklo_pd(a,b);\n"
    code += "  b = _mm256_unpackhi_pd(a,b);\n"
    code += "  a = tmp;\n"
    code += "}\n"

    code += "void transpose4x4_pd(__m256d& a, __m256d& b,__m256d& c,__m256d& d){\n"
    code += "  __m256d tmp0 = _mm256_unpacklo_pd(a, b);\n"
    code += "  __m256d tmp1 = _mm256_unpackhi_pd(a, b);\n"
    code += "  __m256d tmp2 = _mm256_unpacklo_pd(c, d);\n"
    code += "  __m256d tmp3 = _mm256_unpackhi_pd(c, d);\n"
    code += "  a = _mm256_permute2f128_pd(tmp0, tmp2, 0|(2<<4));\n"
    code += "  b = _mm256_permute2f128_pd(tmp1, tmp3, 0|(2<<4));\n"
    code += "  c = _mm256_permute2f128_pd(tmp0, tmp2, 1|(3<<4));\n"
    code += "  d = _mm256_permute2f128_pd(tmp1, tmp3, 1|(3<<4));\n"
    code += "}\n"

    code += "void unpack2x2_pd(__m256d& a,__m256d& b){\n"
    code += "  transpose2x2_pd(a,b);\n"
    code += "}\n"
    code += "void pack2x2_pd(__m256d& a,__m256d& b){\n"
    code += "  transpose2x2_pd(a,b);\n"
    code += "}\n"
    code += "void unpack4x4_pd(__m256d& a,__m256d& b,__m256d& c,__m256d& d){\n"
    code += "  transpose4x4_pd(a,b,c,d);\n"
    code += "}\n"
    code += "void pack4x4_pd(__m256d& a,__m256d& b,__m256d& c,__m256d& d){\n"
    code += "  transpose4x4_pd(a,b,c,d);\n"
    code += "}\n"

    code += "void unpack2x2_ps(__m256& a, __m256& b){\n"
    code += "  __m256 tmp = _mm256_shuffle_ps(a,b,0xdd);\n"
    code += "  b = _mm256_shuffle_ps(a,b,0x88);\n"
    code += "  a = tmp;\n"
    code += "}\n"
    code += "void pack2x2_ps(__m256& a, __m256& b){\n"
    code += "  __m256 tmp = _mm256_unpackhi_ps(a,b);\n"
    code += "  b = _mm256_unpacklo_ps(a,b);\n"
    code += "  a = tmp;\n"
    code += "}\n"

    #code += "void unpack4x4_ps(__m256& a, __m256& b,__m256& c,__m256& d){\n"
    #code += "  __m256 src0 = _mm256_permute2f128_ps(a,b,0x20);\n"
    #code += "  __m256 src1 = _mm256_permute2f128_ps(a,b,0x31);\n"
    #code += "  __m256 src3 = _mm256_permute2f128_ps(c,d,0x20);\n"
    #code += "  __m256 src4 = _mm256_permute2f128_ps(c,d,0x31);\n"
    #code += "  __m256d tmp0 = _mm256_castps_pd(_mm256_unpacklo_ps(a, b)); // x0 x2 y0 y2 x1 y1 x3 y3\n"
    #code += "  __m256d tmp1 = _mm256_castps_pd(_mm256_unpackhi_ps(a, b)); // z0 z2 w0 w2 z1 w1 z3 w3\n"
    #code += "  __m256d tmp2 = _mm256_castps_pd(_mm256_unpacklo_ps(c, d)); // x4 x6 y4 y6 x5 x7 y5 y7\n"
    #code += "  __m256d tmp3 = _mm256_castps_pd(_mm256_unpackhi_ps(c, d)); // z4 z6 w4 w6 z5 z7 w5 w7\n"
    #code += "  a = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp0, tmp2)); // x0 x2 x4 x6 x1 x3 x5 x7\n"
    #code += "  b = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp0, tmp2)); // y0 y2 y4 y6 y1 y3 y5 y7\n"
    #code += "  c = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp1, tmp3)); // z0 z2 z4 z6 z1 z3 z5 z7\n"
    #code += "  d = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp1, tmp3)); // w0 w2 w4 w6 w1 w3 w5 w7\n"
    #code += "}\n"
    code += "void gather4_ps(__m256& a, __m256& b,__m256& c,__m256& d){\n"
    code += "  __m256 src0 = _mm256_permute2f128_ps(a,c,0x20); // x0 y0 z0 w0 x4 y4 z4 w4\n"
    code += "  __m256 src1 = _mm256_permute2f128_ps(a,c,0x31); // x1 y1 z1 w1 x5 y5 z5 w5\n"
    code += "  __m256 src2 = _mm256_permute2f128_ps(b,d,0x20);\n"
    code += "  __m256 src3 = _mm256_permute2f128_ps(b,d,0x31);\n"
    code += "  __m256d tmp0 = _mm256_castps_pd(_mm256_unpacklo_ps(src0, src1)); // x0 x2 y0 y2 x1 y1 x3 y3\n"
    code += "  __m256d tmp1 = _mm256_castps_pd(_mm256_unpackhi_ps(src0, src1)); // z0 z2 w0 w2 z1 w1 z3 w3\n"
    code += "  __m256d tmp2 = _mm256_castps_pd(_mm256_unpacklo_ps(src2, src3)); // x4 x6 y4 y6 x5 x7 y5 y7\n"
    code += "  __m256d tmp3 = _mm256_castps_pd(_mm256_unpackhi_ps(src2, src3)); // z4 z6 w4 w6 z5 z7 w5 w7\n"
    code += "  a = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp0, tmp2)); // x0 x2 x4 x6 x1 x3 x5 x7\n"
    code += "  b = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp0, tmp2)); // y0 y2 y4 y6 y1 y3 y5 y7\n"
    code += "  c = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp1, tmp3)); // z0 z2 z4 z6 z1 z3 z5 z7\n"
    code += "  d = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp1, tmp3)); // w0 w2 w4 w6 w1 w3 w5 w7\n"
    code += "}\n"
    code += "void scatter4_ps(__m256& a, __m256& b,__m256& c,__m256& d){\n"
    code += "  __m256d tmp0 = _mm256_castps_pd(_mm256_unpacklo_ps(a, b)); // x0 x2 y0 y2 x1 y1 x3 y3\n"
    code += "  __m256d tmp1 = _mm256_castps_pd(_mm256_unpackhi_ps(a, b)); // z0 z2 w0 w2 z1 w1 z3 w3\n"
    code += "  __m256d tmp2 = _mm256_castps_pd(_mm256_unpacklo_ps(c, d)); // x4 x6 y4 y6 x5 x7 y5 y7\n"
    code += "  __m256d tmp3 = _mm256_castps_pd(_mm256_unpackhi_ps(c, d)); // z4 z6 w4 w6 z5 z7 w5 w7\n"
    code += "  __m256 dst0 = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp0, tmp2)); // x0 x2 x4 x6 x1 x3 x5 x7\n"
    code += "  __m256 dst1 = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp0, tmp2)); // y0 y2 y4 y6 y1 y3 y5 y7\n"
    code += "  __m256 dst2 = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp1, tmp3)); // z0 z2 z4 z6 z1 z3 z5 z7\n"
    code += "  __m256 dst3 = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp1, tmp3)); // w0 w2 w4 w6 w1 w3 w5 w7\n"
    code += "  a = _mm256_permute2f128_ps(dst0,dst1,0x20);\n"
    code += "  b = _mm256_permute2f128_ps(dst2,dst3,0x20);\n"
    code += "  c = _mm256_permute2f128_ps(dst0,dst1,0x31);\n"
    code += "  d = _mm256_permute2f128_ps(dst2,dst3,0x31);\n"
    code += "}\n"
    code += "void unpack4x4_ps(__m256& a, __m256& b,__m256& c,__m256& d){\n"
    code += "  gather4_ps(a,b,c,d);\n"
    code += "}\n"
    code += "void pack4x4_ps(__m256& a, __m256& b,__m256& c,__m256& d){\n"
    code += "  scatter4_ps(a,b,c,d);\n"
    code += "}\n"
    code += "void unpack4x4_ps(__m256x4& v){\n"
    code += "  unpack4x4_ps(v.v0,v.v1,v.v2,v.v3);\n"
    code += "}\n"
    code += "void pack4x4_ps(__m256x4& v){\n"
    code += "  unpack4x4_ps(v.v0,v.v1,v.v2,v.v3);\n"
    code += "}\n"
    code += "void unpack4x4_pd(__m256dx4& v){\n"
    code += "  unpack4x4_pd(v.v0,v.v1,v.v2,v.v3);\n"
    code += "}\n"
    code += "void pack4x4_pd(__m256dx4& v){\n"
    code += "  unpack4x4_pd(v.v0,v.v1,v.v2,v.v3);\n"
    code += "}\n"

    code += "void print_ps(const __m256 v){ for(int i=0;i<8;i++) printf(\" %f\",v[i]); printf(\"\\n\");}\n"
    code += "void print_ps(const __m256x4 v){\n"
    code += "  print_ps(v.v0);\n"
    code += "  print_ps(v.v1);\n"
    code += "  print_ps(v.v2);\n"
    code += "  print_ps(v.v3);\n"
    code += "}\n"
  end
end

class FuncCall
  def convert_to_code_avx2(conversion_type)
    retval = @name
    retval += "("
    if @ops.length > 0
      @ops.each_with_index{ |op,i|
        retval += "," if i > 0
        retval += op.convert_to_code(conversion_type)
      }
    end
    retval += ")"
    retval
  end
end

class Expression
  def convert_to_code_avx2(conversion_type)
    if [:uminus,:plus,:minus,:mult,:div,:and,:or,:andnot,:dot].index(@operator)
      case self.get_type
      when "F64"
        type = "pd"
      when "F32"
        type = "ps"
      when "F64vec"
        type = "pdx3"
      when "F64vec3"
        type = "pdx3"
      when "F32vec"
        type = "psx3"
      when "F32vec3"
        type = "psx3"
      when "F64vec2"
        type = "pdx2"
      when "F32vec2"
        type = "psx2"
      when "F64vec4"
        type = "pdx4"
      when "F32vec4"
        type = "psx4"
      when "S64"
        type = "epi64"
      when "S32"
        type = "epi32"
      when "S16"
        type = "epi16"
      when "U64"
        type = "epu64"
      when "U32"
        type = "epu32"
      when "U16"
        type = "epu16"
      else
        abort "unsupported return type of operator #{@operator} for AVX2"
      end
    elsif [:lt,:le,:gt,:ge,:eq,:neq,:land,:lor,:not,:landnot].index(@operator)
      case @lop.get_type
      when "F64"
        type = "pd"
      when "F32"
        type = "ps"
      when "F64vec"
        type = "pdx3"
      when "F32vec"
        type = "psx3"
      when "S64"
        type = "pd"
      when "S32"
        type = "ps"
      when "U64"
        type = "pd"
      when "U32"
        type = "ps"
      when /B/
        type = "ps"
      else
        abort "unsupported mask type for AVX2"
      end
    end

    retval = "_mm256"
    set1_suffix = ""
    set1_suffix = "x" if type == "epi64"
    lop = @lop.convert_to_code(conversion_type)
    rop = @rop.convert_to_code(conversion_type) if @rop != nil
    if [:lt,:le,:gt,:ge,:eq,:neq,:land,:lor,:landnot].index(@operator) && @lop.get_type =~ /(S|U)(64|32)/
      lop = "_mm256_castsi256_#{type}(" + lop + ")"
      rop = "_mm256_castsi256_#{type}(" + rop + ")" if @rop != nil
    end
    case @operator
    when :uminus then
      retval += "_sub_#{type}("
      retval += "_mm256_set1_#{type}#{set1_suffix}((PIKG::#{@type})0.0)," + lop + ")"
    when :plus then
      retval += "_add_#{type}("
      retval += lop + "," + rop + ")"
    when :minus then
      retval += "_sub_#{type}("
      retval += lop + "," + rop + ")"
    when :mult then
      retval += "_mul_#{type}("
      retval += lop + "," + rop + ")"
    when :div then
      abort "division of integer value is not supported in AVX2" if @type =~ /(S|U)(32|64)/
      retval += "_div_#{type}("
      retval += lop + "," + rop + ")"
    when :lt then
      retval += "_cmp_#{type}("
      retval += lop + "," + rop + ",_CMP_LT_OS)"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :le then
      retval += "_cmp_#{type}("
      retval += lop + "," + rop + ",_CMP_LE_OS)"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :gt then
      retval += "_cmp_#{type}("
      retval += lop + "," + rop + ",_CMP_GT_OS)"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :ge then
      retval += "_cmp_#{type}("
      retval += lop + "," + rop + ",_CMP_GE_OS)"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :eq then
      retval += "_cmp_#{type}("
      retval += lop + "," + rop + ",_CMP_EQ_OQ)"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :neq then
      retval += "_cmp_#{type}("
      retval += lop + "," + rop + ",_CMP_NEQ_OQ)"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :and then
      retval += "_and_#{type}("
      retval += lop + "," + rop + ")"
    when :or then
      retval += "_or_#{type}("
      retval += lop + "," + rop + ")"
    when :land then
      retval += "_and_#{type}("
      retval += lop + "," + rop + ")"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :landnot then # PIKG assume andnot(a,b) = a & (!b), but AVX2 andnot returns (!a) & b
      retval += "_andnot_#{type}("
      retval += rop + "," + lop + ")"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :lor then
      retval += "_or_#{type}("
      retval += lop + "," + rop + ")"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :not then
      retval += "_xor_#{type}("
      retval += lop + "," + "_mm256_castsi256_#{type}(_mm256_set1_epi32(-1))" + ")"
      retval = "_mm256_cast#{type}_ps(" + retval + ")" if type != "ps"
    when :dot then
      if @rop == "x" || @rop == "y" || @rop == "z" || @rop == "w"
        retval=lop+"."
        retval += ["v0","v1","v2","v3"][["x","y","z","w"].index(@rop)]
        #rop
      elsif @rop == "v0" || @rop == "v1" || @rop == "v2" || @rop == "v3"
        retval=lop+"." + @rop
        #rop
      else
        #retval = "#{lop}.#{rop}"
        retval = "_mm256_set1_#{type}#{set1_suffix}(#{lop}.#{rop})"
      end
    when :func then
      retval = @lop
      #retval += "_mask"
      retval += "(" + rop + ")"
    when :array then
      retval = lop+"[" + rop + "]"
    else
      abort "error: unsupported operator #{@operator} for AVX2"
    end
    retval
  end
end

class TableDecl
  def convert_to_code_avx2(conversion_code)
    ret = ""
    simd_width = get_simd_width_avx2(@type)
    nreg = (@table.vals.length/simd_width)
    abort "size of table must be multiple of SIMD width (#{simd_width})" if @table.vals.length%simd_width != 0 || nreg < 1

    for i in 0...nreg
      tmp = self.dup
      tmp.name += "_tmp#{i}"
      ret += tmp.convert_to_code("reference")
      ret += "#{get_declare_type_avx2(@type)} #{@name};\n"
      ret += "#{name} = _mm256_load_ps(#{tmp.name});\n"
    end
    ret
  end
end

class FloatingPoint
  def convert_to_code_avx2(conversion_type)
    case @type
    when "F64"
      "_mm256_set1_pd(#{@val})"
    when "F32"
      "_mm256_set1_ps(#{@val})"
    else
      abort "error: unsupported floating point type #{@type}"
    end
  end
end

class IntegerValue
  def convert_to_code_avx2(conversion_type)
    case @type
    when "S64"
      "_mm256_set1_epi64x(#{@val})"
    when "S32"
      "_mm256_set1_epi32(#{@val})"
    when "S16"
      "_mm256_set1_epi16(#{@val})"
    when "U64"
      "_mm256_set1_epu64(#{@val})"
    when "U32"
      "_mm256_set1_epu32(#{@val})"
    when "U16"
      "_mm256_set1_epu16(#{@val})"
    end
  end
end

class String
  def convert_to_code_avx2(conversion_type,h=$varhash)
    name = get_name(self)
    #abort "error: undefined reference to #{name} in convert_to_code_a64fx of String"
    s = self
    if h[name] != nil
      iotype = h[name][0]
      if iotype == "MEMBER"
        type = self.get_type
        case type
        when "F64"
          s = "_mm256_set1_pd("+self+")"
        when "F32"
          s = "_mm256_set1_ps("+self+")"
        when "S64"
          s = "_mm256_set1_eps64("+self+")"
        when "S32"
          s = "_mm256_set1_eps32("+self+")"
        when "U64"
          s = "_mm256_set1_epu64("+self+")"
        when "U32"
          s = "_mm256_set1_epu32("+self+")"
        end
      end
    end
    s
  end
end

