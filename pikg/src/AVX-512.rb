require_relative "common.rb"

def get_declare_type_avx512(type)
  case type
  when "F64" then
    decl = "__m512d"
  when "F32" then
    decl = "__m512"
  when "F64vec" then
    decl = "__m512dx3"
  when "F32vec" then
    decl = "__m512x3"
  when "F64vec3" then
    decl = "__m512dx3"
  when "F32vec3" then
    decl = "__m512x3"
  when "F64vec4" then
    decl = "__m512dx4"
  when "F32vec4" then
    decl = "__m512x4"
  when "S64" then
    decl = "__m512i"
  when "S32" then
    decl = "__m512i"
  when "S16" then
    decl = "__m512i"
  when "U64" then
    decl = "__m512i"
  when "U32" then
    decl = "__m512i"
  when "U16" then
    decl = "__m512i"
  when "B64" then
    decl = "__mmask8"
  when "B32" then
    decl = "__mmask16"
  when "B16" then
    decl = "__mmask32"
  else
    abort "error: unsupported declaration of type #{type} for AVX-512"
  end
  decl
end

def get_type_suffix_avx512(type)
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
    suffix = "mask8"
  when "B32" then
    suffix = "mask16"
  when "B16" then
    suffix = "mask32"
  else
    abort "error: unsupported suffix of type #{type} for AVX-512"
  end
  suffix
end

class Kernelprogram
  def reserved_func_def_avx512(conversion_type)
    code = ""
    code += "__m512 rsqrt(__m512 op){\n"
    #code += "__m512 rinv = _mm512_rsqrt14_ps(op);\n"
    #code += "__m512 h = _mm512_mul_ps(op,rinv);\n"
    #code += "h = _mm512_fnmadd_ps(h,rinv,_mm512_set1_ps(1.f));\n"
    #code += "__m512 poly = _mm512_fmadd_ps(h,_mm512_set1_ps(0.375f),_mm512_set1_ps(0.5f));\n"
    #code += "poly = _mm512_mul_ps(poly,h);\n"
    #code += "rinv = _mm512_fmadd_ps(rinv,poly,rinv);\n"
    #code += "return rinv;\n"

    code += "return _mm512_rsqrt14_ps(op);\n"
    #code += "return _mm512_rsqrt28_ps(op);\n" # only for AVX-512ER
    code += "}\n"

    code += "__m512 sqrt(__m512 op){ return _mm512_sqrt_ps(op); }\n"
    # http://math-koshimizu.hatenablog.jp/entry/2017/07/28/083000
    code += "__m512 inv(__m512 op){\n"
    code += "__m512 x1 = _mm512_rcp14_ps(op);\n"
    code += "__m512 x2 = _mm512_fnmadd_ps(op,x1,_mm512_set1_ps(2.f));\n"
    code += "x2 = _mm512_mul_ps(x2,x1);\n"
    code += "__m512 ret = _mm512_fnmadd_ps(op,x2,_mm512_set1_ps(2.f));\n"
    code += "ret = _mm512_mul_ps(ret,x2);\n"
    code += "return ret;\n"
    code += "}\n"

    code += "__m512d rsqrt(__m512d op){\n"
    code += "__m512d rinv = _mm512_rsqrt14_pd(op);\n"
    code += "__m512d h = _mm512_mul_pd(op,rinv);\n"
    code += "h = _mm512_fnmadd_pd(h,rinv,_mm512_set1_pd(1.0));\n"
    code += "__m512d poly = _mm512_fmadd_pd(h,_mm512_set1_pd(0.375),_mm512_set1_pd(0.5));\n"
    code += "poly = _mm512_mul_pd(poly,h);\n"
    code += "return _mm512_fmadd_pd(rinv,poly,rinv);\n"
    code += "}\n"
    
    code += "__m512d max(__m512d a,__m512d b){ return _mm512_max_pd(a,b);}\n"
    code += "__m512d min(__m512d a,__m512d b){ return _mm512_min_pd(a,b);}\n"
    code += "__m512  max(__m512  a,__m512  b){ return _mm512_max_ps(a,b);}\n"
    code += "__m512  min(__m512  a,__m512  b){ return _mm512_min_ps(a,b);}\n"
    code += "__m512i max(__m512i a,__m512i b){ return _mm512_max_epi32(a,b);}\n"
    code += "__m512i min(__m512i a,__m512i b){ return _mm512_min_epi32(a,b);}\n"

    code += "__m512d table(__m512d tab,__m512i index){ return _mm512_permutexvar_pd(index,tab);}\n"
    code += "__m512  table(__m512  tab,__m512i index){ return _mm512_permutexvar_ps(index,tab);}\n"

    code += "__m512d to_double(__m512i op){ return _mm512_cvt_roundepi64_pd(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    code += "__m512  to_float(__m512i op){ return _mm512_cvt_roundepi32_ps(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    #code += "__m512d to_float(__m512i op){ return _mm512_cvt_roundepu64_pd(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    #code += "__m512  to_float(__m512i op){ return _mm512_cvt_roundepu32_ps(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    code += "__m512i  to_long(__m512d op){ return _mm512_cvt_roundpd_epi64(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    code += "__m512i  to_int(__m512  op){ return _mm512_cvt_roundps_epi32(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    code += "__m512i  to_ulong(__m512d op){ return _mm512_cvt_roundpd_epu64(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
    code += "__m512i  to_uint(__m512  op){ return _mm512_cvt_roundps_epu32(op,(_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));}\n"
  end
end

class FuncCall
  def convert_to_code_avx512(conversion_type)
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
  def convert_to_code_avx512(conversion_type)
    if @operator != :array && @operator != :func
      if [:eq,:neq,:gt,:ge,:lt,:le].index(@operator)
        suffix = get_type_suffix_avx512(@lop.get_type)
      else
        suffix = get_type_suffix_avx512(self.get_type)
      end
    end
    retval = "_mm512"
    #suffix = get_type_suffix_avx512(@type) if @operator != :array
    case @operator
    when :uminus then
      retval += "_sub_#{suffix}("
      retval += "_mm512_set1_#{suffix}((PIKG::#{@type})0.0)," + @lop.convert_to_code(conversion_type) + ")"
    when :plus then
      retval += "_add_#{suffix}("
      retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
    when :minus then
      retval += "_sub_#{suffix}("
      retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
    when :mult then
      retval += "_mul_#{suffix}("
      retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
    when :div then
      retval += "_div_#{suffix}("
      retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
    when :lt then
      retval += "_cmp_#{suffix}_mask("
      retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type)
      if suffix =~ /epi/
        retval += ",_MM_CMPINT_LT)"
      else
        retval += ",_CMP_LT_OQ)"
      end
    when :le then
      retval += "_cmp_#{suffix}_mask("
      retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type)
      if suffix =~ /epi/
        retval += ",_MM_CMPINT_LE)"
      else
        retval += ",_CMP_LE_OQ)"
      end
    when :gt then
      retval += "_cmp_#{suffix}_mask("
      retval += @rop.convert_to_code(conversion_type) + "," + @lop.convert_to_code(conversion_type)
      if suffix =~ /epi/
        retval += ",_MM_CMPINT_LT)"
      else
        retval += ",_CMP_LT_OQ)"
      end
    when :ge then
      retval += "_cmp_#{suffix}_mask("
      retval += @rop.convert_to_code(conversion_type) + "," + @lop.convert_to_code(conversion_type)
      if suffix =~ /epi/
        retval += ",_MM_CMPINT_LE)"
      else
        retval += ",_CMP_LE_OQ)"
      end
    when :eq then
      retval += "_cmp_#{suffix}_mask("
      retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type)
      if suffix =~ /epi/
        retval += ",_MM_CMPINT_EQ)"
      else
        retval += ",_CMP_EQ_OQ)"
      end
    when :neq then
      retval += "_cmp_#{suffix}_mask("
        retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type)
      if suffix =~ /epi/
        retval += ",_MM_CMPINT_NE)"
      else
        retval += ",_CMP_NEQ_OQ)"
      end
    when :and then
      if @lop.get_type == @lop.get_type
        retval += "_and_#{suffix}_mask("
        retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
      else
        p self
        abort "unsupported operands for :and in convert_to_code_avx512"
      end
    when :or then
      if @lop.get_type == @lop.get_type
        retval += "_or_#{suffix}_mask("
        retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
      else
        p self
        abort "unsupported operands for :or in convert_to_code_avx512"
      end
    when :land then
      if @lop.get_type =~ /B/ && @rop.get_type =~ /B/
        retval = "_kand_#{suffix}("
        retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
      else
        p self
        abort "unsupported operands for :land in convert_to_code_avx512"
      end
    when :lor then
      if @lop.get_type =~ /B/ && @rop.get_type =~ /B/
        retval = "_kor_#{suffix}("
        retval += @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
      else
        p self
        abort "unsupported operands for :lor in convert_to_code_avx512"
      end
    when :landnot then
      if @lop.get_type =~ /B/ && @rop.get_type =~ /B/
        retval = "_kand_#{suffix}("
        retval += @lop.convert_to_code(conversion_type) + ",_knot_#{suffix}(" + @rop.convert_to_code(conversion_type) + "))"
      else
        p self
        abort "unsupported operands for :land in convert_to_code_avx512"
      end
    when :not then
      if true #@lop.get_type =~ /B/
        retval = "_knot_#{suffix}("
        retval += @lop.convert_to_code(conversion_type) + ")"
      else
        p self
        abort "unsupported operands for :not in convert_to_code_avx512"
      end
    when :dot then
      if @rop == "x" || @rop == "y" || @rop == "z" || @rop == "w"
        retval=@lop.convert_to_code(conversion_type)+"."
        retval += ["v0","v1","v2","v3"][["x","y","z","w"].index(@rop)]
        #@rop.convert_to_code(conversion_type)
      elsif @rop == "v0" || @rop == "v1" || @rop == "v2" || @rop == "v3"
        retval=@lop.convert_to_code(conversion_type)+"." + @rop
        #@rop.convert_to_code(conversion_type)
      else
        #retval = "#{@lop.convert_to_code(conversion_type)}.#{@rop.convert_to_code(conversion_type)}"
        retval = "_mm512_set1_#{suffix}(#{@lop.convert_to_code(conversion_type)}.#{@rop.convert_to_code(conversion_type)})"
      end
    when :func then
      retval = @lop
      retval += "(" + @rop.convert_to_code(conversion_type) + ")"
    when :array then
      retval = @lop+"[" + @rop.convert_to_code(conversion_type) + "]"
    else
      abort "error: unsupported operator #{@operator} for AVX-512"
    end
    retval
  end
end

class TableDecl
  def convert_to_code_avx512(conversion_type)
    ret = ""
    nelem = get_num_elem(@type,conversion_type)
    nreg = (@table.vals.length/nelem)
    abort "size of table must be multiple of # of elements of SIMD register (#{nelem})" if @table.vals.length%nelem != 0 || nreg < 1

    for i in 0...nreg
      tmp = self.dup
      tmp.name += "_tmp#{i}"
      ret += tmp.convert_to_code("reference")
      ret += "#{get_declare_type_avx512(@type)} #{@name};\n"
      ret += "#{name} = _mm512_load_ps(#{tmp.name});\n"
    end
    ret
  end
end

class FloatingPoint
  def convert_to_code_avx512(conversion_type)
    case @type
    when "F64"
      "_mm512_set1_pd(#{@val})"
    when "F32"
      "_mm512_set1_ps(#{@val})"
    else
      abort "error: unsupported floating point type #{@type}"
    end
  end
end

class IntegerValue
  def convert_to_code_avx512(conversion_type)
    case @type
    when "S64"
      "_mm512_set1_epi64(#{@val})"
    when "S32"
      "_mm512_set1_epi32(#{@val})"
    when "S16"
      "_mm512_set1_epi16(#{@val})"
    when "U64"
      "_mm512_set1_epu64(#{@val})"
    when "U32"
      "_mm512_set1_epu32(#{@val})"
    when "U16"
      "_mm512_set1_epu16(#{@val})"
    end
  end
end

class String
  def convert_to_code_avx512(conversion_type,h=$varhash)
    name = get_name(self)
    #abort "error: undefined reference to #{name} in convert_to_code_avx512 of String"
    s = self
    if h[name] != nil
      iotype = h[name][0]
      if iotype == "MEMBER"
        type = self.get_type
        case type
        when "F64"
          s = "_mm512_set1_pd("+self+")"
        when "F32"
          s = "_mm512_set1_ps("+self+")"
        when "S64"
          s = "_mm512_set1_eps64("+self+")"
        when "S32"
          s = "_mm512_set1_eps32("+self+")"
        when "U64"
          s = "_mm512_set1_epu64("+self+")"
        when "U32"
          s = "_mm512_set1_epu32("+self+")"
        end
      end
    end
    s
  end
end

