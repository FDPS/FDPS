require_relative "common.rb"
require_relative "loop_fission.rb"

def get_declare_type_a64fx(type)
  case type
  when "F64" then
    decl = "svfloat64_t"
  when "F32" then
    decl = "svfloat32_t"
  when "S64" then
    decl = "svint64_t"
  when "S32" then
    decl = "svint32_t"
  when "U64" then
    decl = "svuint64_t"
  when "U32" then
    decl = "svuint32_t"
  when "F64vec1" then 
    decl = "svfloat64_t"
  when "F32vec1" then 
    decl = "svfloat32_t"
  when "F64vec2" then 
    decl = "svfloat64x2_t"
  when "F32vec2" then 
    decl = "svfloat32x2_t"
  when "F64vec3" then 
    decl = "svfloat64x3_t"
  when "F32vec3" then 
    decl = "svfloat32x3_t"
  when "F64vec4" then 
    decl = "svfloat64x4_t"
  when "F32vec4" then 
    decl = "svfloat32x4_t"
  when "F64vec" then 
    decl = "svfloat64x3_t"
  when "F32vec" then 
    decl = "svfloat32x3_t"
  when "B32" then
    decl = "svbool_t"
  when "B64" then
    decl = "svbool_t"
  else
    abort "error: unsupported vector type of #{type} for A64FX"
    decl = "UNSUPPORTED_VECTOR_TYPE"
  end
  decl
end

def get_declare_scalar_type(type,h=$varhash)
  case type
  when "F64" then
    decl = "float64_t"
  when "F32" then
    decl = "float32_t"
  when "S64" then
    decl = "int64_t"
  when "S32" then
    decl = "int32_t"
  when "U64" then
    decl = "uint64_t"
  when "U32" then
    decl = "uint32_t"
  when "F64vec1" then 
    decl = "float64_t"
  when "F32vec1" then 
    decl = "float32_t"
  when "F64vec2" then 
    decl = "float64x2_t"
  when "F32vec2" then 
    decl = "float32x2_t"
  when "F64vec3" then 
    decl = "float64x3_t"
  when "F32vec3" then 
    decl = "float32x3_t"
  when "F64vec4" then 
    decl = "float64x4_t"
  when "F32vec4" then 
    decl = "float32x4_t"
  when "F64vec" then 
    decl = "float64x3_t"
  when "F32vec" then 
    decl = "float32x3_t"
  else
    abort "error: unsupported scalar type of #{@type} for A64FX (get_declare_scalar_type)"
  end
  decl
end

def get_type_suffix_a64fx(type)
  case type
  when "F64" then
    suffix = "f64"
  when "F32" then
    suffix = "f32"
  when "S64" then
    suffix = "s64"
  when "S32" then
    suffix = "s32"
  when "U64" then
    suffix = "u64"
  when "U32" then
    suffix = "u32"
  when /F64vec/ then
    suffix = "f64"
  when /F32vec/ then
    suffix = "f32"
  when "B32" then
    suffix = "b32"
  when "B64" then
    suffix = "b64"
  else
    abort "error: unsupported scalar type of #{@type} for A64FX (get_type_suffix)"
  end
  suffix
end

class Kernelprogram
  def reserved_func_def_a64fx(conversion_type)
    code = ""
    code += "svfloat32_t rsqrt(svbool_t pg,svfloat32_t op){\n"
    code += "svfloat32_t rinv = svrsqrte_f32(op);\n"
    code += "svfloat32_t h = svmul_f32_z(pg,op,rinv);\n"
    code += "h = svmsb_n_f32_z(pg,h,rinv,1.f);\n"
    code += "svfloat32_t poly = svmad_n_f32_z(pg,h,svdup_f32(0.375f),0.5f);\n"
    code += "poly = svmul_f32_z(pg,poly,h);\n"
    code += "rinv = svmad_f32_z(pg,rinv,poly,rinv);\n"
    code += "return rinv;\n"
    code += "}\n"

    code += "svfloat64_t rsqrt(svbool_t pg,svfloat64_t op){\n"
    code += "svfloat64_t rinv = svrsqrte_f64(op);\n"
    code += "svfloat64_t h = svmul_f64_z(pg,op,rinv);\n"
    code += "h = svmsb_n_f64_z(pg,h,rinv,1.f);\n"
    code += "svfloat64_t poly = svmad_n_f64_z(pg,h,svdup_f64(0.375f),0.5f);\n"
    code += "poly = svmul_f64_z(pg,poly,h);\n"
    code += "rinv = svmad_f64_z(pg,rinv,poly,rinv);\n"
    code += "return rinv;\n"
    code += "}\n"

    code += "svfloat32x2_t svdup_n_f32x3(PIKG::F32vec2 v){\n"
    code += "  svfloat32x2_t ret;\n"
    code += "  ret.v0 = svdup_n_f32(v.x);\n  ret.v1 = svdup_n_f32(v.y);\n"
    code += "  return ret;\n"
    code +="}\n"
    code += "svfloat32x3_t svdup_n_f32x3(PIKG::F32vec v){\n"
    code += "  svfloat32x3_t ret;\n"
    code += "  ret.v0 = svdup_n_f32(v.x);\n  ret.v1 = svdup_n_f32(v.y);\n  ret.v2 = svdup_n_f32(v.z);\n"
    code += "  return ret;\n"
    code +="}\n"
    code += "svfloat32x4_t svdup_n_f32x4(PIKG::F32vec4 v){\n"
    code += "  svfloat32x4_t ret;\n"
    code += "  ret.v0 = svdup_n_f32(v.x);\n  ret.v1 = svdup_n_f32(v.y);\n  ret.v2 = svdup_n_f32(v.z);\n  ret.v3 = svdup_n_f32(v.w);\n"
    code += "  return ret;\n"
    code +="}\n"

    code += "svfloat64x2_t svdup_n_f64x3(PIKG::F64vec2 v){\n"
    code += "  svfloat64x2_t ret;\n"
    code += "  ret.v0 = svdup_n_f64(v.x);\n  ret.v1 = svdup_n_f64(v.y);\n"
    code += "  return ret;\n"
    code +="}\n"
    code += "svfloat64x3_t svdup_n_f64x3(PIKG::F64vec v){\n"
    code += "  svfloat64x3_t ret;\n"
    code += "  ret.v0 = svdup_n_f64(v.x);\n  ret.v1 = svdup_n_f64(v.y);\n  ret.v2 = svdup_n_f64(v.z);\n"
    code += "  return ret;\n"
    code +="}\n"
    code += "svfloat64x4_t svdup_n_f64x4(PIKG::F64vec4 v){\n"
    code += "  svfloat64x4_t ret;\n"
    code += "  ret.v0 = svdup_n_f64(v.x);\n  ret.v1 = svdup_n_f64(v.y);\n  ret.v2 = svdup_n_f64(v.z);\n  ret.v3 = svdup_n_f64(v.w);\n"
    code += "  return ret;\n"
    code +="}\n"

    code += "svfloat32_t sqrt(svbool_t pg,svfloat32_t op){ return svsqrt_f32_z(pg,op); }\n"
    code += "svfloat64_t sqrt(svbool_t pg,svfloat64_t op){ return svsqrt_f64_z(pg,op); }\n"
    # http://math-koshimizu.hatenablog.jp/entry/2017/07/28/083000
    code += "svfloat32_t inv(svbool_t pg,svfloat32_t op){\n"
    code += "svfloat32_t x1 = svrecpe_f32(op);\n"
    code += "svfloat32_t x2 = svmsb_n_f32_z(pg,op,x1,2.f);\n"
    code += "x2 = svmul_f32_z(pg,x2,x1);\n"
    code += "svfloat32_t ret = svmsb_n_f32_z(pg,op,x2,2.f);\n"
    code += "ret = svmul_f32_z(pg,ret,x2);\n"
    code += "return ret;\n"
    code += "}\n"

    code += "svfloat64_t inv(svbool_t pg,svfloat64_t op){\n"
    code += "svfloat64_t x1 = svrecpe_f64(op);\n"
    code += "svfloat64_t x2 = svmsb_n_f64_z(pg,op,x1,2.f);\n"
    code += "x2 = svmul_f64_z(pg,x2,x1);\n"
    code += "svfloat64_t ret = svmsb_n_f64_z(pg,op,x2,2.f);\n"
    code += "ret = svmul_f64_z(pg,ret,x2);\n"
    code += "return ret;\n"
    code += "}\n"

    code += "svfloat64_t max(svbool_t pg,svfloat64_t a,svfloat64_t b){ return svmax_f64_z(pg,a,b);}\n"
    code += "svfloat64_t min(svbool_t pg,svfloat64_t a,svfloat64_t b){ return svmin_f64_z(pg,a,b);}\n"
    code += "svuint64_t max(svbool_t pg,svuint64_t a,svuint64_t b){ return svmax_u64_z(pg,a,b);}\n"
    code += "svuint64_t min(svbool_t pg,svuint64_t a,svuint64_t b){ return svmin_u64_z(pg,a,b);}\n"
    code += "svint64_t max(svbool_t pg,svint64_t a,svint64_t b){ return svmax_s64_z(pg,a,b);}\n"
    code += "svint64_t min(svbool_t pg,svint64_t a,svint64_t b){ return svmin_s64_z(pg,a,b);}\n"
    code += "svfloat32_t max(svbool_t pg,svfloat32_t a,svfloat32_t b){ return svmax_f32_z(pg,a,b);}\n"
    code += "svfloat32_t min(svbool_t pg,svfloat32_t a,svfloat32_t b){ return svmin_f32_z(pg,a,b);}\n"
    code += "svuint32_t max(svbool_t pg,svuint32_t a,svuint32_t b){ return svmax_u32_z(pg,a,b);}\n"
    code += "svuint32_t min(svbool_t pg,svuint32_t a,svuint32_t b){ return svmin_u32_z(pg,a,b);}\n"
    code += "svint32_t max(svbool_t pg,svint32_t a,svint32_t b){ return svmax_s32_z(pg,a,b);}\n"
    code += "svint32_t min(svbool_t pg,svint32_t a,svint32_t b){ return svmin_s32_z(pg,a,b);}\n"

    code += "svfloat64_t table(svfloat64_t tab,svuint64_t index){ return svtbl_f64(tab,index);}\n"
    code += "svfloat32_t table(svfloat32_t tab,svuint32_t index){ return svtbl_f32(tab,index);}\n"

    code += "svfloat64_t to_float(svbool_t pg, svint64_t op){ return svcvt_f64_s64_z(pg,op);}\n"
    code += "svfloat32_t to_float(svbool_t pg, svint32_t op){ return svcvt_f32_s32_z(pg,op);}\n"
    code += "svfloat64_t to_float(svbool_t pg,svuint64_t op){ return svcvt_f64_u64_z(pg,op);}\n"
    code += "svfloat32_t to_float(svbool_t pg,svuint32_t op){ return svcvt_f32_u32_z(pg,op);}\n"
    code += "svint64_t  to_int(svbool_t pg,svfloat64_t op){ return svcvt_s64_f64_z(pg,op);}\n"
    code += "svint32_t  to_int(svbool_t pg,svfloat32_t op){ return svcvt_s32_f32_z(pg,op);}\n"
    code += "svuint64_t to_uint(svbool_t pg,svfloat64_t op){ return svcvt_u64_f64_z(pg,op);}\n"
    code += "svuint32_t to_uint(svbool_t pg,svfloat32_t op){ return svcvt_u32_f32_z(pg,op);}\n"
    # code += "void transpose4x4(svfloat32x4_t& v){\n"
    # code += "const unsigned int tmp[16] = { 0, 2, 1, 3, 4, 6, 5, 7, 8,10, 9,11,12,14,13,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v.v0 = svtbl_f32(v.v0,index);\n" # { 0, 2, 1, 3, 4, 6, 5, 7, 8,10, 9,11,12,14,13,15}
    # code += "v.v1 = svtbl_f32(v.v1,index);\n" #
    # code += "v.v2 = svtbl_f32(v.v2,index);\n" #
    # code += "v.v3 = svtbl_f32(v.v3,index);\n" #
    # code += "svfloat64_t xy0 = svreinterpret_f64_f32(svtrn1_f32(v.v0,v.v1));\n" # { 0,16, 1,17, 4,20, 5,21, 8,24, 9,25,12,28,13,29}
    # code += "svfloat64_t xy1 = svreinterpret_f64_f32(svtrn2_f32(v.v0,v.v1));\n" # { 2,18, 3,19, 6,22, 7,23,10,26,11,27,14,30,15,31}
    # code += "svfloat64_t zw0 = svreinterpret_f64_f32(svtrn1_f32(v.v2,v.v3));\n" # {32,48,33,49,36,52,37,53,40,56,41,57,44,60,45,61}
    # code += "svfloat64_t zw1 = svreinterpret_f64_f32(svtrn2_f32(v.v2,v.v3));\n" # {34,50,35,51,38,54,39,55,42,58,43,59,46,62,47,63}
    # code += "v.v0 = svreinterpret_f32_f64(svtrn1_f64(xy0,zw0));\n" # { 0,16,32,48, 4,20,36,52, 8,24,40,56,12,28,44,60}
    # code += "v.v1 = svreinterpret_f32_f64(svtrn2_f64(xy0,zw0));\n" # { 1,17,33,49, 5,21,37,53, 9,25,41,57,13,29,45,61}
    # code += "v.v2 = svreinterpret_f32_f64(svtrn1_f64(xy1,zw1));\n" # { 2,18,34,50, 6,22,38,54,10,26,42,58,14,30,46,62}
    # code += "v.v3 = svreinterpret_f32_f64(svtrn2_f64(xy1,zw1));\n" # { 3,19,35,51, 7,23,39,55,11,27,43,59,15,31,47,63}
    # code += "}\n"

    # code += "void gather8(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "const unsigned int tmp[16] = {  0,  8,  1,  9,  2, 10,  3, 11, 4, 12,  5, 13,  6, 14,  7, 15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "svfloat64_t a = svreinterpret_f64_f32(svtbl_f32(v0,index));\n"
    # code += "svfloat64_t b = svreinterpret_f64_f32(svtbl_f32(v1,index));\n"
    # code += "svfloat64_t c = svreinterpret_f64_f32(svtbl_f32(v2,index));\n"
    # code += "svfloat64_t d = svreinterpret_f64_f32(svtbl_f32(v3,index));\n"
    # code += "svfloat64_t e = svreinterpret_f64_f32(svtbl_f32(v4,index));\n"
    # code += "svfloat64_t f = svreinterpret_f64_f32(svtbl_f32(v5,index));\n"
    # code += "svfloat64_t g = svreinterpret_f64_f32(svtbl_f32(v6,index));\n"
    # code += "svfloat64_t h = svreinterpret_f64_f32(svtbl_f32(v7,index));\n"

    # code += "svfloat64_t ae0 = svzip1_f64(a,e);\n" # {  0,  8, 64, 72,  1,  9, 65, 73,  2, 10, 66, 74,  3, 11, 67, 75} is composed of v0,v4
    # code += "svfloat64_t ae1 = svzip2_f64(a,e);\n" # {  4, 12, 68, 76,  5, 13, 69, 77,  6, 14, 70, 78,  7, 15, 71, 79} is composed of v0,v4
    # code += "svfloat64_t bf0 = svzip1_f64(b,f);\n"
    # code += "svfloat64_t bf1 = svzip2_f64(b,f);\n" 
    # code += "svfloat64_t cg0 = svzip1_f64(c,g);\n" # { 32, 40, 96,104, 33, 41, 97,105, 34, 42, 98,106, 35, 43, 99,107} is composed of v2,v6
    # code += "svfloat64_t cg1 = svzip2_f64(c,g);\n" # { 36, 44,100,108, 37, 45,101,109, 38, 46,102,110, 39, 47,103,111} is composed of v2,v6
    # code += "svfloat64_t dh0 = svzip1_f64(d,h);\n"
    # code += "svfloat64_t dh1 = svzip2_f64(d,h);\n" 
    
    # code += "svfloat64_t aceg0 = svzip1_f64(ae0,cg0);\n" # {  0,  8, 32, 40, 64, 72, 96,104,  1,  9, 33, 41, 65, 73, 97,105} is composed of v0,v2,v4,v6
    # code += "svfloat64_t aceg1 = svzip2_f64(ae0,cg0);\n" # {  2, 10, 34, 42, 66, 74, 98,106,  3, 11, 35, 43, 67, 75, 99,107} is composed of v0,v2,v4,v6
    # code += "svfloat64_t aceg2 = svzip1_f64(ae1,cg1);\n"
    # code += "svfloat64_t aceg3 = svzip2_f64(ae1,cg1);\n"
    # code += "svfloat64_t bdfh0 = svzip1_f64(bf0,dh0);\n" # { 16, 24, 48, 56, 80, 88,112,120, 17, 25, 49, 57, 81, 89,113,121} is composed of v1,v3,v5,v7
    # code += "svfloat64_t bdfh1 = svzip2_f64(bf0,dh0);\n" # { 18, 26, 50, 58, 82, 90,114,122, 19, 27, 51, 59, 83, 91,115,123} is composed of v1,v3,v5,v7
    # code += "svfloat64_t bdfh2 = svzip1_f64(bf1,dh1);\n"
    # code += "svfloat64_t bdfh3 = svzip2_f64(bf1,dh1);\n"
    
    # code += "v0 = svreinterpret_f32_f64(svzip1_f64(aceg0,bdfh0));\n" # {  0,  8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96,104,112,120}
    # code += "v1 = svreinterpret_f32_f64(svzip2_f64(aceg0,bdfh0));\n" # {  1,  9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97,105,113,121}
    # code += "v2 = svreinterpret_f32_f64(svzip1_f64(aceg1,bdfh1));\n" # {  2, 10, 18, 26, 34, 42, 50, 58, 66, 74, 82, 90, 98,106,113,122}
    # code += "v3 = svreinterpret_f32_f64(svzip2_f64(aceg1,bdfh1));\n" # {  3, 11, 19, 27}
    # code += "v4 = svreinterpret_f32_f64(svzip1_f64(aceg2,bdfh2));\n" # {  4, 12, 20, 28}
    # code += "v5 = svreinterpret_f32_f64(svzip2_f64(aceg2,bdfh2));\n" # {  5, 13, 21, 29}
    # code += "v6 = svreinterpret_f32_f64(svzip1_f64(aceg3,bdfh3));\n" # {  6, 14, 22, 30}
    # code += "v7 = svreinterpret_f32_f64(svzip2_f64(aceg3,bdfh3));\n" # {  7, 15, 23, 31}
    # code += "}\n"

    # code += "void gather5(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,10,11,12,5,6,7,8,9,13,14,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code ++ "v0 = svtbl_f32(v0,index);\n"
    # code ++ "v1 = svtbl_f32(v1,index);\n"
    # code ++ "v2 = svtbl_f32(v2,index);\n"
    # code ++ "v3 = svtbl_f32(v3,index);\n"
    # code ++ "v4 = svtbl_f32(v4,index);\n"
    # code ++ "v5 = svtbl_f32(v5,index);\n"
    # code ++ "v6 = svtbl_f32(v6,index);\n"
    # code ++ "v7 = svtbl_f32(v7,index);\n"
    # code += "gather8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "}\n"     
    # code += "void gather6(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,5,12,13,6,7,8,9,10,11,14,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code ++ "v0 = svtbl_f32(v0,index);\n"
    # code ++ "v1 = svtbl_f32(v1,index);\n"
    # code ++ "v2 = svtbl_f32(v2,index);\n"
    # code ++ "v3 = svtbl_f32(v3,index);\n"
    # code ++ "v4 = svtbl_f32(v4,index);\n"
    # code ++ "v5 = svtbl_f32(v5,index);\n"
    # code ++ "v6 = svtbl_f32(v6,index);\n"
    # code ++ "v7 = svtbl_f32(v7,index);\n"
    # code += "gather8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "}\n"     
    # code += "void gather7(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,5,6,14,7,8,9,10,11,12,13,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code ++ "v0 = svtbl_f32(v0,index);\n"
    # code ++ "v1 = svtbl_f32(v1,index);\n"
    # code ++ "v2 = svtbl_f32(v2,index);\n"
    # code ++ "v3 = svtbl_f32(v3,index);\n"
    # code ++ "v4 = svtbl_f32(v4,index);\n"
    # code ++ "v5 = svtbl_f32(v5,index);\n"
    # code ++ "v6 = svtbl_f32(v6,index);\n"
    # code ++ "v7 = svtbl_f32(v7,index);\n"
    # code += "gather8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "}\n"     
    
    # code += "void scatter8(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "svfloat32_t ae0 = svzip1_f32(v0,v4);\n" # {  0, 64,  1, 65,  2, 66,  3, 67,  4, 68,  5, 69,  6, 70,  7, 71} is composed of v0,v4
    # code += "svfloat32_t ae1 = svzip2_f32(v0,v4);\n" # {  8, 72,  9, 73, 10, 74, 11, 75, 12, 76, 13, 77, 14, 78, 15, 79} is composed of v0,v4
    # code += "svfloat32_t bf0 = svzip1_f32(v1,v5);\n"
    # code += "svfloat32_t bf1 = svzip2_f32(v1,v5);\n"
    # code += "svfloat32_t cg0 = svzip1_f32(v2,v6);\n"
    # code += "svfloat32_t cg1 = svzip2_f32(v2,v6);\n"
    # code += "svfloat32_t dh0 = svzip1_f32(v3,v7);\n"
    # code += "svfloat32_t dh1 = svzip2_f32(v3,v7);\n"
    # code += "svfloat32_t aceg0 = svzip1_f32(ae0,cg0);\n" # {  0, 32, 64, 96,  1, 33, 65, 97,  2, 34, 66, 98,  3, 35, 67, 99} is composed of v0,v2,v4,v6
    # code += "svfloat32_t aceg1 = svzip2_f32(ae0,cg0);\n" # {  4, 36, 68,100,  5, 37, 69,101,  6, 38, 70,102,  7, 39, 71,103} is composed of v0,v2,v4,v6
    # code += "svfloat32_t aceg2 = svzip1_f32(ae1,cg1);\n"
    # code += "svfloat32_t aceg3 = svzip2_f32(ae1,cg1);\n"
    # code += "svfloat32_t bdfh0 = svzip1_f32(bf0,dh0);\n" # { 16, 48, 80,112, 17, 49, 81,113, 18, 50, 82,114, 19, 51, 83,115} is composed of v1,v3,v5,v7
    # code += "svfloat32_t bdfh1 = svzip2_f32(bf0,dh0);\n" # { 20, 52, 84,116, 21, 53, 85,117, 22, 54, 86,118, 23, 55, 87,119} is composed of v1,v3,v5,v7
    # code += "svfloat32_t bdfh2 = svzip1_f32(bf1,dh1);\n"
    # code += "svfloat32_t bdfh3 = svzip2_f32(bf1,dh1);\n"
    # code += "v0 = svzip1_f32(aceg0,bdfh0);\n" # {  0, 16, 32, 48, 64, 80, 96,112,  1, 17, 33, 49, 65, 81, 97,113}
    # code += "v1 = svzip2_f32(aceg0,bdfh0);\n" # {  2, 18, 34, 50, 66, 82, 98,114,  3, 19, 35, 51, 67, 83, 99,115} 
    # code += "v2 = svzip1_f32(aceg1,bdfh1);\n" # {  4, 20, 36, 52, 68, 84,100,116,  5, 21, 37, 53, 69, 85,101,117}
    # code += "v3 = svzip2_f32(aceg1,bdfh1);\n" # {  6, 22, 38, 54, 70, 86,102,118,  7, 23, 39, 55, 71, 87,103,119}
    # code += "v4 = svzip1_f32(aceg2,bdfh2);\n" # {  8, 24, 40, 56, 72, 88,104,120,  9, 25, 41, 57, 73, 89,105,121}
    # code += "v5 = svzip2_f32(aceg2,bdfh2);\n" # { 10, 26, 42, 58, 74, 90,106,122, 11, 27, 43, 59, 75, 91,107,123}
    # code += "v6 = svzip1_f32(aceg3,bdfh3);\n" # { 12, 28, 44, 60, 76, 92,108,124, 13, 29, 45, 61, 77, 93,109,125}
    # code += "v7 = svzip2_f32(aceg3,bdfh3);\n" # { 14, 30, 46, 62, 78, 94,110,126, 15, 31, 47, 63, 79, 94,111,127
    # code += "}\n"
    
    # code += "void scatter5(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "scatter8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,8,9,10,11,12,5,6,7,13,14,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v0 = svtbl_f32(v0,index);\n"
    # code += "v1 = svtbl_f32(v1,index);\n"
    # code += "v2 = svtbl_f32(v2,index);\n"
    # code += "v3 = svtbl_f32(v3,index);\n"
    # code += "v4 = svtbl_f32(v4,index);\n"
    # code += "v5 = svtbl_f32(v5,index);\n"
    # code += "v6 = svtbl_f32(v6,index);\n"
    # code += "v7 = svtbl_f32(v7,index);\n"
    # code += "}\n"

    # code += "void scatter6(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "scatter8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,5,8,9,10,11,12,13,6,7,14,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v0 = svtbl_f32(v0,index);\n"
    # code += "v1 = svtbl_f32(v1,index);\n"
    # code += "v2 = svtbl_f32(v2,index);\n"
    # code += "v3 = svtbl_f32(v3,index);\n"
    # code += "v4 = svtbl_f32(v4,index);\n"
    # code += "v5 = svtbl_f32(v5,index);\n"
    # code += "v6 = svtbl_f32(v6,index);\n"
    # code += "v7 = svtbl_f32(v7,index);\n"
    # code += "}\n"

    # code += "void scatter7(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){\n"
    # code += "scatter8(v0,v1,v2,v3,v4,v5,v6,v7);\n"
    # code += "const unsigned int tmp[16] = {0,1,2,3,4,5,6,8,9,10,11,12,13,14,7,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v0 = svtbl_f32(v0,index);\n"
    # code += "v1 = svtbl_f32(v1,index);\n"
    # code += "v2 = svtbl_f32(v2,index);\n"
    # code += "v3 = svtbl_f32(v3,index);\n"
    # code += "v4 = svtbl_f32(v4,index);\n"
    # code += "v5 = svtbl_f32(v5,index);\n"
    # code += "v6 = svtbl_f32(v6,index);\n"
    # code += "}\n"

    # code += "void transpose8x8(svfloat32x4_t& v0,svfloat32x4_t v1){\n"
    # code += "const unsigned int tmp[16] = {0, 8, 4,12, 1, 9, 5,13, 2,10, 6,14, 3,11, 7,15};\n"
    # code += "const svuint32_t index = svld1_u32(svptrue_b32(),tmp);\n"
    # code += "v0.v0 = svtbl_f32(v0.v0,index);\n" # {0, 8, 4,12, 1, 9, 5,13, 2,10, 6,14, 3,11, 7,15}
    # code += "v0.v1 = svtbl_f32(v0.v1,index);\n"
    # code += "v0.v2 = svtbl_f32(v0.v2,index);\n"
    # code += "v0.v3 = svtbl_f32(v0.v3,index);\n"
    # code += "v1.v0 = svtbl_f32(v1.v0,index);\n" # {0, 8, 4,12, 1, 9, 5,13, 2,10, 6,14, 3,11, 7,15}
    # code += "v1.v1 = svtbl_f32(v1.v1,index);\n"
    # code += "v1.v2 = svtbl_f32(v1.v2,index);\n"
    # code += "v1.v3 = svtbl_f32(v1.v3,index);\n"
    # code += "svfloat64_t ab0 = svreinterpret_f64_f32(svtrn1_f32(v0.v0,v0.v1));\n" # ab{0, 4, 1, 5, 2, 6, 3, 7}
    # code += "svfloat64_t ab1 = svreinterpret_f64_f32(svtrn2_f32(v0.v0,v0.v1));\n" # ab{8,12, 9,13,10,14,11,15}
    # code += "svfloat64_t cd0 = svreinterpret_f64_f32(svtrn1_f32(v0.v2,v0.v3));\n" # cd{0, 4, 1, 5, 2, 6, 3, 7}
    # code += "svfloat64_t cd1 = svreinterpret_f64_f32(svtrn2_f32(v0.v2,v0.v3));\n" # cd{8,12, 9,13,10,14,11,15}
    # code += "svfloat64_t ef0 = svreinterpret_f64_f32(svtrn1_f32(v1.v0,v1.v1));\n" # ef{0, 4, 1, 5, 2, 6, 3, 7}
    # code += "svfloat64_t ef1 = svreinterpret_f64_f32(svtrn2_f32(v1.v0,v1.v1));\n" # ef{8,12, 9,13,10,14,11,15}
    # code += "svfloat64_t gh0 = svreinterpret_f64_f32(svtrn1_f32(v1.v2,v1.v3));\n" # gh{0, 4, 1, 5, 2, 6, 3, 7}
    # code += "svfloat64_t gh1 = svreinterpret_f64_f32(svtrn2_f32(v1.v2,v1.v3));\n" # gh{8,12, 9,13,10,14,11,15}
    # code += "v0.v0 = svreinterpret_f32_f64(svtrn1_f64(ab0,cd0));\n" # abcd{ 0, 1, 2, 3}
    # code += "v0.v1 = svreinterpret_f32_f64(svtrn2_f64(ab0,cd0));\n" # abcd{ 4, 5, 6, 7}
    # code += "v0.v2 = svreinterpret_f32_f64(svtrn1_f64(ab1,cd1));\n" # abcd{ 8, 9,10,11}
    # code += "v0.v3 = svreinterpret_f32_f64(svtrn2_f64(ab1,cd1));\n" # abcd{12,13,14,15}
    # code += "v1.v0 = svreinterpret_f32_f64(svtrn1_f64(ef0,gh0));\n" # efgh{ 0, 1, 2, 3}
    # code += "v1.v1 = svreinterpret_f32_f64(svtrn2_f64(ef0,gh0));\n" # efgh{ 4, 5, 6, 7}
    # code += "v1.v2 = svreinterpret_f32_f64(svtrn1_f64(ef1,gh1));\n" # efgh{ 8, 9,10,11}
    # code += "v1.v3 = svreinterpret_f32_f64(svtrn2_f64(ef1,gh1));\n" # efgh{12,13,14,15}
    # code += "v0.v0 = svreinterpret_f32_f64(svtrn1_f64(ab0,cd0));\n" # abcdefgh{ 0, 1}
    # code += "v0.v1 = svreinterpret_f32_f64(svtrn2_f64(ab0,cd0));\n" # abcdefgh{ 2, 3}
    # code += "v0.v2 = svreinterpret_f32_f64(svtrn1_f64(ab1,cd1));\n" # abcdefgh{ 4, 5}
    # code += "v0.v3 = svreinterpret_f32_f64(svtrn2_f64(ab1,cd1));\n" # abcdefgh{ 6, 7}
    # code += "v1.v0 = svreinterpret_f32_f64(svtrn1_f64(ef0,gh0));\n" # abcdefgh{ 8, 9}
    # code += "v1.v1 = svreinterpret_f32_f64(svtrn2_f64(ef0,gh0));\n" # abcdefgh{10,11}
    # code += "v1.v2 = svreinterpret_f32_f64(svtrn1_f64(ef1,gh1));\n" # abcdefgh{12,13}
    # code += "v1.v3 = svreinterpret_f32_f64(svtrn2_f64(ef1,gh1));\n" # abcdefgh{14,15}
    # code += ""
    # code += "}\n"

    # code += "void transpose16x16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){\n"
    # code += "svfloat32_t ai0 = svzip1_f32(v0,v8);\n" # {  0, 64,  1, 65,  2, 66,  3, 67,  4, 68,  5, 69,  6, 70,  7, 71} is composed of v0,v8
    # code += "svfloat32_t ai1 = svzip2_f32(v0,v8);\n" # {  8, 72,  9, 73, 10, 74, 11, 75, 12, 76, 13, 77, 14, 78, 15, 79} is composed of v0,v8
    # code += "svfloat32_t bj0 = svzip1_f32(v1,v9);\n"
    # code += "svfloat32_t bj1 = svzip2_f32(v1,v9);\n"
    # code += "svfloat32_t ck0 = svzip1_f32(v2,v10);\n"
    # code += "svfloat32_t ck1 = svzip2_f32(v2,v10);\n"
    # code += "svfloat32_t dl0 = svzip1_f32(v3,v11);\n"
    # code += "svfloat32_t dl1 = svzip2_f32(v3,v11);\n"
    # code += "svfloat32_t em0 = svzip1_f32(v4,v12);\n"
    # code += "svfloat32_t em1 = svzip2_f32(v4,v12);\n"
    # code += "svfloat32_t fn0 = svzip1_f32(v5,v13);\n"
    # code += "svfloat32_t fn1 = svzip2_f32(v5,v13);\n"
    # code += "svfloat32_t go0 = svzip1_f32(v6,v14);\n"
    # code += "svfloat32_t go1 = svzip2_f32(v6,v14);\n"
    # code += "svfloat32_t hp0 = svzip1_f32(v7,v15);\n"
    # code += "svfloat32_t hp1 = svzip2_f32(v7,v15);\n"
    
    # code += "svfloat32_t aeim0 = svzip1_f32(ai0,em0);\n" # {  0, 64,128,192,  1, 65,129,193,  2, 68,130,194,  3, 67,131,195} is composed of v0,v4,v8,v12
    # code += "svfloat32_t aeim1 = svzip2_f32(ai0,em0);\n"
    # code += "svfloat32_t aeim2 = svzip1_f32(ai1,em1);\n"
    # code += "svfloat32_t aeim3 = svzip2_f32(ai1,em1);\n"
    # code += "svfloat32_t bfjn0 = svzip1_f32(bj0,fn0);\n"
    # code += "svfloat32_t bfjn1 = svzip2_f32(bj0,fn0);\n"
    # code += "svfloat32_t bfjn2 = svzip1_f32(bj1,fn1);\n"
    # code += "svfloat32_t bfjn3 = svzip2_f32(bj1,fn1);\n"
    # code += "svfloat32_t cgko0 = svzip1_f32(ck0,go0);\n"
    # code += "svfloat32_t cgko1 = svzip2_f32(ck0,go0);\n"
    # code += "svfloat32_t cgko2 = svzip1_f32(ck1,go1);\n"
    # code += "svfloat32_t cgko3 = svzip2_f32(ck1,go1);\n"
    # code += "svfloat32_t dhlp0 = svzip1_f32(dl0,hp0);\n"
    # code += "svfloat32_t dhlp1 = svzip2_f32(dl0,hp0);\n"
    # code += "svfloat32_t dhlp2 = svzip1_f32(dl1,hp1);\n"
    # code += "svfloat32_t dhlp3 = svzip2_f32(dl1,hp1);\n"
    
    # code += "svfloat32_t acegikmo0 = svzip1_f32(aeim0,cgko0);\n" # {  0, 32, 64, 96,128,160,192,224,  1, 33, 65, 97,129,161,193,225} is composed of v0,v2,v4,v6,v8,v10,v12,v14
    # code += "svfloat32_t acegikmo1 = svzip2_f32(aeim0,cgko0);\n" # {  2, 34, 68, 98,130,162,194,226,  3, 35, 67, 99,131,163,195,227} is composed of v0,v2,v4,v6,v8,v10,v12,v14
    # code += "svfloat32_t acegikmo2 = svzip1_f32(aeim1,cgko1);\n"
    # code += "svfloat32_t acegikmo3 = svzip2_f32(aeim1,cgko1);\n"
    # code += "svfloat32_t acegikmo4 = svzip1_f32(aeim2,cgko2);\n"
    # code += "svfloat32_t acegikmo5 = svzip2_f32(aeim2,cgko2);\n"
    # code += "svfloat32_t acegikmo6 = svzip1_f32(aeim3,cgko3);\n"
    # code += "svfloat32_t acegikmo7 = svzip2_f32(aeim3,cgko3);\n"
    # code += "svfloat32_t bdfhjlnp0 = svzip1_f32(bfjn0,dhlp0);\n" # { 16, 48, 80,112,144,176,208,240, 17, 49, 81,113,145,177,209,241} is composed of v1,v3,v5,v7,v9,v11,v13,v15
    # code += "svfloat32_t bdfhjlnp1 = svzip2_f32(bfjn0,dhlp0);\n" # { 18, 50, 82,114,146,178,210,242, 19, 51, 83,115,147,179,211,243} is composed of v1,v3,v5,v7,v9,v11,v13,v15
    # code += "svfloat32_t bdfhjlnp2 = svzip1_f32(bfjn1,dhlp1);\n"
    # code += "svfloat32_t bdfhjlnp3 = svzip2_f32(bfjn1,dhlp1);\n"
    # code += "svfloat32_t bdfhjlnp4 = svzip1_f32(bfjn2,dhlp2);\n"
    # code += "svfloat32_t bdfhjlnp5 = svzip2_f32(bfjn2,dhlp2);\n"
    # code += "svfloat32_t bdfhjlnp6 = svzip1_f32(bfjn3,dhlp3);\n"
    # code += "svfloat32_t bdfhjlnp7 = svzip2_f32(bfjn3,dhlp3);\n"
    # code += "v0  = svzip1_f32(acegikmo0,bdfhjlnp0);\n" # {  0, 16, 32, 48, 64, 80, 96,112,128,144,160,176,192,208,224,240}
    # code += "v1  = svzip2_f32(acegikmo0,bdfhjlnp0);\n" # {  1, 17, 33, 49, 65, 81, 97,113,129,145,161,177,193,209,225,241}
    # code += "v2  = svzip1_f32(acegikmo1,bdfhjlnp1);\n" # {  2, 18, 34, 50, 66, 82, 98,114,130,146,162,178,194,210,226,242}
    # code += "v3  = svzip2_f32(acegikmo1,bdfhjlnp1);\n" # {  3, 19, 35, 51, 67, 83, 99,115,131,147,163,179,195,211,227,243}
    # code += "v4  = svzip1_f32(acegikmo2,bdfhjlnp2);\n" # {  4, 20, 36, 52, 68, 84,100,116,132,148,164,180,196,212,228,244}
    # code += "v5  = svzip2_f32(acegikmo2,bdfhjlnp2);\n" # {  5, 21, 37, 53, 69, 85,101,117,133,149,165,181,197,213,229,245}
    # code += "v6  = svzip1_f32(acegikmo3,bdfhjlnp3);\n" # {  6, 22, 38, 54, 70, 86,102,118,134,150,166,182,198,214,230,246}
    # code += "v7  = svzip2_f32(acegikmo3,bdfhjlnp3);\n" # {  7, 23, 39, 55, 71, 87,103,119,135,151,167,183,199,215,231,247}
    # code += "v8  = svzip1_f32(acegikmo4,bdfhjlnp4);\n" # {  8, 24, 40, 56, 72, 88,104,120,136,152,168,184,200,216,232,248}
    # code += "v9  = svzip2_f32(acegikmo4,bdfhjlnp4);\n" # {  9, 25, 41, 57, 73, 89,105,121,137,153,169,185,201,217,233,249}
    # code += "v10 = svzip1_f32(acegikmo5,bdfhjlnp5);\n" # { 10, 26, 42, 58, 74, 90,106,122,138,154,170,186,202,218,234,250} 
    # code += "v11 = svzip2_f32(acegikmo5,bdfhjlnp5);\n" # { 11, 27, 43, 59, 75, 91,107,123,139,155,171,187,203,219,235,251}
    # code += "v12 = svzip1_f32(acegikmo6,bdfhjlnp6);\n" # { 12, 28, 44, 60, 76, 92,108,124,140,156,172,188,204,220,236,252}
    # code += "v13 = svzip2_f32(acegikmo6,bdfhjlnp6);\n" # { 13, 29, 45, 61, 77, 93,109,125,141,157,173,189,205,221,237,253}
    # code += "v14 = svzip1_f32(acegikmo7,bdfhjlnp7);\n" # { 14, 30, 46, 62, 78, 94,110,126,142,158,174,190,206,222,238,254}
    # code += "v15 = svzip2_f32(acegikmo7,bdfhjlnp7);\n" # { 15, 31, 47, 63, 79, 95,111,127,143,159,175,191,207,223,239,255}
    # code += "}\n"

    # code += "void gather16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){ transpose16x16(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15); };\n"
    # code += "void scatter16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){ transpose16x16(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15); };\n"
    code
  end
end

class FuncCall
  def convert_to_code_a64fx(conversion_type)
    retval = @name + "("
    retval += "#{$current_predicate}," if @name != "table"
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
  def convert_to_code_a64fx(conversion_type,predicate=$current_predicate)
    predicate = "svptrue_b#{$min_element_size}()" if $current_predicate == nil

    if @operator != :array && @operator != :func
      if [:eq,:neq,:gt,:ge,:lt,:le].index(@operator)
        type = get_type_suffix_a64fx(@lop.get_type)
      else
        type = get_type_suffix_a64fx(self.get_type)
      end
    end
    
    case @operator
    when :uminus then
      retval="svneg_#{type}_z(" + predicate + "," + @lop.convert_to_code(conversion_type) + ")"
    when :plus  then
      retval="svadd_#{type}_z(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :minus then
      retval="svsub_#{type}_z(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :mult  then
      retval="svmul_#{type}_z(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :div   then
      retval="svdiv_#{type}_z(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :lt    then
      retval="svcmplt_#{type}(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :le    then
      retval="svcmple_#{type}(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :gt    then
      retval="svcmpgt_#{type}(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :ge    then
      retval="svcmpge_#{type}(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :eq    then
      retval="svcmpeq_#{type}(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :neq   then
      retval="svcmpne_#{type}(" + predicate + ","+@lop.convert_to_code(conversion_type)+","+@rop.convert_to_code(conversion_type)+")"
    when :not  then
      retval="svnot_b_z(" + predicate + "," + @lop.convert_to_code(conversion_type) + ")"
    when :land  then
      retval="svand_b_z(" + predicate + "," + @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
    when :lor  then
      retval="svorr_b_z(" + predicate + "," + @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
    when :landnot  then
      retval="svbic_b_z(" + predicate + "," + @lop.convert_to_code(conversion_type) + "," + @rop.convert_to_code(conversion_type) + ")"
    when :dot   then
      if (@rop == "x" || @rop == "y" || @rop == "z" || @rop == "w") && @lop.type =~ /vec/
        retval=@lop.convert_to_code(conversion_type)+"."
        retval += ["v0","v1","v2","v3"][["x","y","z","w"].index(@rop)]
        #@rop.convert_to_code(conversion_type)
      elsif @rop == "v0" || @rop == "v1" || @rop == "v2" || @rop == "v3"
        retval=@lop.convert_to_code(conversion_type)+"." + @rop
        #@rop.convert_to_code(conversion_type)
      else
        #retval = "#{@lop.convert_to_code(conversion_type)}.#{@rop.convert_to_code(conversion_type)}"
        retval = "svdup_n_#{type}(#{@lop.convert_to_code(conversion_type)}.#{@rop.convert_to_code(conversion_type)})"
      end
    when :func  then
      retval=@lop+"("+@rop.convert_to_code(conversion_type)+")"
    when :array
      retval = @lop+"["+@rop.convert_to_code(conversion_type)+"]"
    else
      abort "error: unsupported operator #{@operator} for a64fx"
    end
    retval
  end
end

class TableDecl
  def convert_to_code_a64fx(conversion_type)
    ret = ""
    nelem = get_num_elem(@type,conversion_type)
    nreg = (@table.vals.length/nelem)
    abort "size of table must be multiple of # of element of SIMD register (#{nelem})" if @table.vals.length%nelem != 0 || nreg < 1

    abort "table length must be < #{4*nelem}" if nreg > 4

    tmp = self.dup
    tmp.name += "_tmp"
    ret += tmp.convert_to_code("reference")
    if nreg == 1
      ret += "#{get_declare_type_a64fx(@type)} #{@name};\n"
      ret += "#{name} = svld1_f32(svptrue_b32(),#{tmp.name});\n"
    else
      type = "#{@type}vec#{nreg}"
      ret += "#{get_declare_type_a64fx(type)}x#{nreg} #{@name};\n"
      for i in 0...nreg
        ret += "#{name}.v#{i} = svld1_f32(svptrue_b32(),#{tmp.name} + #{nelem*i});\n"
      end
    end
    ret
  end
end

class FloatingPoint
  def convert_to_code_a64fx(conversion_type)
    case @type
    when "F64"
      "svdup_n_f64(#{@val})"
    when "F32"
      "svdup_n_f32(#{@val})"
    when "F16"
      "svdup_n_f16(#{@val})"
    end
  end
end

class IntegerValue
  def convert_to_code_a64fx(conversion_type)
    case @type
    when "S64"
      "svdup_n_s64(#{@val})"
    when "S32"
      "svdup_n_s32(#{@val})"
    when "S16"
      "svdup_n_s16(#{@val})"
    when "U64"
      "svdup_n_u64(#{@val})"
    when "U32"
      "svdup_n_u32(#{@val})"
    when "U16"
      "svdup_n_u16(#{@val})"
    end
  end
end

class String
  def convert_to_code_a64fx(conversion_type,h=$varhash)
    name = get_name(self)
    #abort "error: undefined reference to #{name} in convert_to_code_a64fx of String"
    s = self
    if h[name] != nil
      iotype = h[name][0]
      if iotype == "MEMBER"
        type = self.get_type
        case type
        when "F64"
          s = "svdup_n_f64("+self+")"
        when "F32"
          s = "svdup_n_f32("+self+")"
        when "S64"
          s = "svdup_n_s64("+self+")"
        when "S32"
          s = "svdup_n_s32("+self+")"
        when "U64"
          s = "svdup_n_u64("+self+")"
        when "U32"
          s = "svdup_n_u32("+self+")"
        end
      end
    end
    s
  end
end

def aos2soa_simd(fvars,coversion_type,h=$varhash)
  ret = []

  # count bytes of EPI and FORCE member variable
  epi_tot, epi_max_byte_size, is_epi_uniform = count_class_member("EPI")
  force_tot, force_max_byte_size, is_force_uniform = count_class_member("FORCE")

  tot = [epi_tot,force_tot]
  max_byte_size = [epi_max_byte_size,force_max_byte_size]
  is_uniform = [is_epi_uniform,is_force_uniform]
  nelems = [epi_tot/epi_max_byte_size, force_tot/force_max_byte_size]
  iotypes = ["EPI","FORCE"]

  # declare and load EPI, FORCE and local variables
  ["EPI","FORCE"].each{ |io|
    nelem  = nelems[iotypes.index(io)]
    offset = (tot[iotypes.index(io)] + max_byte_size[iotypes.index(io)] - 1) / max_byte_size[iotypes.index(io)]
    if is_uniform[iotypes.index(io)] && nelem <= 4
      #if false # structure pattern is not supported yet
      vname = ["__fkg_tmp_epi","__fkg_tmp_force"][iotypes.index(io)]
      type = "F#{max_byte_size[iotypes.index(io)]}"
      type += "vec#{nelem}" if nelem > 1
      ret += [Declaration.new([type,vname])]
      ret += [StructureLoad.new([vname,PointerOf.new([type,Expression.new([:array,get_iotype_array(io),"i",type])]),nelem,type])]

      fmem = []
      h.each{ |v|
        if v[1][0] == io
          fmem += [v[0]] 
          ret += [Declaration.new([v[1][1],v[0]])]
        end
      }
      count = 0
      h.each{ |v|
        if v[1][0] == io
          fvars.each{ |name|
            if name == v[0]
              type_single = get_single_element_type(type)
              vdim = get_vector_dim(v[1][1])
              if vdim > 1
                for i in 0...vdim
                  dest = Expression.new([:dot,v[0],["x","y","z","w"][i],type_single])
                  src  = Expression.new([:dot,vname,["v0","v1","v2","v3"][i+count],type_single])
                  ret += [Statement.new([dest,src,type_single])]
                end
                count += vdim
              else
                ret += [Statement.new([v[0],Expression.new([:dot,vname,"v#{count}",type_single]),type_single])]
                count += 1
              end
            end
          }
        end
      }
    else # is not uniform
      fvars.each{ |v|
        name = v
        iotype   = h[v][0]
        type     = h[v][1]
        fdpsname = h[v][2]
        modifier = h[v][3]
        if iotype == io
          ret += [Declaration.new([type,name])]
          if modifier == "local"
            if type =~ /vec/
              vdim = get_vector_dim(type)
              tmp_type = get_single_element_type(type)
              for i in 0...vdim
                dst = Expression.new([:dot,name,["x","y","z","w"][i],tmp_type])
                src = PointerOf.new([tmp_type,Expression.new([:array,"#{name}_tmp_"+["x","y","z","w"][i],"i",tmp_type])])
                ret += [LoadState.new([dst,src,tmp_type])]
              end
            else
              ret += [LoadState.new([name,PointerOf.new([type,Expression.new([:array,"#{name}_tmp","i",type])]),type])]
            end
          else
            tmp_type = get_single_element_type(type)
            offset_tmp = offset * max_byte_size[iotypes.index(io)] / byte_count(tmp_type)
            if type =~ /vec/
              vdim = get_vector_dim(type)
              for i in 0...vdim
                dst = Expression.new([:dot,name,["x","y","z","w"][i],tmp_type])
                src = PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),["x","y","z","w"][i],tmp_type])])
                ret += [GatherLoad.new([dst,src,"0","#{offset_tmp}",tmp_type])]
              end
            else
              ret += [GatherLoad.new([name,PointerOf.new([type,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type])]),"0","#{offset_tmp}",type])]
            end
          end
        end
      }
    end
  }
  ret
end

def soa2aos_simd(fvars,h=$varhash)
  ret = []
  # count bytes of EPI and FORCE member variable
  tot, max_byte_size, is_uniform = count_class_member("FORCE")

  nelem = tot/max_byte_size
  io = "FORCE"

  # declare and load EPI, FORCE and local variables
  offset = (tot + max_byte_size - 1) / max_byte_size
  if is_uniform && nelem <= 4 # structure store does not work with FCCpx
    vname = "__fkg_tmp_force"
    count = 0
    fvars.each{ |v|
      name     = v
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      modifier = h[v][3]
      if iotype == io
        if type =~ /vec/
          ret += [Statement.new([Expression.new([:dot,vname,"v#{count+0}",type]),Expression.new([:dot,v,"x"]),type])]
          ret += [Statement.new([Expression.new([:dot,vname,"v#{count+1}",type]),Expression.new([:dot,v,"y"]),type])]
          ret += [Statement.new([Expression.new([:dot,vname,"v#{count+2}",type]),Expression.new([:dot,v,"z"]),type])]
          count += 3
        else
          ret += [Statement.new([Expression.new([:dot,vname,"v#{count}",type]),v,type])]
          count += 1
        end
      end
    }
    type = "F#{max_byte_size}"
    ret += [StructureStore.new([PointerOf.new([type,Expression.new([:array,get_iotype_array(io),"i",type])]),vname,nelem,type])]
  else # is not uniform
    fvars.each{ |v|
      name = v
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      modifier = h[v][3]
      if iotype == io
        tmp_type = type.delete("vec")
        tmp_offset = offset * max_byte_size / byte_count(tmp_type)
        if type =~ /vec/
          ret += [ScatterStore.new([PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"x",tmp_type])]),Expression.new([:dot,name,"x",tmp_type]),0,"#{tmp_offset}",tmp_type])]
          ret += [ScatterStore.new([PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"y",tmp_type])]),Expression.new([:dot,name,"y",tmp_type]),0,"#{tmp_offset}",tmp_type])]
          ret += [ScatterStore.new([PointerOf.new([type,Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type]),"z",tmp_type])]),Expression.new([:dot,name,"z",tmp_type]),0,"#{tmp_offset}",tmp_type])]
        else
          ret += [ScatterStore.new([PointerOf.new([type,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i"]),fdpsname,type])]),name,0,"#{tmp_offset}",type])]
        end
      end
    }
  end
  ret
end
