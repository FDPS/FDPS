def get_initial_value(op,type)
  if op == :plus || op == :minus
    case type
    when /(S|U)/
    #IntegerValue.new("0")
      "0"
    when "F64"
    #FloatingPoint.new("0.0")
      "0.0"
    when "F32"
      #FloatingPoint.new("0.0f")
      "0.0f"
    when "F16"
      #FloatingPoint.new("0.0h")
      "0.0h"
    end
  elsif op == :mult || op == :div
    case type
    when /(S|U)/
      #IntegerValue.new("1")
      "1"
    when "F64"
      #FloatingPoint.new("1.0")
      "1.0"
    when "F32"
      #FloatingPoint.new("1.0f")
      "1.0f"
    when "F16"
      #FloatingPoint.new("1.0h")
      "1.0h"
    end
  elsif op == "max"
    case type
    when "S64"
      "std::numeric_limits<int64_t>::lowest()"
    when "S32"
      "std::numeric_limits<int32_t>::lowest()"
    when "S16"
      "std::numeric_limits<int16_t>::lowest()"
    when "U64"
      "std::numeric_limits<uint64_t>::lowest()"
    when "U32"
      "std::numeric_limits<uint32_t>::lowest()"
    when "U16"
      "std::numeric_limits<uint16_t>::lowest()"
    when "F64"
      "std::numeric_limits<double>::lowest()"
    when "F32"
      "std::numeric_limits<float>::lowest()"
    when "F16"
      "std::numeric_limits<float16_t>::lowest()"
    end
  elsif op == "min"
    case type
    when "S64"
      "std::numeric_limits<int64_t>::max()"
    when "S32"
      "std::numeric_limits<int32_t>::max()"
    when "S16"
      "std::numeric_limits<int16_t>::max()"
    when "U64"
      "std::numeric_limits<uint64_t>::max()"
    when "U32"
      "std::numeric_limits<uint32_t>::max()"
    when "U16"
      "std::numeric_limits<uint16_t>::max()"
    when "F64"
      "std::numeric_limits<double>::max()"
    when "F32"
      "std::numeric_limits<float>::max()"
    when "F16"
      "std::numeric_limits<float16_t>::max()"
    else
      abort "error: unsupported type #{type} for get_initial_value"
    end
  else
    abort "unsupported accumulate operator #{op} for get_initial_value"
  end
end

def get_iotype_from_hash(v)
  v[1][0]
end
def get_type_from_hash(v)
  v[1][1]
end
def get_fdpsname_from_hash(v)
  v[1][2]
end
def get_modifier_from_hash(v)
  v[1][3]
end

def get_vector_elements(type)
  case type
  when /F(16|32|64)vec2/
    ["x","y"]
  when /F(16|32|64)vec(3)?/
    ["x","y","z"]
  when /F(16|32|64)vec4/
    ["x","y","z","w"]
  else
    [""]
  end
end

def sizeof(type)
  case type
  when /F64vec(3)?/
    64*3
  when /F32vec(3)?/
    32*3
  when /F16vec(3)?/
    16*3
  when "F64vec2"
    64*2
  when "F32vec2"
    32*2
  when "F16vec2"
    16*2
  when "F64vec4"
    64*4
  when "F32vec4"
    32*4
  when "F16vec4"
    16*4
  when /(F|S|U|B)64/
    64
  when /(F|S|U|B)32/
    32
  when /(F|S|U|B)16/
    16
  end
end

def isVector(val)
  if val.get_type =~ /vec/
    true
  else
    false
  end
end

def isVariable(val)
  if val.class == String
    true
  else
    false
  end
end

def isLeaf(val)
  ret = false
  if val.class == String || val.class == FuncCall || val.class == IntegerValue || val.class == FloatingPoint
    ret = true
  else
    if val.class == Expression
      if val.operator == :dot || val.operator == :array
        ret = true
      end
    end
  end
  ret
end

def isStatement(val)
  if val.class == Statement
    true
  else
    false
  end
end

def get_name(x)
  ret = nil
  if x.class == Expression
    if x.operator == :dot || x.operator == :array
      ret = get_name(x.lop)
    else
      abort "get_name for not :dot epxpression"
    end
  elsif x.class ==  Statement || x.class == NonSimdState
    if x.name.class == Expression
      ret = get_name(x.name)
    else
      ret = x.name
    end
  elsif x.class == Duplicate
    ret = get_name(x.name)
  elsif x.class == Declaration
    ret = get_name(x.name)
  elsif x.class == LoadState
    ret = get_name(x.dest)
  elsif x.class == StoreState
    ret = get_name(x.src)
  elsif x.class == PointerOf
    ret = get_name(x.exp)
  elsif x.class == String
    ret = x
  elsif x.class == TableDecl
    ret = get_name(x.name)
  else
    return nil
    abort "get_name is not allowed to use for #{x.class}"
  end
  ret
end

def get_tail(x)
  ret = nil
  if x.class == Expression
    if x.operator == :dot
      ret = x.rop
    else
      abort "get_name for not :dot epxpression"
    end
  elsif x.class == Statement
    if x.name.class == Expression
      if x.name.operator == :dot
        ret = x.name.rop
      else
        ret = x.name
      end
    else
      ret = nil
    end
  elsif x.class == String
    ret = nil
  else
    abort "get_tail is not allowed to use for #{x.class}"
  end
  ret
end

def get_declare_type(type,conversion_type)
  case conversion_type
  when "reference"
    case type
    when "S32"
      decl = "PIKG::S32"
    when "S64"
      decl = "PIKG::S64"
    when "U32"
      decl = "PIKG::U32"
    when "U64"
      decl = "PIKG::U64"
    when "F32"
      decl = "PIKG::F32"
    when "F64"
      decl = "PIKG::F64"
    when "F32vec"
      decl = "PIKG::F32vec"
    when "F64vec"
      decl = "PIKG::F64vec"
    else
      abort "error: unsupported type #{type} for get_declare_type"
    end
  when /A64FX/
    decl = get_declare_type_a64fx(type)
  when /AVX2/
    decl = get_declare_type_avx2(type)
  when /AVX-512/
    decl = get_declare_type_avx512(type)
  when /CUDA/
    case type
    when "S32"
      decl = "int"
    when "S64"
      decl = "long long"
    when "U32"
      decl = "unsigned int"
    when "U64"
      decl = "unsigned long long"
    when "F64"
     decl =  "double"
    when "F32"
      decl = "float"
    when "F64vec"
      decl = "double3"
    when "F32vec"
      decl = "float3"
    when "F64vec4"
      decl = "double4"
    when "F32vec4"
      decl = "float4"
    when "F64vec2"
      decl = "double2"
    when "F32vec2"
      decl = "float2"
    when "F64vec3"
      decl = "double3"
    when "F32vec3"
      decl = "float3"
    end
  else
    abort "error: unsupported conversion_type #{conversion_type} at get_declare_type"
  end
  decl
end

def get_iotype_array(iotype)
  case iotype
  when "EPI"
    "epi"
  when "EPJ"
    "epj"
  when "FORCE"
    "force"
  else
    abort "error: undefined iotype #{iotype}"
  end
end
def get_io_index(iotype)
  case iotype
  when "EPI"
    "i"
  when "EPJ"
    "j"
  when "FORCE"
    "i"
  else
    abort "error: undefined iotype #{iotype}"
  end
end


def get_data_size(type)
  case type
  when /^F64vec2$/
    128
  when /^F64vec(3)?$/
    192
  when /^F64vec4$/
    256
  when /^F32vec2$/
    64
  when /^F32vec(3)?$/
    96
  when /^F32vec4$/
    128
  when /^F16vec2$/
    32
  when /^F16vec(3)?$/
    48
  when /^F16vec4$/
    64
  when /^(F|S|U)64$/
    64
  when /^(F|S|U)32$/
    32
  when /^(F|S|U)16$/
    16
  else
    p type
    abort "unsupported type for get_data_size : #{type}"
  end
end

def get_single_data_size(type)
  case type
  when /64/
    64
  when /32/
    32
  when /16/
    16
  else
    p type
    abort "unsupported type for get_single_data_size : #{type}"
  end
end

def get_single_element_type(type)
  case type
  when /F64/
    "F64"
  when /F32/
    "F32"
  when /F16/
    "F16"
  when /S64/
    "S64"
  when /S32/
    "S32"
  when /S16/
    "S16"
  when /U64/
    "U64"
  when /U32/
    "U32"
  when /U16/
    "U16"
  else
    abort "unsupported type for get_single_element_type : #{type}"
  end
end

def get_simd_width(conversion_type)
  case conversion_type
  when "reference"
    1
  when /A64FX/
    512
  when /AVX2/
    256
  when /AVX-512/
    512
  end
end

def get_num_elem(type,conversion_type)
  if conversion_type == "reference"
    1
  else
    get_simd_width(conversion_type) / get_single_data_size(type)
  end
end

def get_byte_size(type)
  get_single_data_size(type) / 8
end

def get_vector_dim(type)
  case type
  when /vec2$/
    2
  when /vec(3)?$/
    3
  when /vec4$/
    4
  else
    1
  end
end

def add_new_tmpvar(type,h=$varhash)
  ret = "__fkg_tmp#{$tmp_val_count}"
  abort "type of #{ret} is nil at add_new_tmpvar!" if type == nil
  h[ret] = [nil,type,nil]
  $tmp_val_count += 1
  ret
end

def get_pointer_cast(type)
  case type
  when "F64" then
    cast = "(double*)"
  when "F32" then
    cast = "(float*)"
  when "S64" then
    cast = "(long long*)"
  when "S32" then
    cast = "(int*)"
  when "U64" then
    cast = "(unsigned long long*)"
  when "U32" then
    cast = "(unsigned int*)"
  when "F64vec" then
    cast = "(double*)"
  when "F32vec" then
    cast = "(float*)"
  when "F32vec3"
    cast = "(float*)"
  else
    abort "error: unsupported scalar type of #{@type} for get_pointer_cast"
  end
  cast
end

def count_class_member(io,h=$varhash)
  tot = 0
  is_uniform = true
  prev_type = nil
  max_byte_size = 0
  h.each{ |v|
    iotype = v[1][0]
    type   = v[1][1]
    modifier = v[1][3]
    if iotype == io && modifier == nil
      next if v[0] =~ /_swpl/
      prev_type = type.delete("vec") if tot == 0
      byte = byte_count(type)
      tot += byte_count(type) if modifier == nil
      byte = byte / 3 if type =~ /vec/
      max_byte_size = byte if byte > max_byte_size
      is_uniform = false if prev_type != nil && type.delete("vec") != prev_type
    end
  }
  #warn "size of #{io} member is #{max_byte_size}, # of elemennts are #{tot/max_byte_size}" if is_uniform
  [tot,max_byte_size,is_uniform]
end
