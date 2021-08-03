class Kernelprogram
  def process_iodecl(h = $varhash)
    @iodeclarations.each{|x|
      h[x.name] = [x.iotype, x.type, x.fdpsname, x.modifier]
    }
  end

  def process_funcdecl(h = $funchash)
    @functions.each{|x|
      decl = x.decl
      stmt = x.statements
      ret  = x.retval
      h[decl.name] = x
    }
  end

  def fusion_iotag
    ["EPI","EPJ","FORCE"].each{ |iotag|
      @statements.each{|s|
        s.fusion_iotag(iotag)
      }
    }
  end

  def cpp_type_to_pikg_type(type)
    case type
    when /(PS::|PIKG::)?F64vec(3)?/
      "F64vec"
    when /(PS::|PIKG::)?F64vec2/
      "F64vec2"
    when /(PS::|PIKG::)?F64vec4/
      "F64vec4"
    when /(PS::|PIKG::)?F32vec(3)?/
      "F32vec"
    when /(PS::|PIKG::)?F32vec2/
      "F32vec2"
    when /(PS::|PIKG::)?F32vec4/
      "F32vec4"
    when /(PS::|PIKG::)?F16vec(3)?/
      "F16vec"
    when /(PS::|PIKG::)?F16vec2/
      "F16vec2"
    when /(PS::|PIKG::)?F16vec4/
      "F16vec4"
    when /(double|float64_t|(PS::|PIKG::)?F64)/
      "F64"
    when /(float|float32_t|(PS::|PIKG::)?F32)/
      "F32"
    when /(half|float16_t|(PS::|PIKG::)?F16)/
      "F16"
    when /(long\s+long(\s+int)?|int64_t|(PS::|PIKG::)?S64)/
      "S64"
    when /(int|int32_t|(PS::|PIKG::)?S32)/
      "S32"
    when /(short|int16_t|(PS::|PIKG::)?S16)/
      "S16"
    when /(unsigned\s+long\s+long(\s+int)?|uint64_t|(PS::|PIKG::)?U64)/
      "U64"
    when /(unsigned(\s+int)?|uint32_t|(PS::|PIKG::)?U32)/
      "U32"
    when /(unsigned\s+short|uint16_t|(PS::|PIKG::)?U16)/
      "U16"
    else
      abort "unsupported c++ type #{type} for cpp_type_to_pikg_type"
    end
  end
  
  def generate_hash_from_cpp(filename,iotype,class_name,h = $varhash)
    decl = /(const)?\ +((PS::|PIKG::)?((F|S)(64|32|16)(vec(2|3|4)?)?)|(double|float((64|32|16)_t)?|(unsigned\ +)?(long\ +)+?int((64|32|16)_t)?))/
    ident=/[a-zA-Z_][a-zA-Z_0-9]*/
    code = String.new
    File.open(filename){ |f|
      f.each_line{ |line|
        next if line[0] == '#'
        next if line =~ /^\/\//
        code += line.chomp
      }
    }
    code = code.gsub(";","\n").gsub("{","{\n").gsub("}","}\n").gsub("public:",'').gsub("private:",'')
    #warn code

    #warn "test result:\n"
    nest_level = 0
    ignore = 0
    base_level = 1
    read = false
    code.each_line{ |line|
      next if line =~ /using/
      base_level += 1 if line =~ /namespace/
      if line =~ /(class|struct)\ +#{class_name}/
        read = true
      end
      if read
        nest_level += 1 if line =~ /\{/
        if nest_level == base_level
          if line =~ /^\s+#{decl}\s+#{ident}(\s*,\s*#{ident})*\s*\n/
            tmp = line.split(/(\s|,)/).select{ |s| s=~/(#{decl}|#{ident})/}
            type = String.new
            type = cpp_type_to_pikg_type(tmp.shift)
            if type == "const"
              type = tmp.shift
            end
            type = type.gsub("PS::","").gsub("PIKG::","")
            vars = tmp
            vars.each{ |v|
              h[iotype+"."+v] = [iotype,type,v,nil]
            }
          elsif line =~ /#{ident}\ +#{ident}(\ *,\ *#{ident})*\ *\n/
            #warn "undefined type: #{line}"
          end
        end
        nest_level -= 1 if line =~ /\}/
        read = false if nest_level == 0
      end
    }
  end

  def c_type_to_pikg_type(ctype)
    case ctype
    when "double"
      "F64"
    when "float"
      "F32"
    when "long long int"
      "S64"
    when "long long"
      "S64"
    when "int"
      "S32"
    when "unsigned long long int"
      "U64"
    when "unsigned long long"
      "U64"
    when "unsigned int"
      "U32"
    when "unsigned"
      "U32"
    when "pikg_f64vec"
      "F64vec"
    when "pikg_f64vec2"
      "F64vec2"
    when "pikg_f64vec3"
      "F64vec3"
    when "pikg_f64vec4"
      "F64vec4"
    when "pikg_f32vec"
      "F32vec"
    when "pikg_f32vec2"
      "F32vec2"
    when  "pikg_f32vec3"
      "F32vec3"
    when "pikg_f32vec4"
      "F32vec4"
    else
      abort "unsupported type for c_type_to_pikg_type \"#{ctype}\""
    end

  end
  def generate_hash_from_c(filename,iotype,class_name,h = $varhash)
    types = /((fdps|pikg)_(f|s)(64|32|16)(vec(2|3|4)?)?|double|float((64|32|16)_t)?|(unsigned\ +)?(long\ *)+?(int((64|32|16)_t)?)?)/
    decl = /(const)?\ +#{types}/
    ident=/[a-zA-Z_][a-zA-Z_0-9]*/
    code = String.new
    File.open(filename){ |f|
      f.each_line{ |line|
        next if line[0] == '#'
        next if line =~ /^\/\//
        code += line.chomp
      }
    }
    code = code.gsub(";","\n").gsub("{","{\n").gsub("}","}\n")
    #warn code

    #warn "test result:\n"
    nest_level = 0
    ignore = 0
    base_level = 1
    read = false
    h_tmp = Hash.new
    code.each_line{ |line|
      if line =~ /(typedef struct)/
        read = true
        h_tmp = Hash.new
      end
      if line =~ /\s*#{class_name}$/
        h_tmp.each{ |v|
          h[v[0]] = v[1]
        }
      end
      if read
        nest_level += 1 if line =~ /\{/
        if nest_level == base_level
          if line =~ /^\s*#{decl}\s+#{ident}(\s*,\s*#{ident})*\s*\n/
            tmp = line.split(/(\s|,)/).select{ |s| s=~/\S+/}
            type = tmp.shift
            next if type == "static"
            while tmp[0] =~ /(const|#{types})/
              tmp.shift if tmp[0] == "const"
              type = type + " " + tmp.shift
            end

            type = c_type_to_pikg_type(type)
            vars = tmp
            vars.each{ |v|
              h_tmp[iotype+"."+v] = [iotype,type,v,nil]
            }
          elsif line =~ /#{ident}\ +#{ident}(\ *,\ *#{ident})*\ *\n/
            warn "undefined type: #{line}"
          end
        end
        nest_level -= 1 if line =~ /\}/
        read = false if nest_level == 0
      end
    }
  end

  def fortran_type_to_pikg_type(type)
    case type
    when "real(kind=c_double)"
      "F64"
    when "real(kind=c_float)"
      "F32"
    when "integer(kind=c_long_long)"
      "S64"
    when "integer(kind=c_int)"
      "S32"
    when /type\((pikg|fdps)_f64vec4\)/
      "F64vec4"
    when /type\((pikg|fdps)_f64vec(3)?\)/
      "F64vec"
    when /type\((pikg|fdps)_f64vec2\)/
      "F64vec2"
    when /type\((pikg|fdps)_f32vec4\)/
      "F32vec4"
    when /type\((pikg|fdps)_f32vec(3)?\)/
      "F32vec"
    when /type\((pikg|fdps)_f32vec2\)/
      "F32vec2"
    else
      abort "unsupported fortran type #{type}"
    end
  end
  def generate_hash_from_fortran(filename,iotype,class_name,h = $varhash)
    types = /(integer\(kind\s*=\s*(c_int((8|16|32|64)_t)?|c_long(_long)?)\)|real\(kind\s*=\s*c_((long_)?double|float)\)|type\((pikg|fdps)_f(16|32|64)vec(2|3|4)?\))/i
    decl = /#{types}/
    ident=/[a-zA-Z_][a-zA-Z_0-9]*/
    attribute = /(dimension|intent|value|public)/
    fdps_tag = /\$fdps.+/
    code = String.new
    File.open(filename){ |f|
      f.each_line{ |line|
        #next if line[0] == '#'
        #next if line =~ /^\/\//
        code += line
      }
    }
    #code = code.gsub(";","\n").gsub("{","{\n").gsub("}","}\n")
    #warn code

    #warn "test result:\n"
    nest_level = 0
    ignore = 0
    base_level = 1
    read = false
    h_tmp = Hash.new
    code.each_line{ |line|
      if line =~ /\s*type.+\:\:\s+#{class_name}/
        read = true
        h_tmp = Hash.new
        next
      end
      if line =~ /\s*end\s*type\s*#{class_name}/
        h_tmp.each{ |v|
          h[v[0]] = v[1]
        }
        read = false
      end
      if read
        if line =~ /^\s*#{decl}\s*(::)?\s*#{ident}(\s*,\s*#{ident})*\s*(\!#{fdps_tag})?\n/
          type_decl = line.split(/(\s|,|::|\!)/).select{ |s| s=~/\S+/}
          type = type_decl.shift
          type_decl.shift while type_decl[0] =~ attribute
          vars = type_decl.select{ |s| s =~ ident}
          type = fortran_type_to_pikg_type(type)
          vars.each{ |v|
            break if v == "$fdps"
            h_tmp[iotype+"."+v] = [iotype,type,v,nil]
          }
        end
      end
    }
  end

  def append_iotag(statement,iotag,h=$varhash)
    ret = statement

    rexps = ret.expression.get_related_variable
    rexps.each{ |rexp|
      type = nil
      h.each{ |v|
        iotype = get_iotype_from_hash(v)
        fdpsname = get_fdpsname_from_hash(v)
        if iotype == iotag && fdpsname == rexp
          type = get_type_from_hash(v)
        end
      }
      ret.expression = ret.expression.replace_recursive(rexp,iotag+"_"+rexp) if type != nil
    }
    ret
  end
  
  def generate_alias(h=$varhash)
    fusion_iotag
    process_iodecl_alias(@iodeclarations,h)
    new_h = Hash.new
    h.each{ |v|
      replace = v[0].gsub(".","_")
      new_h[replace] = v[1]
      h.delete(v[0])
    }
    new_h.each{ |v|
      h[v[0]] = v[1]
    }
    @statements.each{ |s|
      if isStatement(s)
        lexp = get_name(s)
        type = nil
        h.each{ |v|
          iotype = get_iotype_from_hash(v)
          fdpsname = get_fdpsname_from_hash(v)
          if iotype == "FORCE" && fdpsname == lexp
            type = get_type_from_hash(v)
          end
        }
        s.replace_name(lexp,"FORCE_"+lexp) if type != nil
        s.expression.replace_recursive(lexp,"FORCE_"+lexp) if type != nil

        s = append_iotag(s,"EPI",h)
        s = append_iotag(s,"EPJ",h)
        s = append_iotag(s,"FORCE",h)
      elsif s.class == TableDecl
        s.add_to_varhash
      end
    }
    #p h
    #@statements.each{ |s|
    #  p s.convert_to_code("reference")
    #}
  end

  def process_iodecl_alias(ios,h=$varhash)
    ios.each{|x|
      if x.modifier == nil && ["EPI","EPJ","FORCE"].index(x.iotype)
        isIncluded = false
        h.each{ |v|
          fdpsname = v[1][2]
          isIncluded = true if fdpsname == x.fdpsname
        }
        abort "alias to undefined member #{x.fdpsname}" if !isIncluded
        h[x.name] = [x.iotype, x.type, x.fdpsname, "alias"]
      else
        h[x.name] = [x.iotype, x.type, x.fdpsname, x.modifier]
      end
    }
  end

end


