require_relative 'common.rb'
require_relative "intermediate_exp_class.rb"
require_relative "disassemble_statement.rb"
require_relative "kernelparser.rb"
require_relative "expand_function.rb"
require_relative "reduce_madd.rb"
require_relative "loop_fission.rb"
require_relative "software_pipelining.rb"
require_relative "gen_hash.rb"
require_relative "kernel_body_multi_prec.rb"

require_relative "A64FX.rb"
require_relative "AVX-512.rb"
require_relative "AVX2.rb"
require_relative "CUDA.rb"

$dimension = 3
$reserved_function = ["rsqrt","sqrt","inv","max","min","madd","msub","nmadd","nmsub","table","to_int","to_uint","to_float","to_f16","to_f32","to_f64","to_f16vec","to_f32vec","to_f64vec","to_f16vec2","to_f32vec2","to_f64vec2","to_f16vec3","to_f32vec3","to_f64vec3","to_f16vec4","to_f32vec4","to_f64vec4","to_s16","to_s32","to_s64","to_u16","to_u32","to_u64"]
$tmp_val_count = 0

$iotypes = ["EPI","EPJ","FORCE","MEMBER","TABLE"]
$types = ["F64","F32","F16","S64","S32","S16","U64","U32","U16","F64vec","F32vec","F16vec","F64mat","F32mat","F16mat"]
$float_scalar_types = ["F64","F32","F16"]
$modifiers = ["local"]

def accumulate_related_variable(orig,vars,h)
  ret = []
  #p orig
  if vars != nil
    ret = vars
    vars.each{ |v|
      if  h[v] != nil
        if !($varhash[v][0] =~ /(EPI|EPJ|FORCE)/)
          ret += accumulate_related_variable(v,h[v][0],h) if (v != orig && h[v][1] == true)
        end
        h[v][1] = false
      end
    }
  end
  ret
end

def generate_related_map2(fs,ss,h=$varhash)
  tmp = []
  fs.each{ |f|
    tmp += [f]
    #p f
    ss.reverse_each{ |s|
      var = get_name(s)
      if isStatement(s) && h[var][3] != "local"
        if tmp.find(var)
          tmp += s.expression.get_related_variable
        elsif s.class == CoditionalBranch
          tmp += s.get_related_variable
        elsif s.class == IfElseState
          p s
          abort
        end
      end
    }
    tmp.sort!.uniq!
    #p tmp
  }
  #abort
  tmp
end

def generate_related_map(fs,ss,h=$varhash)
  tmp = []
  ss.each{ |s|
    #p get_name(s),h[get_name(s)] if isStatement(s)
    if isStatement(s) && h[get_name(s)][3] != "local"
      exp = s.expression.get_related_variable
      tmp += [[get_name(s),exp]] if exp != []
    elsif s.class == ConditionalBranch
      tmp += s.get_related_variable
    end
  }
  tmp.sort!
  #tmp.sort_by{ |v| v[0]}
  tmp.uniq!

  tmp2 = Hash.new()
  tmp.each{ |v|
    if tmp2[v[0]] == nil
      tmp2[v[0]] = [v[1],true]
    else
      tmp2[v[0]][0] += v[1]
    end
  }

  tmp3 = []
  fs.each{ |f|
    tmp3 += accumulate_related_variable(f,tmp2[f][0],tmp2)
  }
  tmp3 +=  accumulate_related_variable("",tmp2[""][0],tmp2) if tmp2[""] != nil

  ss.each{ |s|
    tmp3 += s.expression.get_related_variable if s.class == IfElseState && s.expression != nil
    if s.class == ConditionalBranch
      c = s.get_cond_related_variable
      tmp3 += c
      c.each{ |v|
        #tmp3 += tmp2[v][0] if tmp2[v] != nil
      }
    end
  }

  ret = tmp3.sort.uniq

  ret
end

def generate_force_related_map(ss,h=$varhash,fvars=[])
  ss.reverse_each{ |s|
    if isStatement(s)
      name = get_name(s)
      next if h[name][3] == "local"
      if h[name][0] == "FORCE"
        fvars += [name]
      end
      if h[name][0] == "FORCE" || fvars.index(name)
        rexps = s.expression.get_related_variable
        rexps.each{ |rexp|
          if rexp.class == Expression && rexp.operator == :dot
            fvars += [rexp.lop]
          else
            fvars += [rexp]
          end
        }
      end
    elsif s.class == ConditionalBranch
      s.conditions.zip(s.bodies){ |c,b|
        fvars += generate_force_related_map(b,h,fvars)
        fvars += c.get_related_variable
      }
    end
  }
  ret = fvars.sort.uniq
  ret
end

class Kernelprogram
  def check_references(h = $varhash)
    ref_list = []
    h.each{ |v|
      ref_list.push(v[0])
    }
    @statements.each{ |s|
      ref_list.push(get_name(s))
      vars = s.get_related_variable
      vars.each{ |v|
        if !ref_list.index(v)
          lineno = get_lineno(v)
          line   = get_line(v)
          message = "line #{lineno} : #{line}"
          for i in 1..message.index(v)
            message += " "
          end
          message += "^\n"
          message = "error : undefined reference to \"#{v}\" in check_references\n" + message
          abort message
        end
      }
    }

    @functions.each{ |f|
      ref_list = []
      f.decl.vars.each{ |v|
        ref_list.push(v)
      }
      f.statements.each{ |s|
        ref_list.push(get_name(s))
        vars = s.get_related_variable
        vars.each{ |v|
          if !ref_list.index(v)
            lineno = get_lineno(v)
            line   = get_line(v)
            message = "line #{lineno} : #{line}"
            for i in 1..message.index(v)
              message += " "
            end
            message += "^\n"
            message = "error : undefined reference to \"#{v}\"\n" + message
            abort message
          end
        }
      }
    }
  end
  def generate_hash(kerneltype)
    #    print "print kernel\n"
    #p self
    #p $funchash
    #process_iodecl($varhash)
    #process_funcdecl($funchash)
    #p $varhash
    @statements.each{|s|
      #p s
      if isStatement(s)
        if s.name.class == Expression
          if s.name.operator == :dot
            if ["x","y","z"].index(s.name.rop)
              type = s.get_type
              s.type = type + "vec"
              s.add_to_varhash
              s.type = type
            elsif $varhash[s.name.lop + "_" + s.name.rop] != nil
            else
              abort  "left value must be vector or scalar variable"
            end
          else
            abort "only vector expression is allowed for left value"
          end
        else
          s.get_type
          s.add_to_varhash
        end
      elsif s.class == TableDecl
        #p s
        s.add_to_varhash
      end
    }
    #reserved variables
    ["ni","nj","i","j","jj","jjj"].each{ |x|
      $varhash[x] = [nil,"S32",nil]
    }

    #accumhash
    @statements.each{|s|
      if isStatement(s)
        head = get_name(s)
        tail = get_tail(s)
        if $varhash[head][0] == "FORCE"
          type = $varhash[head][1]
          $accumhash[head] = Array.new if $accumhash[head] == nil
          if s.op != nil
            if tail == nil
              get_vector_elements(type).each_with_index{ |dim,i|
                $accumhash[head][i] = s.op
              }
            else
              $accumhash[head][["x","y","z","w"].index(tail)] = s.op
            end
          else
            if s.expression.class == Expression
              operator = s.expression.operator
            elsif s.expression.class == MADD
              case s.expression.operator
              when :madd
                operator = :plus
              when :msub
                operator = :minus
              else
                p s
                abort "unsupported accumulate operator #{s.expression.operator}"
              end
            elsif s.expression.class == FuncCall
              case s.expression.name
              when "max"
                operator = "max"
              when "min"
                operator = "min"
              else
                p s
                abort "unsupported accumulate fuction #{s.expression.name}"
              end
            end
            if tail == nil
              get_vector_elements(type).each_with_index{ |dim,i|
                $accumhash[head][i] = operator
              }
            else
              $accumhash[head][["x","y","z","w"].index(tail)] = operator
            end
          end
        end
      end
    }
    #p $funchash
  end
  def print_function(conversion_type)
    @functions.each{ |f|
      print f.convert_to_code(conversion_type)
    }
  end
  def print_statements(conversion_type)
    @statements.each{|s|
      #p s
      print s.convert_to_code(conversion_type)+"\n"
    }
  end
 
  def split_coefficients(orig)
    #p "split_coefficients",orig
    ret = []
    exp = orig.dup
    if !isLeaf(orig)
      if orig.operator == :mult
        lt = orig.lop.get_type
        rt = orig.rop.get_type
        if ((lt =~ /vec/) && !(rt =~ /vec/))
          #p "split candidate lop: #{orig.rop} with #{orig.lop}"
          if !isLeaf(orig.rop)
            tmp_name = add_new_tmpvar(rt)
            tmp = Statement.new([tmp_name,Expression.new([orig.rop.operator,orig.rop.lop,orig.rop.rop])])
            tmp.type = rt
            tmp.expression.type = orig.rop.get_type
            tmp.add_to_varhash
            orig = Expression.new([orig.operator,orig.lop,tmp_name,orig.type])
            ret += split_coefficients(orig.lop)
            ret.push(tmp)
          end
        elsif (!(lt =~ /vec/) && (rt =~ /vec/))
          #p "split candidate lop: #{orig.lop} with #{orig.rop}"
          if !isLeaf(orig.lop)
            tmp_name = add_new_tmpvar(lt) #"__fkg_tmp#{$tmp_val_count}"
            tmp = Statement.new([tmp_name,Expression.new([orig.lop.operator,orig.lop.lop,orig.lop.rop])])
            tmp.type = lt
            tmp.expression.type = orig.lop.get_type
            tmp.add_to_varhash
            #orig.lop = tmp_name
            orig = Expression.new([orig.operator,tmp_name,orig.rop,orig.type])
            ret += split_coefficients(orig.rop)
            ret.push(tmp)
          end
        else
          ret += split_coefficients(orig.lop)
          ret += split_coefficients(orig.rop) if orig.rop != nil
        end
      else
        ret += split_coefficients(orig.lop)
        ret += split_coefficients(orig.rop) if orig.rop != nil
      end
    else
      #if orig.class == FuncCall
      #tmp_name = "__fkg_tmp#{$tmp_val_count}"
      #tmp = Statement.new([tmp_name,orig])
      #tmp.type = orig.get_type
      #tmp.expression.type = orig.get_type
      #tmp.add_to_varhash
      #p orig.get_type
      #p tmp
      #orig = tmp_name
      #ret.push(tmp)
      #$tmp_val_count += 1
      #end
    end
    #p "return value:", ret, orig
    ret
  end

  def vector_to_scalar(exp,dim)
    #p "vector_to_scalar:"
    ret = exp.dup
    if !isLeaf(ret)
      #p ret
      ret.type = ret.type.delete("vec") if isVector(ret)
      ret.lop = vector_to_scalar(ret.lop,dim)
      ret.rop = vector_to_scalar(ret.rop,dim)
    else
      if isVector(ret)
        ret = Expression.new([:dot,ret,dim])
        ret.type = exp.get_type.delete("vec")
        #ret.type = ret.type.delete("vec")
      end
    end
    ret
  end
  
  def split_vector_expression(orig,dim = $dimension)
    ret = []
    if isVector(orig)
      if !(orig.expression.class == FuncCall)
        type = orig.get_type.split("vec")[0]
        ["x","y","z"].each{ |d|
          val = Expression.new([:dot,orig.name,d,type])
          exp = vector_to_scalar(orig.expression,d)
          #p exp
          tmp = Statement.new([val,exp,type,orig.op])
          tmp.type = type
          ret.push(tmp)
        }
      else
        abort "function returns vector value is not arrowed (must be expanded before splitting)"
        ret.push(orig)
      end
    else
      ret = [orig]
    end
    ret
  end
  
  def expand_vector_statement(orig)
    ret = []
    val = orig.name
    exp = orig.expression
    if exp.get_type =~ /vec/
      ret += split_coefficients(exp)
      ret += split_vector_expression(orig)
    else
      ret += [orig]
    end
    ret
  end

  def expand_tree
    @statements.each{ |s|
      if isStatement(s)
        s.expand_tree
      end
    }
    new_s = []
    @statements.each{ |orig|
      if isStatement(orig)
        expand_vector_statement(orig).each{ |s|
          s.reduce_madd_recursive
          s.reduce_negate_recursive
          new_s.push(s)
        }
      else
        new_s.push(orig)
      end
    }
    @statements = new_s
  end

  def calc_max_predicate_count(ss)
    ret = 1
    count = 1
    ss.each{ |s|
      if s.class == IfElseState
        if s.operator == :if
          count += 1
        elsif s.operator == :endif
          count -= 1
        end
        ret = count if count > ret
      elsif s.class == ConditionalBranch
        tmp = s.calc_max_predicate_count
        ret = tmp if ret < tmp
      end
    }
    ret
  end
  
  def kernel_class_def(conversion_type)
    code = ""
    code += "struct #{$kernel_name}{\n"
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += "PIKG::" + type + " " + name + ";\n"
      end
    }
    code += "#{$kernel_name}(){}\n"
    includeMemberVar = false
    tmp = "#{$kernel_name}("
    member_count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        tmp += "," if member_count > 0
        tmp += "PIKG::" + type + " " + name
        member_count = member_count+1
        includeMemberVar = true
      end
    }
    tmp += ")"

    member_count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        tmp += ":" if member_count == 0
        tmp += "," if member_count > 0
        tmp += name + "(" + name +")"
        member_count = member_count+1
      end
    }
    tmp += "{}\n"
    code += tmp if includeMemberVar

    code += "void initialize("
    count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        code += "," if count > 0
        name = v[0]
        type = v[1][1]
        code += "PIKG::" + type + " " + name +"_"
        count = count + 1
      end
    }
    code += "){\n"

    count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        code += name + " = " + name + "_;\n"
        count = count + 1
      end
    }
    code += "}\n"

    code += "int kernel_id = 0;\n"
    code += "void operator()"
    code += "(const #{$epi_name}* __restrict__ epi,const int ni,const #{$epj_name}* __restrict__ epj,const int nj,#{$force_name}* __restrict__ force,const int kernel_select = 1){\n"
    code
  end

  def reserved_func_def(conversion_type)
    code = ""
    # code += "template<typename Tret,typename Top>\n"
    # code += "Tret rsqrt(Top op){ return (Tret)1.0/std::sqrt(op); }\n"
    # code += "template<typename Tret,typename Top>\n"
    # code += "Tret sqrt(Top op){ return std::sqrt(op); }\n"
    # code += "template<typename Tret,typename Top>\n"
    # code += "Tret inv(Top op){ return 1.0/op; }\n"
    # code += "template<typename Tret,typename Ta,typename Tb>\n"
    # code += "Tret max(Ta a,Tb b){ return std::max(a,b);}\n"
    # code += "template<typename Tret,typename Ta,typename Tb>\n"
    # code += "Tret min(Ta a,Tb b){ return std::min(a,b);}\n"
    code += "PIKG::F64 rsqrt(PIKG::F64 op){ return 1.0/std::sqrt(op); }\n"
    code += "PIKG::F64 sqrt(PIKG::F64 op){ return std::sqrt(op); }\n"
    code += "PIKG::F64 inv(PIKG::F64 op){ return 1.0/op; }\n"
    code += "PIKG::F64 max(PIKG::F64 a,PIKG::F64 b){ return std::max(a,b);}\n"
    code += "PIKG::F64 min(PIKG::F64 a,PIKG::F64 b){ return std::min(a,b);}\n"

    code += "PIKG::F32 rsqrt(PIKG::F32 op){ return 1.f/std::sqrt(op); }\n"
    code += "PIKG::F32 sqrt(PIKG::F32 op){ return std::sqrt(op); }\n"
    code += "PIKG::F32 inv(PIKG::F32 op){ return 1.f/op; }\n"

    code += "PIKG::S64 max(PIKG::S64 a,PIKG::S64 b){ return std::max(a,b);}\n"
    code += "PIKG::S64 min(PIKG::S64 a,PIKG::S64 b){ return std::min(a,b);}\n"
    code += "PIKG::S32 max(PIKG::S32 a,PIKG::S32 b){ return std::max(a,b);}\n"
    code += "PIKG::S32 min(PIKG::S32 a,PIKG::S32 b){ return std::min(a,b);}\n"

    code += "PIKG::F64 table(PIKG::F64 tab[],PIKG::S64 i){ return tab[i]; }\n"
    code += "PIKG::F32 table(PIKG::F32 tab[],PIKG::S32 i){ return tab[i]; }\n"

    code += "PIKG::F64 to_float(PIKG::U64 op){return (PIKG::F64)op;}\n"
    code += "PIKG::F32 to_float(PIKG::U32 op){return (PIKG::F32)op;}\n"
    code += "PIKG::F64 to_float(PIKG::S64 op){return (PIKG::F64)op;}\n"
    code += "PIKG::F32 to_float(PIKG::S32 op){return (PIKG::F32)op;}\n"
    code += "PIKG::S64   to_int(PIKG::F64 op){return (PIKG::S64)op;}\n"
    code += "PIKG::S32   to_int(PIKG::F32 op){return (PIKG::S32)op;}\n"
    code += "PIKG::U64  to_uint(PIKG::F64 op){return (PIKG::U64)op;}\n"
    code += "PIKG::U32  to_uint(PIKG::F32 op){return (PIKG::U32)op;}\n"

    
    code += "template<typename T> PIKG::F64 to_f64(const T& op){return (PIKG::F64)op;}\n"
    code += "template<typename T> PIKG::F32 to_f32(const T& op){return (PIKG::F32)op;}\n"
    #code += "template<typename T> PIKG::F16 to_f16(const T& op){return (PIKG::F16)op;}\n"
    code += "template<typename T> PIKG::S64 to_s64(const T& op){return (PIKG::S64)op;}\n"
    code += "template<typename T> PIKG::S32 to_s32(const T& op){return (PIKG::S32)op;}\n"
    #code += "template<typename T> PIKG::S16 to_s16(const T& op){return (PIKG::S16)op;}\n"
    code += "template<typename T> PIKG::U64 to_u64(const T& op){return (PIKG::U64)op;}\n"
    code += "template<typename T> PIKG::U32 to_u32(const T& op){return (PIKG::U32)op;}\n"
    #code += "template<typename T> PIKG::U16 to_u16(const T& op){return (PIKG::U16)op;}\n"

    #code += "PIKG::F64 table(const PIKG::F64 tab[],const PIKG::U32 index){ return tab[index];);\n"
    #code += "PIKG::F32 table(const PIKG::F32 tab[],const PIKG::U32 index){ return tab[index];);\n"
    case conversion_type
    when /A64FX/ then
      code += self.reserved_func_def_a64fx(conversion_type)
    when /AVX2/ then
      code += self.reserved_func_def_avx2(conversion_type)
    when /AVX-512/ then
      code += self.reserved_func_def_avx512(conversion_type)
    end
    code
  end

  def split_accum_statement(ss,h = $varhash)
    j_list = []
    h.each{ |v|
      iotype = v[1][0]
      if iotype == "EPJ"
        name = v[0]
        j_list += [name]
      end
    }

    init = []
    body = []
    accum = []
    ss.each{|s|
      if isStatement(s)
        name = get_name(s.name)
        if h[name][3] != "local"
          if s.expression.isJRelated(j_list)
            j_list += [name] if j_list.index(name) == nil
            body += [s]
          else
            v = h[name]
            iotype = v[0]
            case iotype
            when "FORCE"
              accum += [s]
            when nil
              init  += [s]
            else
              body  += [s]
            end
          end
        end
      else
        body += [s]
      end
    }
    [init,body,accum]
  end
  
  def aos2soa(fvars,conversion_type,h=$varhash)
    ret = []
    case conversion_type
    when "reference"
      h.each{|v|
        iotype = v[1][0]
        if iotype == "EPI" || iotype == "FORCE"
          name     = v[0]
          type     = v[1][1]
          modifier = v[1][3]
          ret += [Declaration.new([type,name])]
          if modifier == "local"
            if type =~ /vec/
              get_vector_elements(type).each{ |dim|
                ret += [Statement.new([Expression.new([:dot,name,dim]),Expression.new([:array,"#{name}_tmp_"+dim,get_io_index(iotype),type]),type])]
              }
            else
              ret += [Statement.new([name,Expression.new([:array,"#{name}_tmp",get_io_index(iotype),type]),type])]
            end
          elsif fvars.index(name)
            fdpsname = v[1][2]
            ret += [Statement.new([name,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),type])]
          end
        end
      }
    when /(A64FX|AVX)/
      ret += aos2soa_simd(fvars,h)
    else
      abort "error: unsupported conversion_type #{conversion_type} in AoS2SoA"
    end
    ret
  end

  def soa2aos(fvars,conversion_type,h=$varhash)
    ret = []
    case conversion_type
    when "reference"
      h.each{|v|
        iotype = v[1][0]
        if iotype == "FORCE"
          name = v[0]
          type = v[1][1]
          fdpsname = v[1][2]
          if fvars.index(name)
            if type =~ /vec/
              ret += [StoreState.new([Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),"x"]),Expression.new([:dot,name,"x"]),type])]
              ret += [StoreState.new([Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),"y"]),Expression.new([:dot,name,"y"]),type])]
              ret += [StoreState.new([Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),"z"]),Expression.new([:dot,name,"z"]),type])]
            else
              ret += [StoreState.new([Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype),type]),fdpsname,type]),name,type])]
            end
          end
        end
      }
    when /(A64FX|AVX)/
      ret += soa2aos_simd(fvars,h)
    else
      abort "error: unsupported conversion_type #{conversion_type} in SoA2AoS"
    end
    ret
  end

  def generate_iloop_begin(conversion_type,max_size_epi,istart = 0)
        
    case conversion_type
    when "reference"
      ret = Loop.new(["i",istart,"ni",1,[]])
    when  /A64FX/
      nelem = get_simd_width(conversion_type) / max_size_epi
      ret = Loop.new(["i",istart,"((ni+#{nelem-1})/#{nelem})*#{nelem}","#{nelem}",["#{$current_predicate}=svwhilelt_b32_s32(i,ni);"]])
    when  /AVX/
      nelem = get_simd_width(conversion_type) / max_size_epi
      ret = Loop.new(["i",istart,"(ni/#{nelem})*#{nelem}","#{nelem}",[]])
    end
    ret
  end
  
  def kernel_body(conversion_type,istart=0,h=$varhash)
    code = ""
    if conversion_type =~ /(A64FX|AVX)/
      $pg_count = 0
      $current_predicate = "pg#{$pg_count}"
      $max_pg_count = calc_max_predicate_count(@statements)
      if conversion_type =~ /A64FX/
        for i in 0...$max_pg_count
          code += "svbool_t pg#{i};\n"
        end
      elsif conversion_type =~ /AVX2/
      elsif conversion_type =~ /AVX-512/
      end
    end

    ret = []
    ret.push(NonSimdDecl.new(["S32","i"])) if istart != nil
    ret.push(NonSimdDecl.new(["S32","j"]))

    if istart != nil
      # declare "local" variables
      h.each{ |v|
        modify = v[1][3]
        if modify == "local"
          name = v[0]
          iotype = v[1][0]
          type   = v[1][1]
          array_size = ["ni","nj"][["EPI","EPJ"].index(iotype)]
          if type =~ /vec/
            type = type.delete("vec")
            ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_x",array_size,type])])]
            ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_y",array_size,type])])]
            ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_z",array_size,type])])]
          else
            ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp",array_size,type])])]
          end
        end
      }

      # calc or copy local variables from EPI, EPJ, or FORCE
      @statements.each{ |s|
        name = get_name(s) if s.class == Statement
        if name != nil && h[name][3] == "local"
          tail = get_tail(s.name)
          iotype = h[name][0]
          type   = h[name][1]
          exp = s.expression
          index = get_io_index(iotype)
          loop_tmp = Loop.new([index,"0","n#{index}",1,[]])
          if tail != nil
            new_name = Expression.new([:array,"#{name}_tmp_#{tail}",index,type])
          else
            new_name = Expression.new([:array,"#{name}_tmp",index,type]) #"#{name}_tmp[i]"
          end
          new_exp = exp.replace_fdpsname_recursive(h)
          loop_tmp.statements += [Statement.new([new_name,new_exp])]
          ret += [loop_tmp]
        end
      }
    end
    ret.each{ |s|
      code += s.convert_to_code("reference")
    }
    ret = []

    # declare and load TABLE variable
    tmp = []
    @statements.each{ |s|
      if s.class == TableDecl
        ret.push(s)
        tmp.push(s)
      end
    }
    tmp.each{ |s|
      @statements.delete(s)
    }

    # declare temporal variables
    @statements.each{ |s|
      if s.class == Statement || s.class == ConditionalBranch
        ret += s.declare_temporal_var
      end
    }
    # check maximum size of EPI variable
    max_size_epi = -1
    h.each{ |v|
      iotype = v[1][0]
      if iotype == "EPI"
        type = v[1][1]
        size = get_single_data_size(type)
        max_size_epi = size if size > max_size_epi
      end
    }
    max_size_epi = 64 if max_size_epi > 64
    # building i loop
    accum_init,ss,accum_finalize = split_accum_statement(@statements)
    fvars = generate_force_related_map(ss)
    iloop = generate_iloop_begin(conversion_type,max_size_epi,istart)

    # load EPI and FORCE variable
    iloop.statements += aos2soa(fvars,conversion_type)

    accum_init.each{|s|
      iloop.statements += [s]
    }

    iloop.statements += [NonSimdState.new(["j","0"])]
    if $strip_mining != nil
      warn "strip mining is applied"
      loop_fission_vars = find_loop_fission_load_store_vars(ss)
      fission_count = 0
      tmpvars = []
      loop_fission_vars.each{ |vs|
        vs[0].each{ |v|
          tmpvars += [v]
        }
      }
      tmpvars.uniq.each{ |v|
        iotype = $varhash[v][0]
        type = $varhash[v][1]
        if type =~ /vec/
          type = type.delete("vec")
          iloop.statements += [NonSimdDecl.new([type,"#{v}_tmp_x[#{nsimd}*#{$strip_mining}]"])]
          iloop.statements += [NonSimdDecl.new([type,"#{v}_tmp_y[#{nsimd}*#{$strip_mining}]"])]
          iloop.statements += [NonSimdDecl.new([type,"#{v}_tmp_z[#{nsimd}*#{$strip_mining}]"])]
        else
          iloop.statements += [NonSimdDecl.new([type,"#{v}_tmp[#{nsimd}*#{$strip_mining}]"])]
        end
      }

      jloop = Loop.new(["j",nil,"(nj/#{$strip_mining})*#{$strip_mining}","#{$strip_mining}",[NonSimdDecl.new(["S32","jj"])]])
      jjloop = Loop.new(["jj","0","#{$strip_mining}",1,[]])
      jjj = NonSimdExp.new([:plus,"j","jj"])
      first_loop = true
      ss.each{ |s|
        $unroll_stage = s.option[0].to_i if s.class == Pragma && s.name == "unroll"

        if (s.class == Pragma && s.name == "statement" && s.option == ["loop_fission_point"])|| first_loop
          if !first_loop
            loop_fission_vars[fission_count][0].each{ |v|
              type = h[v][1]
              if type =~ /vec/
                jjloop.statements += [StoreState.new([PointerOf.new([type,Expression.new([:array,"#{v}_tmp_x",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),Expression.new([:dot,"#{v}","x"]),type])]
                jjloop.statements += [StoreState.new([PointerOf.new([type,Expression.new([:array,"#{v}_tmp_y",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),Expression.new([:dot,"#{v}","y"]),type])]
                jjloop.statements += [StoreState.new([PointerOf.new([type,Expression.new([:array,"#{v}_tmp_z",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),Expression.new([:dot,"#{v}","z"]),type])]
              else
                jjloop.statements += [StoreState.new([PointerOf.new([type,Expression.new([:array,"#{v}_tmp",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),"#{v}",type])]
              end
            }
            #jjloop = software_pipelining(jjloop) if $swpl_stage > 1
            jjloop = loop_unroll(jjloop,$unroll_stage) if $unroll_stage > 1
            jloop.statements += [jjloop.dup]
          end
          first_loop = false
          jjloop = Loop.new(["jj","0","#{$strip_mining}",1,[]])
          loop_fission_vars[fission_count][1].each{ |v|
            iotype = h[v][0]
            type   = h[v][1]
            if iotype == "declared"
              if type =~ /vec/
                jjloop.statements += [LoadState.new([Expression.new([:dot,"#{v}","x"]),PointerOf.new([type,Expression.new([:array,"#{v}_tmp_x",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),type])]
                jjloop.statements += [LoadState.new([Expression.new([:dot,"#{v}","y"]),PointerOf.new([type,Expression.new([:array,"#{v}_tmp_y",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),type])]
                jjloop.statements += [LoadState.new([Expression.new([:dot,"#{v}","z"]),PointerOf.new([type,Expression.new([:array,"#{v}_tmp_z",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),type])]
              else
                jjloop.statements += [LoadState.new(["#{v}",PointerOf.new([type,Expression.new([:array,"#{v}_tmp",NonSimdExp.new([:mult,"#{nsimd}","jj","S32"])])]),type])]
              end
            elsif iotype == "EPJ"
              name     = v
              fdpsname = h[v][2]
              modifier = h[v][3]

              jjloop.statements += [Declaration.new([type,name])]
              case modifier
              when "local"
                if type =~ /vec/
                  jjloop.statements += [Duplicate.new([Expression.new([:dot,name,"x"]),Expression.new([:array,"#{name}_tmp_x",jjj,type]),type])]
                  jjloop.statements += [Duplicate.new([Expression.new([:dot,name,"y"]),Expression.new([:array,"#{name}_tmp_y",jjj,type]),type])]
                  jjloop.statements += [Duplicate.new([Expression.new([:dot,name,"z"]),Expression.new([:array,"#{name}_tmp_z",jjj,type]),type])] 
                else
                  jjloop.statements += [Duplicate.new([name,Expression.new([:array,"#{name}_tmp",jjj,type]),type])]
                end
              else
                if type =~ /vec/
                  stype = type.delete("vec")
                  jjloop.statements += [Statement.new([Expression.new([:dot,name,"x"]),Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),jjj]),fdpsname,type]),"x"]),stype])]
                  jjloop.statements += [Statement.new([Expression.new([:dot,name,"y"]),Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),jjj]),fdpsname,type]),"y"]),stype])]
                  jjloop.statements += [Statement.new([Expression.new([:dot,name,"z"]),Expression.new([:dot,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),jjj]),fdpsname,type]),"z"]),stype])]
                else
                  jjloop.statements += [Statement.new([name,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),jjj]),fdpsname,type]),type])]
                end
              end
            elsif fvars.index(name)
              fdpsname = v[1][2]
              jjloop.statements += [Duplicate.new([name,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype)]),fdpsname,type]),type])]
            end
          }
          fission_count += 1
        end # s.class != Pragma
        jjloop.statements += [s] if (!isStatement(s) || h[get_name(s)][3] != "local") && s.class != Pragma
      } # ss.each
      #jjloop = software_pipelining(jjloop) if $swpl_stage > 1
      jjloop = loop_unroll(jjloop,$unroll_stage) if $unroll_stage > 1
      jloop.statements += [jjloop.dup]

      iloop.statements += [jloop]
    end # strip_mining

    # tail j loop
    jloop = Loop.new(["j",nil,"nj",1,[]])
    fvars.each{|v|
      iotype = h[v][0]
      if iotype == "EPJ"
        name     = v
        type     = h[v][1]
        modifier = h[v][3]
        jloop.statements += [Declaration.new([type,name])]
        if modifier == "local"
          if type =~ /vec/
            jloop.statements += [Duplicate.new([Expression.new([:dot,name,"x"]),Expression.new([:array,"#{name}_tmp_x",get_io_index(iotype),type]),type])]
            jloop.statements += [Duplicate.new([Expression.new([:dot,name,"y"]),Expression.new([:array,"#{name}_tmp_y",get_io_index(iotype),type]),type])]
            jloop.statements += [Duplicate.new([Expression.new([:dot,name,"z"]),Expression.new([:array,"#{name}_tmp_z",get_io_index(iotype),type]),type])]
          else
            jloop.statements += [Duplicate.new([name,Expression.new([:array,"#{name}_tmp",get_io_index(iotype),type]),type])]
          end
        elsif fvars.index(name)
          fdpsname = h[v][2]
          jloop.statements += [Statement.new([name,Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),get_io_index(iotype)]),fdpsname,type]),type])]
        end
      end
    }
    ss.each{|s|
      jloop.statements += [s] if s.class != Pragma
    }
    iloop.statements += [jloop]

    accum_finalize.each{|s|
      iloop.statements += [s]
    }
    iloop.statements += soa2aos(fvars,conversion_type)

    ret += [iloop]
    ret.each{ |s|
      code += s.convert_to_code(conversion_type)
    }

    # tail i loop
    if conversion_type =~ /AVX/
      h.each{|v|
        iotype   = v[1][0]
        v[1][0] = nil if iotype == "declared"
      }
      code += "{\n"
      code += kernel_body("reference",nil)
      code += "}\n"
    end
    
    code
  end
  def generate_optimized_code(conversion_type,output=$output_file)
    code = "#include<pikg_vector.hpp>\n"
    code +="#include<cmath>\n"
    code +="#include<limits>\n"
    code +="#include<chrono>\n"
    code += $additional_text if $additional_text != nil
    code += "\n"
    case conversion_type
    when /A64FX/
      code += "#include <arm_sve.h>\n"
    when /AVX2/
      code += "#include <pikg_avx2.hpp>\n"
    when /AVX-512/
      code += "#include <pikg_avx512.hpp>\n"
    end

    if $c_interface_impl
      struct_list = ["EPI"]
      struct_list.push("EPJ") if $epi_name != $epj_name
      struct_list.push("FORCE") if $epi_name != $force_name && $epj_name != $force_name
      struct_list.zip([$epi_name,$epj_name,$force_name]){ |c,n|
        next if n.index("PS::")
        code +="struct #{n}{\n"
        $varhash.each{|v|
          iotype = get_iotype_from_hash(v);
          modifier = get_modifier_from_hash(v)
          if iotype == c && modifier == nil
            type = get_type_from_hash(v)
            fdpsname = get_fdpsname_from_hash(v)
            code += "PIKG::#{type} #{fdpsname};\n"
          end
        }
        code += "};\n"
      }
    end

    code += kernel_class_def(conversion_type)
    code += kernel_body(conversion_type)
    code += "}\n"
    code += reserved_func_def(conversion_type)
    code += "}; // kernel functor definition\n"


    if $c_interface_impl
      code += "#{$kernel_name} __pikg_#{$interface_name};\n"

      code += "extern \"C\"{\n"
      code += "void #{$initializer_name}("
      count = 0
      $varhash.each{|v|
        iotype = v[1][0]
        if iotype == "MEMBER"
          code += "," if count > 0
          name = v[0]
          type = v[1][1]
          code += "PIKG::" + type + " " + name +"_"
          count = count + 1
        end
      }
      code += "){\n"
      code += "__pikg_#{$interface_name}.initialize("
      count = 0
      $varhash.each{|v|
        iotype = v[1][0]
        if iotype == "MEMBER"
          code += "," if count > 0
          name = v[0]
          type = v[1][1]
          code += name +"_"
          count = count + 1
        end
      }
      code += ");\n"
      code += "} // intialize_#{$interface_name}\n"

      code += "void #{$interface_name}(const #{$epi_name}* epi, const PIKG::S32 ni, const #{$epj_name}* epj,const PIKG::S32 nj, #{$force_name}* force){\n"
      code += "__pikg_#{$interface_name}(epi,ni,epj,nj,force);\n"
      code += "}\n"
      code += "} // extern \"C\"\n"
    end


    if output == nil
      print code
    else
      File.open(output, mode = 'w'){ |f|
        f.write(code)
      }
    end
  end


  def c_interface_type_decl(type)
    case type
    when "F64"
      "double"
    when "F32"
      "float"
    when "S64"
      "long long int"
    when "S32"
      "int"
    when "U64"
      "unsigned long long int"
    when "U32"
      "unsigned int"
    when "F64vec"
      "pikg_f64vec"
    when "F64vec2"
      "pikg_f64vec2"
    when "F64vec3"
      "pikg_f64vec3"
    when "F64vec4"
      "pikg_f64vec4"
    when "F32vec"
      "pikg_f32vec"
    when "F32vec2"
      "pikg_f32vec2"
    when "F32vec3"
      "pikg_f32vec3"
    when "F32vec4"
      "pikg_f32vec4"
    else
      abort "unsupported type for c_interface_type_decl"
    end
  end

  def convert_fdps_defined_type(name)
    case name
    when "PS::SPJMonopole"
      "fdps_spj_monopole"
    when "PS::SPJQuadrupole"
      "fdps_spj_quadrupole"
    when "PS::SPJMonopoleGeometricCenter"
      "fdps_spj_monopole_geomcen"
    when "PS::SPJDipoleGeometricCenter"
      "fdps_spj_dipole_geomcen"
    when "PS::SPJQuadrupoleGeometricCenter"
      "fdps_spj_quadrupole_geomcen"
    when "PS::SPJMonopoleScatter"
      "fdps_spj_monopole_scatter"
    when "PS::SPJQuadrupoleScatter"
      "fdps_spj_quadrupole_scatter"
    when "PS::SPJMonopoleSymmetry"
      "fdps_spj_quadrupole_symmetry"
    when "PS::SPJQuadrupoleSymmetry"
      "fdps_spj_quadrupole_symmetry"
    when "PS::SPJMonopoleCutoff"
      "fdps_spj_monopole_cutoff"
    else
      name
    end
  end

  def generate_prototype_decl_file()
    code =  ""
    code += "#ifndef H_PROTOTYPE_DECL_#{$interface_name.capitalize}\n"
    code += "#define H_PROTOTYPE_DECL_#{$interface_name.capitalize}\n"
    code += "#include <pikg_vector.h>\n"
    #code += "void #{$interface_name}(const struct #{$epi_name}*, const int, const struct #{convert_fdps_defined_type($epj_name)}*,const int, struct #{$force_name}*);\n"
    code += "void #{$interface_name}(void*, const int, void*,const int, void*);\n"
    code += "void #{$initializer_name}("
    count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        code += "," if count > 0
        name = v[0]
        type = v[1][1]
        code += c_interface_type_decl(type) + " " + name +"_"
        count = count + 1
      end
    }
    code += ");\n"
    code += "#endif\n"
    File.open($prototype_decl_name, mode = 'w'){ |f|
      f.write(code)
    }
  end

  def fortran_type(type)
    case type
    when "F64"
      "real(kind=c_double)"
    when "F32"
      "real(kind=c_float)"
    when "S64"
      "integer(kind=c_long_long)"
    when "S32"
      "integer(kind=c_int)"
    when "F64vec"
      "pikg_f64vec3"
    when "F64vec2"
      "pikg_f64vec2"
    when "F64vec3"
      "pikg_f64vec3"
    when "F32vec"
      "pikg_f32vec3"
    when "F32vec2"
      "pikg_f32vec2"
    when "F32vec3"
      "pikg_f32vec"
    else
      abort "unsupported fortran type"
    end
  end

  def generate_fortran_module()
    indent = "   "
    code =  ""
    code += "module #{$module_name}\n"
    code += indent + "use, intrinsic :: iso_c_binding\n"
    code += "\n"
    code += indent + "interface\n"
    code += "\n"
    code += indent * 2 + "subroutine #{$interface_name}(epi, ni, epj, nj, force) &\n"
    code += indent * 4 + "bind(c,name='#{$interface_name}')\n"
    code += indent * 3 + "use, intrinsic :: iso_c_binding\n"
    code += indent * 3 + "implicit none\n"
    code += indent * 3 + "integer(kind=c_int), value :: ni,nj\n"
    code += indent * 3 + "type(c_ptr), value :: epi, epj, force\n"
    code += indent * 2 + "end subroutine\n"
    code += "\n"
    code += indent * 2 + "subroutine #{$initializer_name} &\n"
    code += indent * 3 + "( &\n"
    count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        code += "," if count > 0
        name = v[0]
        code += indent * 3 + name + " & \n"
        count = count + 1
      end
    }
    code += indent * 3 + ") &\n"
    code += indent * 4 + "bind(c,name='#{$initializer_name}')\n"
    code += indent * 3 + "use, intrinsic :: iso_c_binding\n"
    code += indent * 3 + "implicit none\n"
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += indent * 3 + fortran_type(type) + ", value :: " + name + "\n"
      end
    }
    code += indent * 2 + "end subroutine\n\n" 
    code += indent + "end interface\n\n"
    code += "end module #{$module_name}\n"
    $module_file_name = $module_name + ".F90"
    File.open($module_file_name, mode = 'w'){ |f|
      f.write(code)
    }
  end

  def make_conditional_branch_block(h = $varhash)
    fvars = generate_force_related_map(@statements)
    @statements,dummy = make_conditional_branch_block_recursive3(@statements,h)
  end

  def make_conditional_branch_block_recursive3(ss,h = $varhash,related_vars = [])
    ret = Array.new
    nest_level = 0
    cbb = nil

    ss.reverse_each{ |s|
      if nest_level == 0
        if s.class == IfElseState && s.operator == :endif
          cbb = ConditionalBranch.new([[],[]])
          cbb.init_new_body
        else
          if isStatement(s)
            name = get_name(s)
            if related_vars.index(name) || (h[name] != nil && h[name][0] == "FORCE")
              related_vars.push(name)
              related_vars += s.expression.get_related_variable
            end
          elsif s.class == Pragma
            # do nothing
          elsif s.class == TableDecl
            # do nothing
          else
            related_vars += s.expression.get_related_variable
          end
          related_vars.sort!.uniq!
          ret.push(s)
        end
      elsif nest_level == 1
        if s.class == IfElseState
          case s.operator
          when :if
            cbb.push_condition(s)
          when :elsif
            cbb.push_condition(s)
            cbb.init_new_body
          when :else
            cbb.push_condition(s)
            cbb.init_new_body
          when :endif
            cbb.push_body(s)
          end
        else
          cbb.push_body(s)
        end
      else
        cbb.push_body(s)
      end

      nest_level += 1 if s.class == IfElseState && s.operator == :endif
      if s.class == IfElseState && s.operator == :if
        if nest_level == 1
          new_b = []
          new_c = []
          related_vars_tmp = nil
          cbb.bodies.reverse_each{ |b|
            b.reverse!
            tmp,related_vars_tmp = make_conditional_branch_block_recursive3(b,h,related_vars)
            new_b.push(tmp)
          }
          cbb.conditions.reverse_each{ |c|
            new_c.push(c)
          }
          cbb.bodies = new_b
          cbb.conditions.reverse!

          ret.push(cbb)

          # make tmp var hash
          tmp_name_hash = Hash.new()
          merge_state = []
          cbb.bodies.each{ |bss|
            bss.each { |bs|
              if isStatement(bs) && bs.expression.class != Merge
                name = get_name(bs)
                tail = get_tail(bs)
                if related_vars.index(name) || h[name][0] == "FORCE"
                  tmp_name_hash[name] = add_new_tmpvar(bs.type) if tmp_name_hash[name] == nil
                  if ["x","y","z","w"].index(tail)
                    src = Expression.new([:dot,tmp_name_hash[name],tail,bs.type])
                    dst = Expression.new([:dot,name,tail,bs.type])
                    bss.push(Statement.new([dst,Merge.new([src,dst,bs.type]),bs.type]))
                  else
                    bss.push(Statement.new([name,Merge.new([tmp_name_hash[name],name,bs.type]),bs.type]))
                  end
                end
              end
            }
          }
          related_vars = related_vars_tmp if related_vars_tmp != nil
          #p tmp_name_hash
          cbb.bodies.each{  |bss|
            computed_list = []
            replaced_list = []
            bss.each{ |bs|
              if isStatement(bs) && bs.expression.class != Merge
                name = get_name(bs)
                bs.replace_name(name,tmp_name_hash[name]) if tmp_name_hash[name] != nil
                bs.expression.replace_by_list(computed_list,replaced_list)
                if tmp_name_hash[name] != nil
                  computed_list.push(bs.name)
                  replaced_list.push(tmp_name_hash[name])
                end
              end
            }
          }
          cbb = nil
        end
        nest_level -= 1
      end
      abort "nest_level < 0" if nest_level < 0
    }
    [ret.reverse, related_vars]
  end

  def make_conditional_branch_block_recursive2(ss,h = $varhash,related_vars = [])
    #p h
    new_s = []
    nest_level = 0
    cbb = ConditionalBranch.new([[],[]])

    ss.reverse_each { |s|
      if s.class == IfElseState
        nest_level -= 1 if s.operator == :if
        nest_level += 1 if s.operator == :endif
        abort "nest level < 0" if nest_level < 0
      end
      #p s
      if nest_level == 0
        related_vars.push(get_name(s)) if isStatement(s)
        related_vars += s.expression.get_related_variable

        #p related_vars
        related_vars.sort!.uniq!
        # push new conditional branch block
        if s.class == IfElseState
          abort "operator #{s.operator} must be :if" if s.operator != :if
          cbb.push_condition(s)

          new_b = []
          new_c = []
          cbb.bodies.reverse_each{ |b|
            b.reverse!
            new_b.push(make_conditional_branch_block_recursive2(b,h,related_vars))
          }
          cbb.conditions.reverse_each{ |c|
            new_c.push(c)
          }
          cbb.bodies = new_b
          cbb.conditions.reverse!

         # make tmp var hash
          tmp_name_hash = Hash.new()
          merge_state = []
          tmp_cbb = ConditionalBranch.new([[],[]])
          tmp_cbb.conditions = cbb.conditions
          cbb.bodies.each{ |bss|
            tmp_bss = Array.new
            bss.each { |bs|
              tmp_bss.push(bs)
              if isStatement(bs) && bs.expression.class != Merge
                name = get_name(bs)
                tail = get_tail(bs)
                if related_vars.find(){ |n| n == name } || h[name][0] == "FORCE"
                  tmp_name_hash[name] = add_new_tmpvar(bs.type) if tmp_name_hash[name] == nil
                  tmp_bss = [Statement.new([tmp_name_hash[name],name,h[name][1],nil])] + tmp_bss
                  #p bss
                  if ["x","y","z","w"].index(tail)
                    dst = Expression.new([:dot,tmp_name_hash[name],tail,bs.type])
                    src = Expression.new([:dot,name,tail,bs.type])
                    tmp_bss.push(Statement.new([dst,Merge.new([src,dst,bs.type]),bs.type]))
                  else
                    tmp_bss.push(Statement.new([name,Merge.new([name,tmp_name_hash[name],bs.type]),bs.type]))
                  end
                end
              end
            }
            tmp_cbb.bodies.push(tmp_bss)
          }
          cbb = tmp_cbb
          new_s.push(cbb)
          #p tmp_name_hash
          cbb.bodies.each{  |bss|
            computed_list = []
            replaced_list = []
            bss.each{ |bs|
              if isStatement(bs) && bs.expression.class != Merge
                name = get_name(bs)
                #bs.replace_name(name,tmp_name_hash[name]) if tmp_name_hash[name] != nil
                #bs.expression.replace_by_list(computed_list,replaced_list)
                if tmp_name_hash[name] != nil
                  computed_list.push(bs.name)
                  replaced_list.push(tmp_name_hash[name])
                end
              end
            }
          }
          cbb = ConditionalBranch.new([[],[]])
        else
          new_s.push(s)
        end

      elsif nest_level > 0
        if s.class == IfElseState
          case s.operator
          when :if
            cbb.push_condition(s)
          when :elsif
            cbb.push_condition(s)
            cbb.init_new_body
          when :endif
            cbb.init_new_body
          end
        else
          cbb.push_body(s)
        end
      end
    }
    abort "ConditionalBranch is not terminated" if nest_level > 0
    new_s.reverse
  end
end

class String
  def get_type(h = $varhash)
    propaties=h[self]
    if propaties
      propaties[1]
    else
      p h
      p h[self]
      p self
      nil[0]
      abort "undefined reference to #{self} at get_type of String class"
    end
  end

  def get_related_variable
    [self]
  end
  def isJRelated(list)
    if list.index(self)
      true
    else
      false
    end
  end

  def fusion_iotag(iotag)
    abort "fusion_iotag for #{iotag} failed" if self == iotag
    self
  end

  def replace_recursive(orig,replaced)
    if self == orig
      replaced
    else
      self
    end
  end
  def replace_fdpsname_recursive(h=$varhash)
    name = self.dup
    ret = name
    return ret if h[name] == nil
    iotype = h[name][0]
    type = h[name][1]
    fdpsname = h[name][2]
    if iotype != nil && fdpsname != nil
      op = "i" if iotype == "EPI" || iotype == "FORCE"
      op = "j" if iotype == "EPJ"
      ret = Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),op]),fdpsname,type])
    end
    ret
  end

  def replace_by_list(name_list,replaced_list)
    ret = self
    name_list.zip(replaced_list){ |n,r|
      ret = r if self == n
    }
    ret
  end
  
  def convert_to_code(conversion_type="reference",h=$varhash)
    s=self
    #print "convert to code #{s}\n"
    #p $varhash[self]
    case conversion_type
    when "reference" then
      # nothing to do
    when /A64FX/
      s = self.convert_to_code_a64fx(h)
    when /AVX2/
      s = self.convert_to_code_avx2(h)
    when /AVX-512/
      s = self.convert_to_code_avx512(h)
    end
    #    print "result=", s, "\n"
    s
  end
end



parser=KernelParser.new
$kernel_name="Kernel"
$epi_name="EPI"
$epj_name="EPJ"
$force_name="FORCE"
$conversion_type = "reference"
$swpl_stage = 1
$unroll_stage = 1
$output_file = "kernel.hpp"
while true
  opt = ARGV.shift
  break if opt == nil
  case opt
  when "-i"
    filename = ARGV.shift
    warn "input file: #{filename}\n"
  when "--input"
    filename = ARGV.shift
    warn "input file: #{filename}\n"
  when "--kernel-name"
    $kernel_name = ARGV.shift
    warn "kernel name: #{$kernel_name}\n"
  when "--epi-name"
    $epi_name = ARGV.shift
    warn "epi name: #{$epi_name}\n"
  when "--epj-name"
    $epj_name = ARGV.shift
    warn "epj name: #{$epj_name}\n"
  when "--spj-name"
    $spj_name = ARGV.shift
    warn "spj name: #{$spj_name}\n"
  when "--force-name"
    $force_name = ARGV.shift
    warn "force name: #{$force_name}\n"
  when "--conversion-type"
    $conversion_type = ARGV.shift
    warn "conversion type: #{$conversion_type}\n"
  when "--strip-mining"
    $strip_mining = ARGV.shift.to_i
    warn "strip mining size: #{$strip_mining}\n"
  when "--software-pipelining"
    $swpl_stage = ARGV.shift.to_i
    abort "software pipelining is not available"
    warn "software pipelining stage: #{$swpl_stage}\n"
  when "--unroll"
    $unroll_stage = ARGV.shift.to_i
    warn "software pipelining stage: #{$unroll_stage}\n"
  when "-o"
    $output_file = ARGV.shift
    warn "output file name: #{$output_file}\n"
  when "--output"
    $output_file = ARGV.shift
    warn "output file name: #{$output_file}\n"
  when "--initializer-name"
    $initializer_name = ARGV.shift
    warn "kernel initializer name: #{$initializer_name}"
  when "--additional-text"
     $additional_text = ARGV.shift
    warn "additional text: #{$additional_text}"
  when "--multiwalk"
    $is_multi_walk = true
    warn "multi walk mode\n"
  when "--c-interface"
    $c_interface = true
    warn "c interface mode on\n"
    if !ARGV.empty? && ARGV[0][0] != "-"
      $prototype_decl_name = ARGV.shift
      warn "prototype decl file name: #{$prototype_decl_name}"
    end
  when "--fortran-interface"
    $fortran_interface = true
    warn "fortran interface mode on\n"
    $module_name = ARGV.shift
    warn "module name: #{$module_name}"
  when "--class-file"
    $epi_file = ARGV.shift
    if ARGV[0][0] != "-"
      warn "epi class file: #{$epi_file}"
      $epj_file = ARGV.shift
      warn "epj class file: #{$epj_file}"
      $force_file = ARGV.shift
      warn "force class file: #{$force_file}"
    else
      $epj_file = $force_file = $epi_file
      warn "class file: #{$epi_file}"
    end
  when "--version"
    warn "pikg version 0.1d"
    abort
  when "--help"
    help_message = "available options:\n"
    help_message += "--input | -i file_name : input file name\n"
    help_message += "--output | -o file_name : output file name (default: kernel.hpp)\n"
    help_message += "--kernel-name kernel_name : kernel name (default: Kernel)\n"
    help_message += "--epi-name epi_name : c++ class name of EPI (default: EPI)\n"
    help_message += "--epj-name epj_name : c++ class name of EPJ (default: EPJ)\n"
    help_message += "--force-name force_name : c++ class name of FORCE (default: FORCE)\n"
    help_message += "--class-file epi_file_name [epj_file_name force_file_name] : file name which include EPI/EPJ/FORCE class definition. If this option is not enabled, you need to declare alias of all the member variable of EPI/EPJ/FORCE in input file. If you specified a sigle file_name, the file_name is used for all EPI/EPJ/FORCE class. (default: nil)\n"
    help_message += "--conversion-type type : target architecture (reference, AVX2, AVX-512, or A64FX)\n"
    help_message += "--c-iterface [file_name] : enable c-interface mode. header file name of prototype definition can be specified.\n"
    help_message += "--fortran-iterface module_name : enable fortran-interface mode. c-interface mode is automatically enabled. specify kernel module name as module_name. module is output to module_name + \".F90\"\n"
    help_message += "--initializer-name [func_name] : function name of kernel initializer for c-interface (default: kernel_name + \"_initialize\")\n"
    help_message += "--version : show version info\n"
    help_message += "--help : show this help message\n"
    
    warn help_message
    abort
  else
    abort "error: unsupported option #{opt}"
  end
end
#abort "output file must be specified with --output option" if $output_file == nil

if $c_interface
  $c_interface_impl = true
  $c_interface_decl = true
end
if $fortran_interface
  $c_interface_impl = true
  $c_interface_decl = false
end

if $c_interface_impl || $c_interface_decl
  $interface_name = $kernel_name
  $kernel_name = $kernel_name + "_"
  $initializer_name = $interface_name + "_initialize" if $initializer_name == nil
end

if $c_interface_decl
  if $prototype_decl_name == nil
    tmp =  $output_file.split ('.')

    if tmp.length > 1
      tmp[-1] = "h"
      $prototype_decl_name = tmp.join('.')
    else
      $prototype_decl_name = tmp.join + ".h"
    end
    warn "prototype decl file name: #{$prototype_decl_name}"
  end
end

src = ""
program=parser.parse(filename)
$varhash = Hash.new
$funchash = Hash.new
$accumhash = Hash.new
program.process_funcdecl($funchash)
if $epi_file != nil && $epj_file != nil && $force_file != nil
  ["EPI","EPJ","FORCE"].each{ |iotype|
    class_file, class_name = [[$epi_file,$epi_name],[$epi_file,$epj_name],[$force_file,$force_name]][["EPI","EPJ","FORCE"].index(iotype)]
    if $fortran_interface
      program.generate_hash_from_fortran(class_file,iotype,class_name,$varhash)
    elsif $c_interface
      program.generate_hash_from_c(class_file,iotype,class_name,$varhash)
    else
      program.generate_hash_from_cpp(class_file,iotype,class_name,$varhash)
    end
  }
  program.generate_alias($varhash)
else
  program.process_iodecl($varhash)
end
program.check_references($varhash)

program.generate_hash("noconversion")
program.expand_function
program.expand_tree
program.make_conditional_branch_block
program.disassemble_statement
program.generate_prototype_decl_file if $c_interface_decl
program.generate_fortran_module if $fortran_interface
if $is_multi_walk
  program.generate_optimized_code_multi_walk($conversion_type);
else
  case $conversion_type
  when "CUDA"
    program.generate_optimized_cuda_kernel($conversion_type)
#  when "A64FX"
#    program.generate_optimized_code($conversion_type)
  else
    program.generate_optimized_code_multi_prec($conversion_type)
  end
end

__END__
