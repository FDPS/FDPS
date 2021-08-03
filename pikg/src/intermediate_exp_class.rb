# class definition
class Kernelprogram
  attr_accessor :iodeclarations, :functions, :statements
  def initialize(x)
    #@iodeclarations, @functions, @statements =*x
    if x[0][0].class == Iodeclaration
      @iodeclarations = x.shift
    else
      @iodeclarations = []
    end
    if x[0][0].class == Function
      @functions = x.shift
    else
      @functions = []
    end
    @statements = x.shift
  end
end

class Loop
  attr_accessor :index,:loop_beg,:loop_end,:interval,:statements,:option
  def initialize(x)
    @index,@loop_beg,@loop_end,@interval,@statements,@option = x
  end
  def convert_to_code(conversion_type="reference")
    ret = ""
    ret += "for("
    ret += "#{@index} = #{@loop_beg}" if @loop_beg  != nil
    ret += ";"
    case @option
    when :down
      ret += "#{@index} < (#{@loop_end}/#{@interval})*#{@interval};"
    when :up
      ret += "#{@index} < ((#{@loop_end} + #{@interval} - 1)/#{@interval})*#{@interval};"
    else
      ret += "#{@index} < #{@loop_end};"
    end
    if @interval == 1
      ret += "++#{@index}"
    else
      ret += "#{@index} += #{@interval}"
    end
    ret += "){\n"

    predicate = nil
    count = $pg_count
    if conversion_type == "A64FX"
      if @interval != 1 && @option == :up
        predicate = $current_predicate
        $current_predicate = "pg#{$pg_count}"
        if predicate =~ /svptrue_b/
          ret += "svbool_t #{$current_predicate} = svwhilelt_b#{$min_element_size}_s#{$min_element_size}(#{@index},#{@loop_end});\n"
        else
          ret += "svbool_t #{$current_predicate} = svand_b_z(svptrue_b#{$min_element_size}(),svwhilelt_b#{$min_element_size}_s#{$min_element_size}(#{@index},#{@loop_end}),#{predicate});\n"
        end
        $pg_count += 1
      end
    end
    statements.each{|s|
      ret += s.convert_to_code(conversion_type) + "\n"
    }
    ret += "} // loop of #{@index}\n"
    if conversion_type == "A64FX"
      if @interval != 1 && @option == :up
        $current_predicate = predicate
        $current_predicate = "svptrue_b#{$min_element_size}()" if $current_predicate == nil
      end
    end
    $pg_count = count
    ret
  end
end

class ILoop
  attr_accessor :statements
  def initialize(x)
    @statements = x
  end
  def convert_to_code(conversion_type="reference")
    ret = ""
    case conversion_type
    when "reference"
      ret += "for(i=0;i<ni;i++){\n"
    when /A64FX/
      ret += "for(i=0;i<((ni+15)/16)*16;i+=16){\n"
      ret += "#{$current_predicate} = svwhilelt_b32(i,ni);\n"
    when /AVX-512/
      ret += "for(i=0;i<(ni/16)*16;i+=16){\n"
    else
      abort "unsupporeted conversion_type (#{conversion_type}) in ILoop conversion"
    end
    statements.each{|s|
      ret += s.convert_to_code(conversion_type) + "\n"
    }
    ret += "}\n"
  end
end

class ConditionalBranch
  attr_accessor :conditions,:bodies
  def initialize(x)
    @conditions,@bodies = x
    abort "# of bodies is not correspond to # of conditions" if @conditions.length != @bodies.length
  end

  def init_new_body
    @bodies.push([])
  end
  def push_body(x)
    @bodies.last.push(x)
  end
  def push_condition(x)
    @conditions.push(x)
  end

  def get_related_variable
    ret = []
    @conditions.zip(@bodies){ |c,b|
      ret.push(["",c.expression.get_related_variable]) if c.expression != nil
      b.each{|s|
        if isStatement(s)
          exp = s.expression.get_related_variable
          #exp = [get_name(s)] if exp == []
          exp += [get_name(s)]
          ret.push([get_name(s),exp]) if exp != []
        elsif s.class == ConditionalBranch
          ret += s.get_related_variable
        end
      }
    }
    ret
  end

  def get_cond_related_variable
    ret = []
    @conditions.zip(@bodies){ |c,b|
      ret += c.expression.get_related_variable if c.expression != nil
      b.each{|s|
        ret += s.get_cond_related_variable if s.class == ConditionalBranch
      }
    }
    ret
  end

  def declare_temporal_var(h=$varhash)
    ret = []
    @bodies.each{ |b|
      b.each{ |s|
        ret += s.declare_temporal_var(h) if s.class == Statement || s.class == ConditionalBranch
      }
    }
    ret
  end

  def calc_max_predicate_count(count = 1)
    ret = count + 1
    @bodies.each{ |b|
      b.each{ |s|
        if s.class == ConditionalBranch
          tmp = s.calc_max_predicate_count(ret)
          ret = tmp if tmp > ret
        end
      }
    }
    ret
  end

  def convert_to_code(conversion_type="reference")
    ret = ""
    conditions.zip(bodies){ |c,b|
      #p [c,b]
      ret += c.convert_to_code(conversion_type) + "\n"
      b.each{ |s|
        ret += s.convert_to_code(conversion_type) + "\n"
      }
    }
    ret += IfElseState.new([:endif,nil]).convert_to_code(conversion_type)
  end
end

class Table
  attr_accessor :vals, :type
  def initialize(x)
    @vals = x
    get_type
  end

  def get_type
    @type = vals[0].get_type
    @vals.each{ |v|
      abort "type of table member must be same type #{@type}" if @type != v.get_type
      abort "type of table member must be scalar type" if v.get_type =~ /(vec|mat)/
    }
    @type
  end
end

class TableDecl
  attr_accessor :name, :table, :type
  def initialize(x)
    @name, @table = x
    @type = @table.get_type
  end

  def add_to_varhash(h=$varhash)
    #p $varhash
    if h[name] == nil
      h[name] = [nil, @type, "table"]
    end
  end

  def convert_to_code(conversion_type="reference")
    case conversion_type
    when "reference"
      ret = "PIKG::#{@type} #{@name}[] = {"
      @table.vals.each_with_index{ |v,i|
        ret += "," if i != 0
        ret += v.convert_to_code("reference")
      }
      ret += "};\n"
    when /A64FX/
      ret = convert_to_code_a64fx(conversion_type)
    when /AVX-512/
      ret = convert_to_code_avx512(conversion_type)
    end
    ret
  end

  def get_related_variable
    []
  end
end

class FloatingPoint
  attr_accessor :val, :type
  def initialize(x)
    @val = x
    @type = "F64"
    case @val
    when /f/
      @type = "F32"
    when /h/
      @type = "F16"
    end
  end

  def get_type(h=$varhash)
    @type
  end

  def fusion_iotag(iotag)
    self
  end

  def replace_fdpsname_recursive(h=$varhash)
    self
  end
  def replace_recursive(orig,replaced)
    self.dup
  end
  def replace_by_list(n,l)
  end

  def isJRelated(list)
    false
  end

  def get_related_variable
    []
  end

  def zero(type)
    case type
    when /F64/
      FloatingPoint.new("0.0")
    when /F32/
      FloatingPoint.new("0.0f")
    when /F16/
      FloatingPoint.new("0.0h")
    end
  end

  def convert_to_code(conversion_type="reference")
    case conversion_type
    when "reference"
      @val
    when /A64FX/
      convert_to_code_a64fx(conversion_type)
    when /AVX2/
      convert_to_code_avx2(conversion_type)
    when /AVX-512/
      convert_to_code_avx512(conversion_type)
    end
  end
end

class IntegerValue
  attr_accessor :val, :type
  def initialize(x)
    @val = x
    @type = "S32"
    case @val
    when /ul/
      @type = "U64"
    when /us/
      @type = "U16"
    when /l/
      @type = "S64"
    when /s/
      @type = "S16"
    when /u/
      @type = "U32"
    end
  end

  def fusion_iotag(iotag)
    self
  end
  def get_type(h = $varhash)
    @type
  end

  def replace_recursive(orig,replaced)
    self.dup
  end
  def replace_fdpsname_recursive(h=$varhash)
    self.dup
  end

  def isJRelated(list)
    false
  end

  def get_related_variable
    []
  end

  def convert_to_code(conversion_type="reference")
    case conversion_type
    when "reference"
      case @type
      when "F64"
        @val.to_s
      when "F32"
        @val.to_s
      when "F16"
        @val.to_s
      when "S64"
        @val.to_s + "l"
      when "S32"
        @val.to_s
      when "S16"
        @val.to_s.delete("s")
      when "U64"
        @val.to_s+"l"
      when "U32"
        @val.to_s.delete("u")
      when "U16"
        @val.to_s.delete("us")
      end
    when /A64FX/
      convert_to_code_a64fx(conversion_type)
    when /AVX2/
      convert_to_code_avx2(conversion_type)
    when /AVX-512/
      convert_to_code_avx512(conversion_type)
    end
  end
end

class Iodeclaration
  attr_accessor :iotype,:type,:name,:fdpsname,:modifier
  def initialize(x)
    #p x
    @iotype,@type,@name,@fdpsname,@modifier=x
  end
end

class Declaration
  attr_accessor :type, :name
  def initialize(x)
    #p x
    @type,@name = x
  end
  def append(x)
    @name += [x]
  end
  def get_type(h = $varhash)
    @type = @name.get_type(h) if @type==nil
    @type
  end
  def convert_to_code(conversion_type="reference")
    @type = @name.get_type if @type == nil
    "#{get_declare_type(@type,conversion_type)} #{@name.convert_to_code(conversion_type)};\n"
  end
end

class NonSimdDecl < Declaration
  def convert_to_code(conversion_type="reference")
    #p self
    "#{get_declare_type(@type,"reference")} #{@name.convert_to_code("reference")};\n"
  end
end

class Function
  attr_accessor :decl, :statements, :retval
  def initialize(x)
    @decl, @statements, @retval = x
  end
end

class Funcdeclaration
  attr_accessor :name, :vars, :attr
  def initialize(x)
    #p x
    @name,@vars,@attr = x
  end
end

class ReturnState
  attr_accessor :ret
  def initialize(x)
    #p x
    @ret = x
  end

  def convert_to_code(conversion_type="reference")
    retval = "return " + ret.convert_to_code(conversion_type) + ";\n"
    retval
  end
end

class Statement
  attr_accessor :name, :expression, :type, :op
  def initialize(x)
    #p x
    @name, @expression, @type, @op=x
  end
  def get_type(h = $varhash)
    @type = @expression.get_type(h) if @type==nil
    @type
  end
  def add_to_varhash(h=$varhash)
    #p $varhash
    name = get_name(@name)
    if h[name] == nil
      h[name] = [nil, @type, nil]
    end
  end

  def fusion_iotag(iotag)
    @name = @name.lop + "_" + @name.rop if @name.class == Expression && @name.lop == iotag
    @expression = @expression.fusion_iotag(iotag)
  end
  
  def get_related_variable
    [get_name(self)] + @expression.get_related_variable
  end
  def expand_inner_prod(orig)
    #p "expand_inner_prod:"
    exp = orig
    if !isLeaf(orig)
      lt = orig.lop.get_type
      if !(lt =~ /vec/)
        orig.lop = expand_inner_prod(orig.lop)
      end
      if orig.rop != nil
        rt = orig.rop.get_type
        if !(rt =~ /vec/)
          orig.rop = expand_inner_prod(orig.rop)
        end
        if lt =~ /vec/ && rt =~ /vec/
          if orig.operator == :mult
            type = lt.delete("vec")
            xx  = Expression.new([:mult,Expression.new([:dot,orig.lop,"x",type]),Expression.new([:dot,orig.rop,"x",type]),type])
            yy  = Expression.new([:mult,Expression.new([:dot,orig.lop,"y",type]),Expression.new([:dot,orig.rop,"y",type]),type])
            zz  = Expression.new([:mult,Expression.new([:dot,orig.lop,"z",type]),Expression.new([:dot,orig.rop,"z",type]),type])
            exp = Expression.new([:plus,xx,Expression.new([:plus,yy,zz,type]),type])
          else
            exp = orig
          end
        else
          exp = orig
        end
      end
    end
    exp
  end

  def expand_tree
    if @op !=nil
      @expression = Expression.new([@op,@name,@expression,@type])
      @op = nil
    end
    ret = expand_inner_prod(@expression)
    @expression = ret
  end

  def replace_name(orig,replaced)
    abort "name must be Strinng class to replace" if orig.class != String || replaced.class != String
    if @name.class == String && @name == orig
      @name = replaced
    elsif @name.class == Expression && @name.operator == :dot && @name.lop == orig
      @name.lop = replaced
    end
  end

  def declare_temporal_var(h = $varhash)
    ret = []
    name = get_name(self)
    tail = get_tail(self)
    if h[name][0] == nil && h[name][0] != "declared"
      type = @type
      type += "vec" if tail != nil
      ret += [Declaration.new([type,name])]
      h[name][0] = "declared"
    end
    ret
  end

  def convert_to_code(conversion_type="reference")
    #p self
    #p @name
    #p $varhash[@name]
    if @op == nil
      ret = @name.convert_to_code(conversion_type) + " = " + @expression.convert_to_code(conversion_type) + ";"
    else
      #warn "convert_to_code for non simple assigned operator is used"
      #p self
      ret = @name.convert_to_code(conversion_type) + " = " + Expression.new([@op,@name,@expression,@type]).convert_to_code(conversion_type) + ";"
    end
    ret
  end
end

class NonSimdState < Statement
  def convert_to_code(conversion_type="reference")
    @name.convert_to_code("reference") + " = " + @expression.convert_to_code("reference") + ";"
  end
end

class FuncCall
  attr_accessor :name, :ops, :type
  def initialize(x)
    @name, @ops = x
    @type = nil
  end
  def get_type(h = $varhash)
    ret_type = nil
    if $reserved_function.index(@name)
      type = @ops[0].get_type(h)
      #ret_type = h[@ops[0]][1]
      case @name
      when "to_uint"
        if  type =~ /64/
          ret_type = "U64"
        elsif  type =~ /32/
          ret_type = "U32"
        elsif  type =~ /16/
          ret_type = "U16"
        end
      when "to_int"
        if  type =~ /64/
          ret_type = "S64"
        elsif  type =~ /32/
          ret_type = "S32"
        elsif  type =~ /16/
          ret_type = "S16"
        else
        end
      when "to_float"
        if  type =~ /64/
          ret_type = "F64"
        elsif  type =~ /32/
          ret_type = "F32"
        elsif  type =~ /32/
          ret_type = "F16"
        end
      when "to_f16"
        ret_type = "F16"
      when "to_f32"
        ret_type = "F32"
      when "to_f64"
        ret_type = "F64"
      when "to_f16vec"
        ret_type = "F16vec"
      when "to_f32vec"
        ret_type = "F32vec"
      when "to_f64vec"
        ret_type = "F64vec"
      when "to_f16vec2"
        ret_type = "F16vec2"
      when "to_f32vec2"
        ret_type = "F32vec2"
      when "to_f64vec2"
        ret_type = "F64vec2"
      when "to_f16vec3"
        ret_type = "F16vec3"
      when "to_f32vec3"
        ret_type = "F32vec3"
      when "to_f64vec3"
        ret_type = "F64vec3"
      when "to_f16vec4"
        ret_type = "F16vec4"
      when "to_f32vec4"
        ret_type = "F32vec4"
      when "to_f64vec4"
        ret_type = "F64vec4"
      when "to_s16"
        ret_type = "S16"
      when "to_s32"
        ret_type = "S32"
      when "to_s64"
        ret_type = "S64"
      when "to_u16"
        ret_type = "U16"
      when "to_u32"
        ret_type = "U32"
      when "to_u64"
        ret_type = "U64"
      else
        ret_type = @ops[0].get_type(h)
      end
    else
      if $funchash[@name] != nil
        func = $funchash[@name]
        #p func
        abort "error: number of operands(= #{@ops.length}) for function #{@name} is wrong, #{func.decl.vars.length}" if @ops.length != func.decl.vars.length

        a = []
        func.decl.vars.zip(@ops).each{|tmp, orig|
          type = nil
          if orig.class == Expression
            type = [nil,orig.get_type,nil]
          else
            type = h[orig]
          end
          abort "error: undefined reference to #{orig} in get_type of FuncCall" if type == nil
          a += [tmp, type]
        }
        tmphash = Hash[*a]
        #p tmphash
        if func.statements != nil
          func.statements.each{ |stmt|
            #p stmt
            type = stmt.get_type(tmphash)
            if tmphash[stmt.name] == nil
              tmphash[stmt.name] = [nil,type,nil]
            end
          }
        end
        ret_type = func.retval.ret.get_type(tmphash)
        abort "error: function returns vector type is not allowed" if ret_type =~ /vec/
      else
        abort "error: undefined reference to function #{@name}"
      end
    end
    abort "error: function returns vector variable is not supported" if ret_type =~ /vec/
    abort "error: type inference failed for #{self}" if ret_type == nil
    ret_type
  end

  def fusion_iotag(iotag)
    tmp = Array.new
    @ops.each{ |op| tmp += [op.fusion_iotag(iotag)] }
    @ops = tmp
    self
  end
  
  def find_function
    [self]
  end

  
  def expand_tree
  end

  def get_related_variable
    ret = []
    @ops.each{ |op|
      ret += op.get_related_variable
    }
    ret.sort!
    ret.uniq!
    ret
  end
  
  def isJRelated(list)
    ret = false
    @ops.each{ |op|
      ret |= op.isJRelated(list)
    }
    ret
  end

  def replace_recursive(orig,replaced)
    ops = Array.new
    @ops.each{ |op|
      if orig.class == String
        ops.push(op.replace_recursive(orig,replaced))
      elsif orig.class == Expression && orig.operator == :dot
        ops.push(op.replace_recursive(orig.lop,replaced))
      else
        abort "error: #{orig} cannot be replaced with #{replaced} in FuncCall"
      end
    }
    FuncCall.new([@name,ops,@type])
  end

  def replace_fdpsname_recursive(h=$varhash)
    ops = Array.new
    @ops.each{ |op|
      ops.push(op.replace_fdpsname_recursive(h))
    }
    FuncCall.new([@name,ops,@type])
  end
  
  def replace_by_list(name_list,replaced_list)
    name_list.zip(replaced_list){ |n,r|
      self.replace_recursive(n,r)
    }
  end
  
  def convert_to_code(conversion_type="reference")
    case conversion_type
    when "reference" then
      retval = @name
      # retval += "<PIKG::" + get_type
      # @ops.each{ |op|
      #   retval += ",PIKG::" + op.get_type
      # }
      # retval += ">"
      retval += "("
      if @ops.length > 0
        @ops.each_with_index{ |op,i|
          retval += "," if i > 0
          retval += op.convert_to_code(conversion_type)
        }
      end
      retval += ")"
    when /A64FX/
      retval = self.convert_to_code_a64fx(conversion_type)
    when /AVX2/
      retval = self.convert_to_code_avx2(conversion_type)
    when /AVX-512/
      retval = self.convert_to_code_avx512(conversion_type)
    else
      abort "error: unsupported conversion_type #{conversion_type}"
    end
    retval
  end
end

class Pragma
  attr_accessor :name, :option
  def initialize(x)
    @name,@option = x
  end

  def get_related_variable
    []
  end

  def declare_temporal_var
    []
  end

  def convert_to_code(conversion_type="reference")
    ret = String.new
    if conversion_type != "reference"
      ret = "#pragma #{@name}"
      if @option != nil then
        @option.each{ |x|
          ret += " #{x}"
        }
      end
    end
    ret
  end
end

class IfElseState
  attr_accessor :operator, :expression
  def initialize(x)
    @operator, @expression = x
  end
  def get_related_variable
    ret = []
    ret += @expression.get_related_variable if @expression != nil
    ret
  end
  def fusion_iotag(iotag)
    @expression = @expression.fusion_iotag(iotag) if @operator == :if || @operator == :elsif
  end
  def convert_to_code(conversion_type="reference")
    ret = ""
    case conversion_type
    when "reference"
      case @operator
      when :if
        ret = "if(" + @expression.convert_to_code(conversion_type) + "){"
      when :elsif
        ret = "}else if(" + @expression.convert_to_code(conversion_type) + "){"
      when :else
        ret = "}else{"        
      when :endif
        ret = "}"
      else
        abort "error: undefined if operator #{@operator}"
      end
    when /(A64FX|AVX2|AVX-512)/
      if $predicate_queue == nil
        $predicate_queue = Array.new()
      end
      if $nest_queue == nil
        $nest_queue = Array.new()
      end
      if @expression != nil
        abort if !(@expression.get_type =~ /(F|S|B|U)(64|32|16)/)
        size = sizeof(@expression.get_type)
        type = "B#{size}"
      end
      ret = ""
      pg_accum = "pg_accum"
      case @operator
      when :if
        $predicate_queue.push($current_predicate)
        ret += "{\n"
        $accumulate_predicate = "pg#{$pg_count}"
        $varhash[$accumulate_predicate] = [nil,type,nil,nil]
        $pg_count += 1
        tmp_predicate = "pg#{$pg_count}"
        $varhash[tmp_predicate] = [nil,type,nil,nil]
        $pg_count += 1
        #p size
        #p $min_element_size
        nsplit = size/$min_element_size
        nsplit = 1 if conversion_type == "reference"
        pred_ops = Array.new
        accum_ops = Array.new
        for i in 0...nsplit
          exp = @expression
          suffix = String.new
          suffix = "_#{i}" if nsplit > 1
          ret += Declaration.new([type,tmp_predicate+suffix]).convert_to_code(conversion_type)
          ret += Declaration.new([type,$accumulate_predicate+suffix]).convert_to_code(conversion_type)
          exp.get_related_variable.each{ |v|
            exp = exp.replace_recursive(v,v+suffix)
          }
          ret += Statement.new([tmp_predicate+suffix,exp,type]).convert_to_code(conversion_type) + "\n"
          ret += Statement.new([$accumulate_predicate+suffix,tmp_predicate+suffix,type,nil]).convert_to_code(conversion_type) + "\n"
          pred_ops.push(tmp_predicate + suffix)
          accum_ops.push($accumulate_predicate + suffix)
        end
        if nsplit > 1
          ret += Declaration.new([type,tmp_predicate]).convert_to_code(conversion_type)
          ret += Declaration.new([type,$accumulate_predicate]).convert_to_code(conversion_type)
          ret += Statement.new([tmp_predicate,Fusion.new([pred_ops,type,"B#{$min_element_size}"])]).convert_to_code(conversion_type) + "\n" if
          ret += Statement.new([$accumulate_predicate,Fusion.new([accum_ops,type,"B#{$min_element_size}"])]).convert_to_code(conversion_type) + "\n"
        end
        if $nest_queue.empty?
        else
          if conversion_type =~ /AVX/
            $nest_queue.each{ |predicate|
              ret += Statement.new([tmp_predicate,Expression.new([:land,tmp_predicate,predicate,type]),type,nil]).convert_to_code(conversion_type) + "\n"
            }
          end
        end
        $predicate_queue.push(tmp_predicate)
        $nest_queue.push(tmp_predicate)
        $current_predicate = tmp_predicate
      when :elsif
        $nest_queue.pop
        tmp_predicate = $predicate_queue.pop
        $current_predicate = "pg#{$pg_count}"
        $varhash[$current_predicate] = [nil,type,nil,nil]
        $pg_count += 1
        nsplit = size/$min_element_size
        nsplit = 1 if conversion_type == "reference"
        pred_ops = Array.new
        accum_ops = Array.new
        for i in 0...nsplit
          exp = @expression
          suffix = String.new
          suffix = "_#{i}" if nsplit > 1
          ret += Declaration.new([type,$current_predicate+suffix]).convert_to_code(conversion_type)
          ret += Declaration.new([type,$accumulate_predicate+suffix]).convert_to_code(conversion_type)
          exp.get_related_variable.each{ |v|
            exp = exp.replace_recursive(v,v+suffix)
          }
          ret += Statement.new([$current_predicate+suffix,exp,type]).convert_to_code(conversion_type) + "\n"
          ret += Statement.new([$accumulate_predicate+suffix,$current_predicate+suffix,type,nil]).convert_to_code(conversion_type) + "\n"
          pred_ops.push($current_predicate + suffix)
          accum_ops.push($accumulate_predicate + suffix)
        end
        if nsplit > 1
          ret += Declaration.new([type,$current_predicate]).convert_to_code(conversion_type)
          ret += Declaration.new([type,$accumulate_predicate]).convert_to_code(conversion_type)
          ret += Statement.new([$current_predicate,Fusion.new([pred_ops,type,"B#{$min_element_size}"])]).convert_to_code(conversion_type) + "\n" if
          ret += Statement.new([$accumulate_predicate,Fusion.new([accum_ops,type,"B#{$min_element_size}"])]).convert_to_code(conversion_type) + "\n"
        end
        $nest_queue.each{ |predicate|
          ret += Statement.new([$current_predicate,Expression.new([:land,$current_predicate,predicate,type]),type,nil]).convert_to_code(conversion_type) + "\n"
        }
        $predicate_queue.push(tmp_predicate)
        $nest_queue.push($current_predicate)
      when :else
        $nest_queue.pop
        $current_predicate = $predicate_queue.pop
        ret += Statement.new([$current_predicate,Expression.new([:not,$accumulate_predicate,nil,type]),type,nil]).convert_to_code(conversion_type) + "\n"
        $nest_queue.each{ |predicate|
          ret += Statement.new([$current_predicate,Expression.new([:land,$current_predicate,predicate,type]),type,nil]).convert_to_code(conversion_type) + "\n"
        }
        $nest_queue.push($current_predicate)
        $predicate_queue.push($current_predicate)
      when :endif
        $predicate_queue.pop
        $nest_queue.pop
        abort "pg_count < 0" if $pg_count < 0
        $current_predicate = $predicate_queue.pop #"pg#{$pg_count}"
        $current_predicate = "svptrue_b#{$min_element_size}()" if $current_predicate == nil
        ret += "}\n"
      else
        abort "undefined operator of IfElseState: #{@operator}"
      end
    else
      abort "error: unsupported conversion_type #{conversion_type}"
    end
    ret
  end
end

class StoreState
  attr_accessor :dest,:src,:type
  def initialize(x)
    @dest, @src, @type = x
  end
  def get_type(h = $varhash)
    @type = @dest.get_type if @type == nil
    @type
  end
  def convert_to_code(conversion_type="reference")
    case conversion_type
    when "reference"
      ret = @dest.convert_to_code(conversion_type) + "=" + @src.convert_to_code(conversion_type) + ";"
    when /A64FX/
      ret = "svst1_#{get_type_suffix_a64fx(@type)}(#{$current_predicate},#{@dest.convert_to_code(conversion_type)},#{@src.convert_to_code(conversion_type)});"
    when /AVX2/
      if @type =~ /^(S|U)(64|32)/
        ret = "_mm256_storeu_si256((__m256i*)#{@dest.convert_to_code(conversion_type)},#{@src.convert_to_code(conversion_type)});"
      else
        ret = "_mm256_storeu_#{get_type_suffix_avx2(@type)}(#{@dest.convert_to_code(conversion_type)},#{@src.convert_to_code(conversion_type)});"
      end
    when /AVX-512/
      ret = "_mm512_store_#{get_type_suffix_avx512(@type)}(#{@dest.convert_to_code(conversion_type)},#{@src.convert_to_code(conversion_type)});"
    end
  end
end

class LoadState
  attr_accessor :dest,:src,:type
  def initialize(x)
    @dest, @src, @type = x
  end
  def get_type(h = $varhash)
    @type = @dest.get_type if @type == nil
    @type
  end
  def convert_to_code(conversion_type="reference")
    case conversion_type
    when "reference"
      ret = @dest.convert_to_code(conversion_type) + "=" + @src.convert_to_code(conversion_type) + ";"
    when /A64FX/
      ret = "#{@dest.convert_to_code(conversion_type)} = svld1_#{get_type_suffix_a64fx(type)}(#{$current_predicate},#{@src.convert_to_code(conversion_type)});"
    when /AVX2/
      if @type =~ /^(S|U)(64|32)/
        ret = "#{@dest.convert_to_code(conversion_type)} = _mm256_loadu_si256((__256i const*)#{@src.convert_to_code(conversion_type)});"
      else
        ret = "#{@dest.convert_to_code(conversion_type)} = _mm256_loadu_#{get_type_suffix_avx512(type)}(#{@src.convert_to_code(conversion_type)});"
      end
    when /AVX-512/
      ret = "#{@dest.convert_to_code(conversion_type)} = _mm512_load_#{get_type_suffix_avx512(type)}(#{@src.convert_to_code(conversion_type)});"
    end
  end
end

class PointerOf
  attr_accessor :type,:exp
  def initialize(x)
    @type,@exp = x
  end
  def get_pointer_type(type)
    case type
    when /F64/
      "double*"
    when /F32/
      "float*"
    when /S64/
      "long long int*"
      #"int64_t*"
    when /S32/
      "int*"
      #"int32_t*"
    when /U64/
      "unsigned long long int*"
      #"uint64_t*"
    when /U32/
      "unsigned int*"
      #"uint32_t*"
    end
  end

  def get_type(h = $varhash)
    @exp.get_type(h)
  end
  def convert_to_code(conversion_type="reference")
    case conversion_type
    when "reference"
      ret = @exp.convert_to_code(conversion_type)
    when /(A64FX|AVX2|AVX-512)/
      ret = "((#{get_pointer_type(@type)})&#{@exp.convert_to_code("reference")})"
    end
    ret
  end
end

class GatherLoad
  attr_accessor :dest,:src,:offset,:interval,:type
  def initialize(x)
    @dest,@src,@offset,@interval,@type = x
  end
  def convert_to_code(conversion_type,index = nil)
    $gather_load_count = 0 if $gather_load_count == nil
    ret = ""
    case conversion_type
    when /A64FX/
      index_name = "index_gather_load#{$gather_load_count}"
      vindex_name = "v" + index_name
      nelem = get_num_elem(@type,conversion_type)
      if index == nil
        index = "#{@offset.to_i}"
        for i in 1...nelem
          index +=  ",#{i*@interval.to_i + @offset.to_i}"
        end
      end
      size_single = get_single_data_size(@type)
      ret += "uint#{size_single}_t #{index_name}[#{nelem}] = {#{index}};\n"
      ret += "svuint#{size_single}_t #{vindex_name} = svld1_u#{size_single}(#{$current_predicate},#{index_name});\n"
      ret += "#{@dest.convert_to_code(conversion_type)} = svld1_gather_u#{size_single}index_#{get_type_suffix_a64fx(@type)}(#{$current_predicate},#{@src.convert_to_code(conversion_type)},#{vindex_name});"
    when /AVX2/
      index_name = "index_gather_load#{$gather_load_count}"
      vindex_name = "v" + index_name
      nelem = get_num_elem(@type,conversion_type)
      scale = get_byte_size(@type)
      if index == nil
        index = "#{@offset.to_i}"
        for i in 1...nelem
          index +=  ",#{i*@interval.to_i + @offset.to_i}"
        end
      end
      ret += "alignas(32) int #{index_name}[#{nelem}] = {#{index}};\n"
      index_simd_width = 32 * nelem
      ret += "__m#{index_simd_width}i #{vindex_name} = "
      case index_simd_width
      when 128
        ret += "_mm_load_si128((const __m128i*)#{index_name});\n"
      when 256
        ret += "_mm256_load_si256((const __m256i*)#{index_name});\n"
      end
      ret += "#{@dest.convert_to_code(conversion_type)} = _mm256_i32gather_#{get_type_suffix_avx2(@type)}(#{@src.convert_to_code(conversion_type)},#{vindex_name},#{scale});"
    when /AVX-512/
      index_name = "index_gather_load#{$gather_load_count}"
      vindex_name = "v" + index_name
      nelem = get_num_elem(type,conversion_type)
      size = get_single_data_size(@type)
      scale = get_byte_size(@type)
      if index == nil
        index = ""
        for i in 0...nelem
          index += "," if i > 0
          index +=  "#{i*@interval.to_i + @offset.to_i}"
        end
      end
      ret += "alignas(#{size}) int#{size}_t #{index_name}[#{nelem}] = {#{index}};\n"
      ret += "__m512i #{vindex_name} = _mm512_load_epi#{size}(#{index_name});\n"
      ret += "#{@dest.convert_to_code(conversion_type)} = _mm512_i#{size}gather_#{get_type_suffix_avx512(@type)}(#{vindex_name},#{@src.convert_to_code(conversion_type)},#{scale});"
    else
      abort "unsupported conversion type for GatherLoad"
    end
    $gather_load_count += 1
    ret
  end
end

class ScatterStore
  attr_accessor :dest,:src,:offset,:interval,:type
  def initialize(x)
    @dest,@src,@offset,@interval,@type = x
  end
  def convert_to_code(conversion_type,index = nil)
    $scatter_store_count = 0 if $scatter_store_count == nil
    ret = ""
    case conversion_type
    when /A64FX/
      index_name = "index_scatter_store#{$scatter_store_count}"
      vindex_name = "v" + index_name
      nelem = get_num_elem(@type,conversion_type)
      if index == nil
        index = ""
        for i in 0...nelem
          index += "," if i > 0
          index += "#{@offset.to_i + i*@interval.to_i}"
        end
      end
      size_single = get_single_data_size(@type)
      ret += "uint#{size_single}_t #{index_name}[#{nelem}] = {#{index}};\n"
      ret += "svuint#{size_single}_t #{vindex_name} = svld1_u#{size_single}(#{$current_predicate},#{index_name});\n"
      ret += "svst1_scatter_u#{size_single}index_#{get_type_suffix_a64fx(@type)}(#{$current_predicate},#{@dest.convert_to_code(conversion_type)},#{vindex_name},#{@src.convert_to_code(conversion_type)});"
    when /AVX2/
      nelem = get_num_elem(@type,conversion_type)
      ret += "{\n"
      ret += NonSimdDecl.new([@type,"__fkg_store_tmp[#{nelem}]"]).convert_to_code(conversion_type)
      ret += StoreState.new(["__fkg_store_tmp",@src,@type]).convert_to_code(conversion_type) + "\n"

      if index == nil
        index = Array.new
        for i in 0...nelem
          index.push((@offset.to_i+@interval.to_i*i).to_s)
        end
      else
        index = index.split(",")
      end

      index.each_with_index{ |a,b|
        dest = Expression.new([:array,@dest,a,@type]).convert_to_code(conversion_type)
        ret += NonSimdState.new([dest,"__fkg_store_tmp[#{b}]",@type]).convert_to_code(conversion_type) + "\n"
      }
      ret += "}\n"
    when /AVX-512/
      index_name = "index_scatter_store#{$scatter_store_count}"
      vindex_name = "v" + index_name
      size = get_single_data_size(@type)
      nelem = get_num_elem(@type,conversion_type)
      case @type
      when /64/
        scale = 8
      when /32/
        scale = 4
      when /16/
        scale = 2
      end
      if index == nil
        index = ""
        for i in 0...nelem
          index += "," if i > 0
          index +=  "#{@offset.to_i + i*@interval.to_i}"
        end
      end
      ret += "int#{size}_t #{index_name}[#{nelem}] = {#{index}};\n"
      ret += "__m512i #{vindex_name} = _mm512_load_epi#{size}(#{index_name});\n"
      ret += "_mm512_i#{size}scatter_#{get_type_suffix_avx512(@type)}(#{@dest.convert_to_code(conversion_type)},#{vindex_name},#{@src.convert_to_code(conversion_type)},#{scale});\n"
    else
      abort "unsupported conversion type for ScatterStore"
    end
    $scatter_store_count += 1
    ret
  end
end

class StructureLoad
  attr_accessor :dest,:src,:nelem,:type
  def initialize(x)
    @dest,@src,@nelem,@type = x
    abort "nelem must be 1 to 4" if nelem > 4 || nelem <= 0
  end
  def convert_to_code(conversion_type="reference")
    ret = ""
    case conversion_type
    when /A64FX/
      ret = "#{@dest.convert_to_code(conversion_type)} = svld#{nelem}_#{get_type_suffix_a64fx(@type)}(#{$current_predicate},#{@src.convert_to_code(conversion_type)});"
    when /AVX-512/
      ret = ""
      type_single = get_single_element_type(@type)
      for i in 0...@nelem
        dest = Expression.new([:dot,@dest,["x","y","z","w"][i],type_single])
        ret += GatherLoad.new([dest,@src,i,@nelem,type_single]).convert_to_code(conversion_type) + "\n"
      end
    when /AVX2/
      #abort "nelem of StructureLoad must == 2 or == 4"
      if true #@nelem != 2 && @nelem != 4
        type_single = get_single_element_type(@type)
        for i in 0...@nelem
          dest = Expression.new([:dot,@dest,["x","y","z","w"][i],type_single])
          ret += GatherLoad.new([dest,@src,i,@nelem,type_single]).convert_to_code(conversion_type) + "\n"
        end
      else
        type_single = get_single_element_type(@type)
        suffix = get_type_suffix_avx(type_single)
        decl = get_declare_type_avx2(type_single)
        #src  = @src.convert_to_code(conversion_type)
        dest = @dest.convert_to_code(conversion_type)
        width = get_num_elem(@type,conversion_type)

        for i in 0...@nelem
          tmp = Expression.new([:array,@src.exp.lop,NonSimdExp.new([:plus,@src.exp.rop,(width*i/@nelem).to_s,@src.exp.type]),@src.exp.type])
          src = PointerOf.new([@src.type,tmp])
          tmp = LoadState.new([Expression.new([:dot,@dest,["x","y","z","w"][i],type_single]),src,type_single])
          ret += tmp.convert_to_code(conversion_type) + "\n"
        end
        ret += "unpack#{@nelem}x#{@nelem}_#{suffix}(#{dest});\n"
      end
    else
      abort "error: unsupported conversion_type #{conversion_type} for StructureLoad"
    end
    ret
  end
end

class StructureStore
  attr_accessor :dest,:src,:nelem,:type
  def initialize(x)
    @dest,@src,@nelem,@type = x
    abort "nelem must be 1 to 4" if nelem > 4 || nelem <= 0
  end
  def convert_to_code(conversion_type="reference")
    ret = ""
    case conversion_type
    when /A64FX/
      #ret = "svst#{nelem}_#{get_type_suffix_a64fx(@type)}(#{$current_predicate},#{@dest.convert_to_code(conversion_type)},#{@src.convert_to_code(conversion_type)});"
      for i in 0...@nelem
        src = Expression.new([:dot,@src,["x","y","z","w"][i],@type])
        ret += ScatterStore.new([@dest,src,i,@nelem,@type]).convert_to_code(conversion_type)
      end
    when /AVX-512/
      for i in 0...@nelem
        src = Expression.new([:dot,@src,["x","y","z","w"][i],@type])
        ret += ScatterStore.new([@dest,src,i,@nelem,@type]).convert_to_code(conversion_type)
      end
    when /AVX2/
      #abort "nelem of StructureStore must == 2 or == 4" if @nelem != 2 && @nelem != 4
      if true #@nelem != 2 && @nelem != 4
        for i in 0...@nelem
          src = Expression.new([:dot,@src,["x","y","z","w"][i],@type])
          abort if @type == nil
          ret += ScatterStore.new([@dest,src,i,@nelem,@type]).convert_to_code(conversion_type)
        end
      else
        suffix = get_type_suffix_avx2(@type)
        decl = get_declare_type_avx2(@type)
        src  = @src.convert_to_code(conversion_type)
        width = get_num_elem(@type,conversion_type)

        ret += "pack#{@nelem}x#{@nelem}_#{suffix}(#{src});\n"
        for i in 0...@nelem
          tmp = Expression.new([:array,@dest.exp.lop,NonSimdExp.new([:plus,@dest.exp.rop,(width*i/@nelem).to_s,@dest.exp.type]),@dest.exp.type])
          dest = PointerOf.new([@dest.type,tmp])
          ret += StoreState.new([dest,Expression.new([:dot,@src,["x","y","z","w"][i],@type]),@type]).convert_to_code(conversion_type) + "\n"
        end
      end
    else
      abort "error: unsupporeted conversion_type #{conversion_type} for StructureStore"
    end
    ret
  end
end

class Duplicate
  attr_accessor :name,:expression,:type
  def initialize(x)
    @name,@expression,@type = x
  end
  def get_type(h = $varhash)
    @type = @expression.get_type(h) if @type == nil
    @type
  end
  def convert_to_code(conversion_type="reference")
    case conversion_type
    when "reference"
      ret = "#{@name.convert_to_code(conversion_type)} = #{@expression.convert_to_code(conversion_type)};"
    when /A64FX/
      ret = "#{@name.convert_to_code(conversion_type)} = svdup_n_#{get_type_suffix_a64fx(@type)}(#{@expression.convert_to_code(conversion_type)});"
    when /AVX2/
      set1_suffix = ""
      set1_suffix = "x" if @type =~ /(S|U)64/
      ret = "#{@name.convert_to_code(conversion_type)} = _mm256_set1_#{get_type_suffix_avx2(@type)}#{set1_suffix}(#{@expression.convert_to_code(conversion_type)});"
    when /AVX-512/
      ret = "#{@name.convert_to_code(conversion_type)} = _mm512_set1_#{get_type_suffix_avx512(@type)}(#{@expression.convert_to_code(conversion_type)});"
    end
    ret
  end
end

class MADD
  attr_accessor :operator, :aop, :bop, :cop, :type
  def initialize(x)
    @operator, @aop, @bop, @cop, @type = x
  end

  def derive_type(op,aop,bop,cop,h=$varhash)
    type=nil
    at = aop.get_type(h)
    bt = bop.get_type(h)
    ct = cop.get_type(h)
    if at == "F64" && bt == "F64" && ct == "F64"
      type="F64"
    elsif at == "F32" && bt == "F32" && ct == "F32"
      type="F32"
    elsif at == "F16" && bt == "F16" && ct == "F16"
      type="F16"
    else
      abort "unsupported operand types for MADD operation"
    end
    type
  end
  
  def get_type(h = $varhash)
    if @type == nil
      @type = derive_type(@operator,@aop,@bop,@cop,h)
    end
    @type
  end

  def get_related_variable
    ret = []
    ret += @aop.get_related_variable
    ret += @bop.get_related_variable
    ret += @cop.get_related_variable
    ret.sort!
    ret.uniq!
    ret
  end
  def isJRelated(list)
    ret = false
    ret |= @aop.isJRelated(list)
    ret |= @bop.isJRelated(list)
    ret |= @cop.isJRelated(list)
    ret
  end

  def replace_recursive(orig,replaced)
    aop = @aop.replace_recursive(orig,replaced)
    bop = @bop.replace_recursive(orig,replaced)
    cop = @cop.replace_recursive(orig,replaced)
    MADD.new([@operator,aop,bop,cop,@type])
  end

  def replace_fdpsname_recursive(h=$varhash)
    aop = @aop.replace_fdpsname_recursive(h)
    bop = @bop.replace_fdpsname_recursive(h)
    cop = @cop.replace_fdpsname_recursive(h)
    MADD.new([@operator,aop,bop,cop,@type])
  end

  def replace_by_list(name_list,replaced_list)
    name_list.zip(replaced_list){ |n,r|
      self.replace_recursive(n,r)
    }
  end
  def replace_by_hash(hash)
    hash.each{ |k,v|
      self.replace_recursive(k,v)
    }
  end
  
  def convert_to_code(conversion_type="reference")
    retval=""
    case conversion_type
    when "reference" then
      case @operator
      when :madd then
      #retval = "std::fma(#{@aop.convert_to_code(conversion_type)},#{@bop.convert_to_code(conversion_type)},#{@cop.convert_to_code(conversion_type)})"
        retval = "(#{@aop.convert_to_code(conversion_type)}*#{@bop.convert_to_code(conversion_type)}+#{@cop.convert_to_code(conversion_type)})"
      when :msub then
        retval = "(#{@cop.convert_to_code(conversion_type)} - #{@aop.convert_to_code(conversion_type)}*#{@bop.convert_to_code(conversion_type)})"
      when :nmadd then
        retval = "(-(#{@cop.convert_to_code(conversion_type)} + #{@aop.convert_to_code(conversion_type)}*#{@bop.convert_to_code(conversion_type)}))"
      when :nmsub then
        retval = "(#{@aop.convert_to_code(conversion_type)}*#{@bop.convert_to_code(conversion_type)}-#{@cop.convert_to_code(conversion_type)})"
      else
        abort "error: unsupported operator for MADD class"
      end
    when /A64FX/
      self.get_type
      case @operator
      when :madd
        retval += "svmad_#{get_type_suffix_a64fx(@type)}_z("
      when :msub
        retval += "svmsb_#{get_type_suffix_a64fx(@type)}_z("
      when :nmadd
        retval += "svnmad_#{get_type_suffix_a64fx(@type)}_z("
      when :nmsub
        retval += "svnmsb_#{get_type_suffix_a64fx(@type)}_z("
      end
      $current_predicate = "svptrue_b#{$min_element_size}()" if $current_predicate == nil
      retval += $current_predicate + ","
      retval += @aop.convert_to_code(conversion_type) + ","
      retval += @bop.convert_to_code(conversion_type) + ","
      retval += @cop.convert_to_code(conversion_type) + ")"
    when /AVX2/
      self.get_type
      suffix = get_type_suffix_avx2(@type)
      if @type =~ /F/
        case @operator
        when :madd
          retval += "_mm256_fmadd_#{suffix}("
        when :msub
          retval += "_mm256_fnmadd_#{suffix}("
        when :nmadd
          retval += "_mm256_fnmsub_#{suffix}("
        when :nmsub
          retval += "_mm256_fmsub_#{suffix}("
        end
        retval += @aop.convert_to_code(conversion_type) + ","
        retval += @bop.convert_to_code(conversion_type) + ","
        retval += @cop.convert_to_code(conversion_type) + ")"
      else
        case @operator
        when :madd
          retval += "_mm256_add_#{suffix}(_mm256_mul_#{suffix}(" + @aop.convert_to_code(conversion_type) + "," + @bop.convert_to_code(conversion_type) + ")," + @cop.convert_to_code(conversion_type) + ")"
        when :msub
          retval += "_mm256_sub_#{suffix}(" + @cop.convert_to_code(conversion_type) + ",_mm256_mul_#{suffix}(" + @aop.convert_to_code(conversion_type) + "," + @bop.convert_to_code(conversion_type) + "))"
        when :nmadd
          retval += "_mm256_sub_#{suffix}(_ mm256_set1_#{suffix}(0)," + MADD.new([:madd,@aop,@bop,@cop,@type]).convert_to_code(conversion_type) + ")"
        when :nmsub
          retval += "_mm256_sub_#{suffix}(_mm256_mul_#{suffix}(" + @aop.convert_to_code(conversion_type) + "," + @bop.convert_to_code(conversion_type) + ")," + @cop.convert_to_code(conversion_type) + ")"
        end
      end
    when /AVX-512/
      self.get_type
      case @operator
      when :madd
        retval += "_mm512_fmadd_#{get_type_suffix_avx512(@type)}("
      when :msub
        retval += "_mm512_fnmadd_#{get_type_suffix_avx512(@type)}("
      when :nmadd
        retval += "_mm512_fnmsub_#{get_type_suffix_avx512(@type)}("
      when :nmsub
        retval += "_mm512_fmsub_#{get_type_suffix_avx512(@type)}("
      end
      retval += @aop.convert_to_code(conversion_type) + ","
      retval += @bop.convert_to_code(conversion_type) + ","
      retval += @cop.convert_to_code(conversion_type) + ")"
    end
    retval
  end
end

class Merge
  attr_accessor :op1,:op2,:type, :split_index
  def initialize(x)
    @op1,@op2,@type,@split_index = x
  end
  def get_related_variable
    [@op1,@op2]
  end
  def replace_recursive(orig,replaced)
    @op1 = @op1.replace_recursive(orig,replaced)
    @op2 = @op2.replace_recursive(orig,replaced)
    self.dup
  end
  def convert_to_code(conversion_type="reference")
    ret = ""
    case conversion_type
    when /reference/
      ret = "#{@op1.convert_to_code(conversion_type)}"
    when /A64FX/
      predicate = $current_predicate
      predicate += "_#{@split_index}" if @split_index != nil
      ret = "svsel_#{get_type_suffix_a64fx(@type)}(#{predicate},#{@op1.convert_to_code(conversion_type)},#{@op2.convert_to_code(conversion_type)});"
    when /AVX2/
      suffix = get_type_suffix_avx2(@type)
      predicate = $current_predicate
      predicate += "_#{@split_index}" if @split_index != nil
      # inactive elements come from first input
      if suffix == "epi32"
        ret = "_mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(#{@op2.convert_to_code(conversion_type)}),_mm256_castsi256_ps(#{@op1.convert_to_code(conversion_type)}),#{predicate}));"
      elsif suffix == "epi64"
        #predicate = "_mm256_castsi256_pd(" + predicate + ")"
        predicate = "_mm256_castps_pd(" + predicate + ")"
        ret = "_mm256_blendv_pd(_mm256_castsi256_pd(#{@op2.convert_to_code(conversion_type)}),_mm256_castsi256_pd(#{@op1.convert_to_code(conversion_type)}),#{predicate})"
        ret = "_mm256_castpd_si256(#{ret});"
      else
        predicate = "_mm256_castps_pd(" + predicate + ")" if suffix == "pd"
        ret = "_mm256_blendv_#{suffix}(#{@op2.convert_to_code(conversion_type)},#{@op1.convert_to_code(conversion_type)},#{predicate});"
      end
    when /AVX-512/
      predicate = $current_predicate
      predicate += "_#{@split_index}" if @split_index != nil
     ret = "_mm512_mask_blend_#{get_type_suffix_avx512(@type)}(#{predicate},#{@op2.convert_to_code(conversion_type)},#{@op1.convert_to_code(conversion_type)});" # inactive elements come from first input
    end
    ret
  end
end

class Expression
  attr_accessor :operator, :lop, :rop, :type
  def initialize(x)
    @operator, @lop, @rop, @type=x
  end

  def derive_type (operator,lop,rop,h=$varhash)
    type=nil
    lt = lop.get_type(h)
    rt = rop.get_type(h) if rop != nil && !["x","y","z","w"].index(rop) && operator != :array
    abort "single element size of #{lop.convert_to_code("reference")}(#{lt}) and #{rop.convert_to_code("reference")}(#{rt}) are different" if rt != nil && get_single_data_size(lt) != get_single_data_size(rt)
    #print "derive type ", operator," ", lop," ", rop, "\n"
    #p self
    if [:plus, :minus, :mult, :div,:and,:or].index(operator)
      rt = rop.get_type(h)
      abort "type is not derived for #{lop} = #{lt}, #{rop} = #{rt} in derive_type" if lt == nil || rt == nil
      #p lop,rop
      #print "#{lt},#{rt}\n"
      #p h
      if operator == :mult && lt.index("vec")  && rt.index("vec")
        if lt == "F64vec" || rt == "F64vec"
          type="F64"
        elsif lt == "F32vec" || rt == "F32vec"
          type="F32"
        else
          type="F16"
        end
      else
        type="U"
        type="S" if lt.index("S") || rt.index("S")
        type="F" if lt.index("F") || rt.index("F")
        if lt.index("64") ||  rt.index("64")  
          type += "64"
        elsif lt.index("32") || rt.index("32")
          type += "32"
        else
          p lop,rop,h
          nil[0]
          abort "U16 is currently not supported"
          type += "16"
        end
        type += "vec" if lt.index("vec") ||  rt.index("vec")  
      end
    elsif [:land,:lor,:eq,:neq,:lt,:le,:gt,:ge,:landnot].index(operator)
      rt = rop.get_type(h)
      abort "type is not derived for #{lop} = #{lt}, #{rop} = #{rt} in derive_type" if lt == nil || rt == nil
      abort "vector type comparison is not supported. (#{lt} #{rt})" if lt =~ /(vec|mat)/ || rt =~ /(vec|mat)/
      if lt =~ /64/ || rt =~ /64/
        type = "B64"
      elsif lt =~ /32/ || rt =~ /32/
        type = "B32"
      else
        type = "B16"
      end
    elsif operator == :uminus || operator == :array || operator == :not
      type = lt
    elsif operator == :dot
      if ["x","y","z"].index(rop)
        if lt == "F64vec"
          type = "F64"
        elsif lt == "F32vec"
          type = "F32"
        elsif lt == "F16vec"
          type = "F16"
        else
          nil[0]
          warn "error: #{lop} is not vector type!"
          abort
        end
      elsif lop =~ /^\d+$/ && rop =~ /^\d+(f|h)?$/
        if rop =~ /f/ then
          type = "F32"
        elsif rop =~ /h/ then
          type = "F16"
        else
          type = "F64"
        end
      else
        warn rop
        abort "error: unsupported vector member"
      end
    else
      p self
      abort "error: undefined reference to #{self} in derive_type"
      type=rt
    end
    #print "final type=", type, "\n"
    type
  end

  def  get_type(h = $varhash)
    if @operator == :dot && @lop =~ /^\d+$/
      if @rop =~ /^\d+h$/
        @type = "F16"
      elsif @rop =~ /^\d+f$/
        @type = "F32"
      elsif @rop =~ /^\d+$/
        @type = "F64"
      else
        p self
        abort "unsupported floatinng point number! (must be (\d+).(\d+[f,h]))"
      end
    end
    if @type == nil
      @type = derive_type(@operator, @lop, @rop, h)
    end
    abort "get_type failed" if @type == nil
    @type
  end

  def fusion_iotag(iotag)
    if @lop == iotag && @operator == :dot
      ret = @lop + "_" + @rop
    else
      @lop = @lop.fusion_iotag(iotag)
      @rop = @rop.fusion_iotag(iotag) if @rop != nil
      ret = self
    end
    ret
  end

  def find_function
    ret = []
    if @lop.class == FuncCall
      ret += [@lop]
    elsif !isLeaf(@lop)
      ret += @lop.find_function
    end
    if @rop != nil
      if @rop.class == FuncCall
        ret += [@rop]
      elsif !isLeaf(@rop)
        ret += @rop.find_function
      end
    end
    ret
  end

  def get_related_variable
    ret = []
    case @operator
    when :dot
      ret += [@lop] if @lop =~ /[a-z_].*/
    when :array
      ret += [@lop]
    else
      ret += @lop.get_related_variable
      ret += @rop.get_related_variable if @rop != nil
    end
    ret.sort!
    ret.uniq!
    ret
  end
  def isJRelated(list)
    ret = false
    if isLeaf(self)
      abort "only . and [] operator can be Leaf for Expression class" if @operator != :dot && @operator != :array
      if @lop =~ /\d+/
      else
        ret |= @lop.isJRelated(list)
      end
    else
      ret |= @lop.isJRelated(list)
      ret |= @rop.isJRelated(list) if @rop != nil
    end
    ret
  end

  def get_opcode(op)
    case op
    when :plus
      "+"
    when :minus
      "-"
    when :mult
      "*"
    when :div
      "/"
    when :lt
      "<"
    when :le
      "<="
    when :gt
      ">"
    when :ge
      ">="
    when :eq
      "=="
    when :neq
      "!="
    when :eq
      "=="
    when :and
      "&"
    when :or
      "|"
    when :land
      "&&"
    when :lor
      "||"
    else
      abort "error: undefined opcode #{op}"
    end
  end

  def replace_recursive(orig,replaced)
    lop = @lop.dup
    rop = @rop.dup if @rop != nil
    type = @type.dup
    if @operator == :dot
      if orig.class == Expression && orig.operator == :dot
        if lop.class == String
          lop = replaced if lop == orig.lop && rop == orig.rop
        end
      else
        lop = lop.replace_recursive(orig,replaced)
        rop = rop.replace_recursive(orig,replaced) if rop != nil
      end
    else
      lop = lop.replace_recursive(orig,replaced)
      rop = rop.replace_recursive(orig,replaced) if rop != nil
    end
    Expression.new([@operator,lop,rop,type])
  end
  def replace_fdpsname_recursive(h=$varhash)
    #p self
    ret = self.dup
    if @operator == :array
      name = ret.lop.dup
      iotype = h[name][0]
      type   = h[name][1]
      fdpsname = h[name][2]
      abort "array expression must be used for EPI, EPJ, or FORCE variable" if fdpsname == nil || iotype == nil
      ret = Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),ret.rop]),fdpsname,type])
    elsif @operator == :dot
      ret.lop = ret.lop.replace_fdpsname_recursive(h)
    else
      ret.lop = ret.lop.replace_fdpsname_recursive(h)
      ret.rop = ret.rop.replace_fdpsname_recursive(h) if ret.rop != nil
    end
    ret
  end

  def replace_by_list(name_list,replaced_list)
    name_list.zip(replaced_list){ |n,r|
      self.replace_recursive(n,r)
    }
  end
  def replace_by_hash(hash)
    hash.each{ |k,v|
      self.replace_recursive(k,v)
    }
  end

  def convert_to_code(conversion_type="reference")
    #p self
    retval=""
    case conversion_type
    when "reference"
      if type=[:plus,:minus,:mult,:div,:eq,:neq,:lt,:le,:gt,:ge,:and,:or,:land,:lor].index(@operator)
        opcode=get_opcode(@operator)
        retval="("+@lop.convert_to_code(conversion_type)+opcode+
               @rop.convert_to_code(conversion_type)+")"
      elsif @operator == :uminus
        retval="-(" + @lop.convert_to_code(conversion_type) + ")"
      elsif @operator == :dot
        retval=@lop.convert_to_code(conversion_type)+"."+
               @rop.convert_to_code(conversion_type)
      elsif @operator == :array
        retval=lop.convert_to_code(conversion_type) + "[" + @rop.convert_to_code("reference") + "]"
      elsif @operator == :func
        retval= @lop+"("+ @rop.convert_to_code(conversion_type)+")"
      else
        abort "error: unsupported operator #{@operator} in expression code conversion"
      end
    when /A64FX/
      retval = self.convert_to_code_a64fx(conversion_type)
    when /AVX2/
      retval = self.convert_to_code_avx2(conversion_type)
    when /AVX-512/
      retval = self.convert_to_code_avx512(conversion_type)
    end
    retval
  end
end

class NonSimdExp < Expression
  def convert_to_code(conversion_type="reference")
    super("reference")
  end
end

class Fusion
  attr_accessor :ops,:from, :to
  def initialize(x)
    @ops, @from, @to = x
    abort "type of Fusion must be specified" if @from == nil || @to == nil
  end
  def convert_to_code(conversion_type)
    ret = String.new
    n = @ops.length
    case conversion_type
    when "reference"
    # do nothing
    when "AVX2"
      case n
      when 1
        abort "Fusion does not happen for single operand"
      when 2
        if @from == "F64" && @to == "F32"
          ret = "_mm256_permute2f128_ps("
          for i in 0...n
            ret += "," if i>0
            ret += "_mm256_castps128_ps256(_mm256_cvtpd_ps(#{ops[i]}))"
          end
          ret += ",0x20)"
        elsif @from == "B64" && @to == "B32"
          ret = "_mm256_permute2f128_ps("
          for i in 0...n
            ret += "," if i>0
            ret += "_mm256_castsi256_ps(_mm256_castsi128_si256(_mm256_cvtepi64_epi32(_mm256_castps_si256(#{ops[i]}))))"
          end
          ret += ",0x20)"
        else
          abort "Fusion from #{@from} to #{@to} with AVX2 under construction"
        end
      else
        abort "Fusion of #{n} does not happen for AVX2"
      end

    when "AVX-512"
      case n
      when 1
        abort "Fusion does not happen for single operand"
      when 2
        if @from == "F64" && @to == "F32"
          ret = "_mm512_insertf32x8(_mm512_castps256_ps512(_mm512_cvtpd_ps(#{ops[0]})),_mm512_cvtpd_ps(#{ops[1]}),1)"
        end
      else
        abort "Fusion #{n} for AVX2 under construction"
      end
    else
      abort "unsupported conversion_type #{coversion_type} for Fusion"
    end
    ret
  end
end

class Fission
  attr_accessor :op,:from, :to, :index
  def initialize(x)
    @op, @from, @to, @index = x
    abort "type of Fision must be specified" if @from == nil || @to == nil
    abort "index must be specified" if @index == nil
  end
  def convert_to_code(conversion_type)
    ret = String.new
    case conversion_type
    when "reference"
    # do nothing
    when "AVX2"
      if @from == "F32" && @to == "F64"
        ret = @op.convert_to_code(conversion_type)
        ret = "_mm256_castps256_ps128(#{ret})"  if @index == 0
        ret = "_mm256_extractf128_ps(#{ret},1)" if @index == 1
        ret = "_mm256_cvtps_pd(#{ret})"
      else
        abort "under costruction for AVX2 : Fission"
      end
    when "AVX-512"
      if @from == "F32" && @to == "F64"
        ret = @op.convert_to_code(conversion_type)
        ret = "_mm512_castps512_ps256(#{ret})"  if @index == 0
        ret = "_mm512_extractf32x8_ps(#{ret},1)" if @index == 1
        ret = "_mm512_cvtps_pd(#{ret})"
      else
        abort "under costruction for AVX2 : Fission"
      end
    else
      abort "unsupported conversion_type #{conversion_type} for Fission"
    end
    ret
  end
end
