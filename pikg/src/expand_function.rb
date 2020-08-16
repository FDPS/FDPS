class Kernelprogram
  def expand_function
    new_s = []
    @statements.each{|s|
      if isStatement(s)
        new_s += s.expand_function
      else
        new_s += [s]
      end
    }
    @statements = new_s
  end
end

class Statement
  def expand_function
    exp,statements = @expression.expand_function
    tmp = Statement.new([@name,exp])
    statements += [tmp]
    #p exp
    # statements.each{|s| p s}
    statements
  end
  def replace_variable(h)
    Statement.new([@name.dup.replace_variable(h),@expression.replace_variable(h)])
  end
end

class Expression
  def expand_function
    ret = self.dup
    statements = []
    if !isLeaf(self)
      ret.lop,lstat = @lop.expand_function
      statements += lstat
      if @rop != nil
        ret.rop,rstat = @rop.expand_function
        statements += rstat
      end
    end
    [ret,statements]
  end

  def replace_variable(h)
    ret = Expression.new([@operator,@lop.dup,@rop.dup])
    ret.lop = ret.lop.replace_variable(h) if ret.lop != nil
    ret.rop = ret.rop.replace_variable(h) if ret.rop != nil
    ret
  end
end

class FuncCall
  def expand_function
    ret = self.dup
    statements = []
    if !$reserved_function.index(@name)
      abort "undefined reference to function #{@name}" if $funchash[@name] == nil
      function = $funchash[@name]
      #abort "piecewise polynomial approximation is tested" if function.decl.attr == "__ppa__"
      h = []
      @ops.zip(function.decl.vars).each{ |a,b|
        tmp = add_new_tmpvar(a.get_type) #"__fkg_tmp#{$tmp_val_count}"
        h += [b,tmp]
        stmp = Statement.new([tmp,a,a.get_type])
        statements += [stmp]
      }
      function.statements.each{ |s|
        tmp = add_new_tmpvar(s.get_type) #"__fkg_tmp#{$tmp_val_count}"
        h += [s.name,tmp]
      }
      h = Hash[*h]
      #retval = Statement.new([tmp,function.retval.ret.replace_variable(h)])
      ret = function.retval.ret.replace_variable(h)
      s_tmp = []
      function.statements.each{ |s|
        replaced = s.replace_variable(h)
        s_tmp += replaced.expand_function
      }
      #s_tmp += [retval]
      #s_tmp.each{ |s| p s}

      statements += s_tmp
    end
    [ret,statements]
  end

  def replace_variable(h)
    if h[@name] != nil
      ret = h[@name]
    else
      ops = []
      @ops.each{ |x|
        ops += [x.replace_variable(h)]
      }
      ret = FuncCall.new([@name,ops])
      ret.type = @type
    end
    ret
  end
end

class IfElseState
  def expand_function
    [self]
  end
end

class FloatingPoint
  def expand_function
    [self,[]]
  end
  def replace_variable(h)
    self
  end
end
class IntegerValue
  def expand_function
    [self,[]]
  end
  def replace_variable(h)
    self
  end
end

class String
  def expand_function
    [self,[]]
  end

  def replace_variable(h)
    if h[self] != nil
      #print "#{self} -> #{h[self]}\n"
      h[self]
    else
      self.dup
    end
  end
end

