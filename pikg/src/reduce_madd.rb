require_relative "common.rb"

class Statement
  def reduce_madd
    @expression = @expression.reduce_madd
    self
  end
  def reduce_madd_recursive
    @expression = @expression.reduce_madd_recursive if !isLeaf(@expression)
    self
  end

  def reduce_negate_recursive
    @expression = @expression.reduce_negate_recursive if !isLeaf(@expression)
  end
end

class Expression
  def find_mul(replaced)
    ret = self
    #p "finding mul:",self,@operator
    case @operator
    when :mult
      found = true
      leaf = true
    when :plus
      ret,found,leaf = @rop.find_mul(replaced) if !isLeaf(@rop)
      @rop = replaced if found && leaf
      if !found
        ret,found,leaf = @lop.find_mul(replaced) if !isLeaf(@lop)
        @lop = replaced if found && leaf
      end
      leaf = false if found
    when :minus
      ret,found,leaf = @rop.find_mul(replaced) if !isLeaf(@rop)
      if found && leaf
        @rop = replaced
        @operator = :plus
      end
      ret = Expression.new([:uminus,ret,nil]) if found
      if !found
        ret,found,leaf = @lop.find_mul(replaced) if !isLeaf(@lop)
        @lop = replaced if found && leaf
      end
      leaf = false if found
    when :uminus
      ret,found,leaf = @lop.find_mul(replaced) if !isLeaf(@lop)
      ret = Expression.new([:uminus,ret,nil]) if found
      leaf = false if found
    else
      found = false
      leaf = false
    end
    [ret,found,leaf]
  end

  def reduce_madd_recursive2
    tmp = self
    if @lop.operator == :uminus
      case @operator
      when :plus
        tmp = Expression.new([:minus,@rop,@lop.lop,@type])
      when :minus
        tmp = Expression.new([:uminus,Expression.new([:plus,@lop.lop,@rop]),nil,@type])
      when :mult
      when :div
      end
    end

    case @opreator
    when :plus
      if @lop.operator == :mult
        ret = MADD.new([:madd,@lop.lop,@lop.rop,@rop,@type])
      end
    when :minus
    else
    end
  end

  def reduce_madd_recursive
    ret = self
    if @operator == :plus || @operator == :minus
      if !isLeaf(@rop)
        #p "reducing mad right:",self
        mul,found,leaf = @rop.find_mul(@lop.dup) if @operator == :plus
        mul,found,leaf = @rop.find_mul(Expression.new([:uminus,@lop,nil])) if @operator == :minus
        @rop = @lop if leaf && @operator == :plus
        @rop = Expression.new([:uminus,@lop,nil]) if leaf && @operator == :minus
        if found
          #p "found mul:",mul
          case @operator
          when :plus
            case mul.operator
            when :mult
              ret = MADD.new([:madd,mul.lop,mul.rop,@rop])
            when :uminus
              abort "mult is not found actually" if mul.lop.operator != :mult
              ret = MADD.new([:msub,mul.lop.lop,mul.lop.rop,@rop])
            end
          when :minus
            case mul.operator
            when :mult
              ret = MADD.new([:nmadd,mul.lop,mul.rop,rop])
            when :uminus
              abort "mult is not found actually" if mul.lop.operator != :mult
              ret = MADD.new([:nmsub,mul.lop.lop,mul.lop.rop,@rop])
            end
          end
          ret.get_type
          #p "temporaly madd:",ret
          ret.aop = ret.aop.reduce_madd_recursive if !isLeaf(ret.aop)
          ret.bop = ret.bop.reduce_madd_recursive if !isLeaf(ret.bop)
          ret.cop = ret.cop.reduce_madd_recursive if !isLeaf(ret.cop)
        end
      end
      if !isLeaf(@lop) && !found
        #p "reducing mad left:",self
        mul,found,leaf = @lop.find_mul(@rop)
        @lop = @rop if leaf
        if found
          #p "found mul:",mul
          case @operator
          when :plus
            case mul.operator
            when :mult
              ret = MADD.new([:madd,  mul.lop,mul.rop,@lop,@type])
            when :uminus
              abort "mult is not found actually" if mul.lop.operator != :mult
              ret = MADD.new([:msub, mul.lop.lop,mul.lop.rop,@lop,@type])
            end
          when :minus
            case mul.operator
            when :mult
              ret = MADD.new([:nmsub,mul.lop,mul.rop,@lop,@type])
            when :uminus
              abort "mult is not found actually" if mul.lop.operator != :mult
              ret = MADD.new([:nmadd, mul.lop.lop,mul.lop.rop,@lop,@type])
            end
          end
         #p "temporaly madd:",ret
          ret.aop = ret.aop.reduce_madd_recursive if !isLeaf(ret.aop)
          ret.bop = ret.bop.reduce_madd_recursive if !isLeaf(ret.bop)
          ret.cop = ret.cop.reduce_madd_recursive if !isLeaf(ret.cop)
        end
      end
      if !found
        ret.lop = ret.lop.reduce_madd_recursive if !isLeaf(ret.lop)
        ret.rop = ret.rop.reduce_madd_recursive if ret.rop != nil && !isLeaf(ret.rop)
      end
    else
      ret.lop = ret.lop.reduce_madd_recursive if !isLeaf(ret.lop)
      ret.rop = ret.rop.reduce_madd_recursive if ret.rop != nil && !isLeaf(ret.rop)
    end
    ret
  end
  
  def reduce_madd
    ret = self
    #p ret
    if ret.operator == :plus || ret.operator == :minus
      if !isLeaf(ret.lop) && ret.lop.operator == :mult
        if ret.operator == :plus
          op = :madd
        else
          op = :nmsub
        end
        ret = MADD.new([op,ret.lop.lop,ret.lop.rop,ret.rop,ret.get_type])
        ret.aop = ret.aop.reduce_madd
        ret.bop = ret.bop.reduce_madd
        ret.cop = ret.cop.reduce_madd
      elsif !isLeaf(ret.rop) && ret.rop.operator == :mult
        if ret.operator == :plus
          op = :madd
        else
          op = :msub
        end
        ret = MADD.new([op,ret.rop.lop,ret.rop.rop,ret.lop,ret.get_type])
        ret.aop = ret.aop.reduce_madd
        ret.bop = ret.bop.reduce_madd
        ret.cop = ret.cop.reduce_madd
      else
        ret.lop = ret.lop.reduce_madd
        ret.rop = ret.rop.reduce_madd
      end
    else
      ret.lop = ret.lop.reduce_madd
      ret.rop = ret.rop.reduce_madd if ret.rop != nil
    end
    ret
  end

  def reduce_negate_recursive
    ret = self.dup
    #p "reduce negate:",ret
    case ret.operator
    when :plus
      if !isLeaf(ret.lop) && ret.lop.operator == :uminus
        ret = Expression.new([:minus,ret.rop,ret.lop.lop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.rop) && ret.rop.operator == :uminus
        ret = Expression.new([:minus,ret.lop,ret.rop.lop])
        ret = ret.reduce_negate_recursive
      end
    when :minus
      if !isLeaf(ret.rop) && ret.rop.operator == :uminus
        ret = Expression.new([:plus,ret.lop,ret.rop.lop])
        ret = ret.reduce_negate_recursive
      end
    when :uminus
      if !isLeaf(ret.lop) && ret.lop.operator == :uminus
        ret = ret.lop
        ret = ret.reduce_negate_recursive
      end
    else
      ret.lop = ret.lop.reduce_negate_recursive if !isLeaf(ret.lop)
      ret.rop = ret.rop.reduce_negate_recursive if !isLeaf(ret.rop)
    end
    #p "reduced expression:",ret
    ret
  end
end

class FuncCall
  def reduce_madd
    self
  end
end

class String
  def reduce_madd
    self
  end
end

class MADD
  def find_mul
    ret,found = @cop.find_mul if !isLeaf(@cop)
    [ret,found]
  end

  def reduce_madd_recursive
    @aop = @aop.reduce_madd_recursive if !isLeaf(@aop)
    @bop = @bop.reduce_madd_recursive if !isLeaf(@bop)
    @cop = @cop.reduce_madd_recursive if !isLeaf(@cop)
  end

  def reduce_madd
    @aop = @aop.reduce_madd
    @bop = @bop.reduce_madd
    @cop = @cop.reduce_madd
  end

  def reduce_negate_recursive
    ret = self.dup
    #p "reduce negate:",ret
    case ret.operator
    when :madd
      if !isLeaf(ret.aop) && ret.aop.operator == :uminus
        ret = MADD.new([:msub,ret.aop.lop,ret.bop,ret.cop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.bop) && ret.bop.operator == :uminus
        ret = MADD.new([:msub,ret.aop,ret.bop.lop,ret.cop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.cop) && ret.cop.operator == :uminus
        ret = MADD.new([:nmsub,ret.aop,ret.bop,ret.cop.lop])
        ret = ret.reduce_negate_recursive
      end
    when :msub
      if !isLeaf(ret.aop) && ret.aop.operator == :uminus
        ret = MADD.new([:madd,ret.aop.lop,ret.bop,ret.cop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.bop) && ret.bop.operator == :uminus
        ret = MADD.new([:madd,ret.aop,ret.bop.lop,ret.cop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.cop) && ret.cop.operator == :uminus
        ret = MADD.new([:nmadd,ret.aop,ret.bop,ret.cop.lop])
        ret = ret.reduce_negate_recursive
      end
    when :nmadd
      if !isLeaf(ret.aop) && ret.aop.operator == :uminus
        ret = MADD.new([:nmsub,ret.aop.lop,ret.bop,ret.cop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.bop) && ret.bop.operator == :uminus
        ret = MADD.new([:nmsub,ret.aop,ret.bop.lop,ret.cop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.cop) && ret.cop.operator == :uminus
        ret = MADD.new([:msub,ret.aop,ret.bop,ret.cop.lop])
        ret = ret.reduce_negate_recursive
      end
    when :nmsub
      if !isLeaf(ret.aop) && ret.aop.operator == :uminus
        ret = MADD.new([:nmadd,ret.aop.lop,ret.bop,ret.cop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.bop) && ret.bop.operator == :uminus
        ret = MADD.new([:nmadd,ret.aop,ret.bop.lop,ret.cop])
        ret = ret.reduce_negate_recursive
      end
      if !isLeaf(ret.cop) && ret.cop.operator == :uminus
        ret = MADD.new([:madd,ret.aop,ret.bop,ret.cop.lop])
        ret = ret.reduce_negate_recursive
      end
    end
    ret.aop = ret.aop.reduce_negate_recursive if !isLeaf(ret.aop)
    ret.bop = ret.bop.reduce_negate_recursive if !isLeaf(ret.bop)
    ret.cop = ret.cop.reduce_negate_recursive if !isLeaf(ret.cop)
    #p "reduced expression:",ret
    ret
  end
end
