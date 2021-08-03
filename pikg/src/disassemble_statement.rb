require_relative "common.rb"

class Kernelprogram
  def disassemble_statement(h = $varhash)
    new_s = []
    if @statements != []
      @statements.each{ |s|
        if isStatement(s) && h[get_name(s)][3] != "local"
          new_s += s.disassemble
        else
          new_s += [s.dup]
        end
      }
      @statements = new_s
    end
  end
end

class Statement
  def disassemble
    new_s = []
    if isLeaf(@expression)
      new_s += [Statement.new([@name.dup,@expression.dup,@type.dup])]
    else
      get_type() if @type == nil
      new_s += @expression.disassemble(@name,@type)
    end
    #p new_s
    new_s
  end
end

class Expression
  def disassemble(name,type)
    abort "type is nil at disassemble" if type == nil
    ops = []
    new_s = []
    if !(@operator == :dot || @operator == :array)
      [@lop,@rop].each{ |op|
        next if op == nil 
        if !isLeaf(op)
          new_name = add_new_tmpvar(type)
          new_s += Statement.new([new_name,op.dup,op.get_type.dup]).disassemble if op != nil
        #p new_s
        else
          new_name = op.dup
        end
        ops += [new_name]
      }
      new_s += [Statement.new([name,Expression.new([@operator,ops.shift,ops.shift,type]),type])]
    end
    new_s
  end
end

class MADD
  def disassemble(name,type)
    abort "type is nil at disassemble" if type == nil
    ops = []
    new_s = []
    [@aop,@bop,@cop].each{ |op|
      if !isLeaf(op)
        new_name = add_new_tmpvar(type)
        new_s += Statement.new([new_name,op.dup,op.type]).disassemble
      else
        new_name = op.dup
      end
      ops += [new_name]
    }
    new_s += [Statement.new([name,MADD.new([@operator,ops.shift,ops.shift,ops.shift,type]),type])]
    new_s
  end
end


