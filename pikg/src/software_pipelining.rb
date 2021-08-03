class Statement
  def copy_for_swpl(n,index,map)
    name = @name.copy_for_swpl(n,index,map)
    exp  = @expression.copy_for_swpl(n,index,map)
    Statement.new([name,exp,@type])
  end

  def find_unloaded_var(map)
    ret = []
    ret += @name.find_unloaded_var(map)
    ret += @expression.find_unloaded_var(map)
    ret.uniq
  end
end

class Expression
  def copy_for_swpl(n,index,map)
    lop = @lop.copy_for_swpl(n,index,map) if @lop
    rop = @rop.copy_for_swpl(n,index,map) if @rop != nil
    Expression.new([@operator,lop,rop,@type])
  end
  def find_unloaded_var(map)
    ret = []
    ret += @lop.find_unloaded_var(map)
    ret += @rop.find_unloaded_var(map)
    ret.uniq
  end
end

class NonSimdExp
  def copy_for_swpl(n,index,map)
    lop = @lop.copy_for_swpl(n,index,map)
    rop = @rop.copy_for_swpl(n,index,map) if @rop != nil
    NonSimdExp.new([@operator,lop,rop,@type])
  end
end

class MADD
  def copy_for_swpl(n,index,map)
    aop = @aop.copy_for_swpl(n,index,map)
    bop = @bop.copy_for_swpl(n,index,map)
    cop = @cop.copy_for_swpl(n,index,map)
    MADD.new([@operator,aop,bop,cop,@type])
  end
  def find_unloaded_var(map)
    ret = []
    ret += @aop.find_unloaded_var(map)
    ret += @bop.find_unloaded_var(map)
    ret += @cop.find_unloaded_var(map)
    ret.uniq
  end
end

class FuncCall
  def copy_for_swpl(n,index,map)
    ops = []
    if @name == "table"
      ops.push(@ops[0])
      ops.push(@ops[1].copy_for_swpl(n,index,map))
    else
      @ops.each{ |op|
        ops.push(op.copy_for_swpl(n,index,map))
      }
    end
    FuncCall.new([@name,ops,@type])
  end
  def find_unloaded_var(map)
    ret = []
    @ops.each{|op|
      ret += op.find_unloaded_var(map)
    }
    ret.uniq
  end
end

class FloatingPoint
  def copy_for_swpl(n,index,map)
    self
  end
end

class IntegerValue
  def copy_for_swpl(n,index,map)
    self
  end
end

class String
  def copy_for_swpl(n,index,map)
    if (map.index(self) || self == index) && $varhash[self][0] != "EPI" && $varhash[self][0] != "MEMBER"
      self + "_swpl#{n}"
    elsif self == index
      NonSimdExp.new([:plus,self,"#{n}","S32"])
    else
      self
    end
  end
  def find_unloaded_var(map)
    ret = []
    if map[self] != nil
      ret.push(map[self])
      map[self] = nil
    end
    ret
  end
end

class Declaration
  def copy_for_swpl(n,index,map)
    Declaration.new([@type,@name.copy_for_swpl(n,index,map)])
  end
end

class Duplicate
  def copy_for_swpl(n,index,map)
    name = @name.copy_for_swpl(n,index,map)
    exp  = @expression.copy_for_swpl(n,index,map)
    Duplicate.new([name,exp,@type])
  end
  def find_unloaded_var(map)
    []
  end
end

class IfElseState
  def copy_for_swpl(n,index,map)
    exp = nil
    exp = @expression.copy_for_swpl(n,index,map) if @expression != nil
    IfElseState.new([@operator,exp])
  end
  def find_unloaded_var(map)
    ret = []
    ret += @expression.find_unloaded_var(map) if @expression != nil
    ret
  end
end

class Merge
  def copy_for_swpl(n,index,map)
    op1 = @op1.copy_for_swpl(n,index,map)
    op2 = @op2.copy_for_swpl(n,index,map)
    Merge.new([op1,op2,@type])
  end
  def find_unloaded_var(map)
    []
  end
end

class ConditionalBranch
  def copy_for_swpl(n,index,map)
    conditions = []
    @conditions.each{ |c|
      conditions.push(c.copy_for_swpl(n,index,map))
    }
    bodies = []
    @bodies.each{ |b|
      body = []
      b.each{ |s|
        body.push(s.copy_for_swpl(n,index,map))
      }
      bodies.push(body)
    }
    ConditionalBranch.new([conditions,bodies])
  end

  def find_unloaded_var(map)
    ret = []
    @conditions.each{ |c|
      ret += c.find_unloaded_var(map)
    }
    @bodies.each{ |b|
      b.each{ |s|
        ret += s.find_unloaded_var(map)
      }
    }
    ret
  end
end

class StoreState
  def copy_for_swpl(n,index,map)
    dest = @dest.copy_for_swpl(n,index,map)
    src = @src.copy_for_swpl(n,index,map)
    StoreState.new([dest,src,@type])
  end
  def find_unloaded_var(map)
    []
  end
end

class LoadState
  def copy_for_swpl(n,index,map)
    dest = @dest.copy_for_swpl(n,index,map)
    src  = @src.copy_for_swpl(n,index,map)
    LoadState.new([dest,src,@type])
  end
end

class PointerOf
  def copy_for_swpl(n,index,map)
    exp = @exp.copy_for_swpl(n,index,map)
    PointerOf.new([@type,exp])
  end
  def find_unloaded_var(map)
    ret = []
    ret += @exp.find_unloaded_map(map)
    ret
  end
end

class Nest
  attr_accessor :body
  def initialize(x)
    @body = x
  end
  def convert_to_code(conversion_type)
    ret = "{\n"
    body.each{|s|
      ret += s.convert_to_code(conversion_type) + "\n"
    }
    ret += "}\n"
    ret
  end
end

def software_pipelining(orig,nstage = $swpl_stage)
  warn "software pipelining is applied"
  result = []
  orig.statements.each{ |s|
    result.push(get_name(s)) if s.class == StoreState || (isStatement(s) && $varhash[get_name(s)][0] == "FORCE")
  }
  result.uniq!

  tmpvar_map = generate_related_map(result,orig.statements)
  tmpvar_map += result
  tmpvar_map.uniq!

  pipeline = []
  decl = []

  for i in 0...nstage
    tmpvar_map.each{ |v|
      iotype = $varhash[v][0]
      type = $varhash[v][1]
      decl.push(Declaration.new([type,v]).copy_for_swpl(i,orig.index,tmpvar_map)) if iotype != "EPI"
    }
  end

  for i in 0...nstage
    index = orig.index+"_swpl#{i}"
    $varhash[index] = [nil,"S32",nil,nil]
    pipeline[i] = []
    decl.push(NonSimdDecl.new(["S32",index]))
    lstate = Hash.new()
    orig.statements.each{ |s|
      s.name.get_type if isStatement(s)
      tmp = s.copy_for_swpl(i,orig.index,tmpvar_map)
      iotype = $varhash[get_name(s)][0]
      $varhash[get_name(tmp)] = $varhash[get_name(s)].dup if $varhash[get_name(tmp)] == nil
      $varhash[get_name(tmp)][0] = "SWPL_TMP" if iotype == "FORCE" || iotype == "EPI"
      pipeline[i].push(tmp) if s.class != Declaration
    }
    pipeline[i].push(NonSimdState.new([index,Expression.new([:plus,index,"#{nstage}","S32"]),"S32"]))
  end

  length = pipeline[0].length
  if pipeline[0].length < nstage
    warn "too short loop for #{nstage} stage SWPL. SWPL is not applied"
    return orig
  end

  head = []
  heads = []
  for i in 0...nstage-1
    for j in 0..i
      heads[j] = [] if heads[j] == nil
      heads[j].push(pipeline[j][0])
      head.push(pipeline[j].shift)
    end
  end
  heads[nstage-1] = []
  abort "making head failed. length #{pipeline[nstage-1].length} != #{length}" if pipeline[nstage-1].length != length

  body = []
  for i in 0...length
    for j in 0...nstage
      if i < pipeline[j].length
        body.push(pipeline[j][i])
      else
        body.push(heads[j][ i - pipeline[j].length ])
      end
    end
  end

  tail = []
  for i in 0...nstage
    tail.push(pipeline[i].shift) if pipeline[i] != []
  end

  ret = Nest.new([])
  ret.body += decl
  for i in 0...nstage
    result.each{ |v|
      iotype = $varhash[v][0]
      type = $varhash[v][1]
      name = v.copy_for_swpl(i,orig.index,tmpvar_map)
      if iotype == "FORCE"
        name = v.copy_for_swpl(i,orig.index,tmpvar_map)
        zero = FloatingPoint.zero(type)
        if type =~ /vec/
          type = type.delete("vec")
          ret.body += [Statement.new([Expression.new([:dot,name,"x",type]),zero,type])]
          ret.body += [Statement.new([Expression.new([:dot,name,"y",type]),zero,type])]
          ret.body += [Statement.new([Expression.new([:dot,name,"z",type]),zero,type])]
        else
          ret.body += [Statement.new([name,zero,type,type])]
        end
      end
    }
  end
  for i in 0...nstage
    ret.body.push(NonSimdState.new([orig.index+"_swpl#{i}","#{i}","S32"]))
  end
  ret.body += head
  ret.body.push("//body")
  nloop = (orig.loop_end.to_i / nstage)
  for i in 0...nloop
    ret.body.push("//loop #{i} / #{nloop}")
    ret.body += body
  end
  #ret.body.push("//tail")
  #ret.body += tail

  for i in 0...nstage
    result.each{ |v|
      iotype = $varhash[v][0]
      if iotype == "FORCE"
        type = $varhash[v][1]
        name = v.copy_for_swpl(i,orig.index,tmpvar_map)
        if type =~ /vec/
          type = type.delete("vec")
          ret.body += [Statement.new([Expression.new([:dot,v,"x",type]),Expression.new([:plus,Expression.new([:dot,v,"x",type]),Expression.new([:dot,name,"x",type]),type]),type])]
          ret.body += [Statement.new([Expression.new([:dot,v,"y",type]),Expression.new([:plus,Expression.new([:dot,v,"y",type]),Expression.new([:dot,name,"y",type]),type]),type])]
          ret.body += [Statement.new([Expression.new([:dot,v,"z",type]),Expression.new([:plus,Expression.new([:dot,v,"z",type]),Expression.new([:dot,name,"z",type]),type]),type])]
        else
          ret.body += [Statement.new([v,Expression.new([:plus,v,name,type]),type])]
        end
      end
    }
  end
  
  #print ret.convert_to_code("reference") + "\n"

  #abort "software_pipelining is testing. remove --software-pipelining option to ignore"

  ret
end


def loop_unroll(orig,accum_hash,nstage = $swpl_stage)
  #warn "loop unrolling is applied"
  result = []
  orig.statements.each{ |s|
    result.push(get_name(s)) if s.class == StoreState || (isStatement(s) && $varhash[get_name(s)][0] == "FORCE")
    if s.class == ConditionalBranch
      s.get_related_variable.each{|v|
        result.push(v[0]) if $varhash[v[0]] != nil && $varhash[v[0]][0] == "FORCE"
      }
    end
  }
  result.uniq!

  tmpvar_map = generate_related_map(result,orig.statements)
  tmpvar_map += result
  tmpvar_map.uniq!

  pipeline = []
  decl = []

  for i in 0...nstage
    tmpvar_map.each{ |v|
      iotype = $varhash[v][0]
      type = $varhash[v][1]
      decl.push(Declaration.new([type,v]).copy_for_swpl(i,orig.index,tmpvar_map)) if iotype != "EPI" && iotype != "MEMBER"
    }
  end

  loop_tmp = Loop.new([orig.index,orig.loop_beg,orig.loop_end,orig.interval*nstage,[],orig.option])
  for i in 0...nstage
    index = orig.index+"_swpl#{i}"
    $varhash[index] = [nil,"S32",nil,nil]
    pipeline[i] = []
    decl.push(NonSimdDecl.new(["S32",index]))
    loop_tmp.statements.push(NonSimdState.new([index,NonSimdExp.new([:plus,orig.index,"#{i}","S32"])]))
    lstate = Hash.new()
    orig.statements.each{ |s|
      pipeline[i] += unroll_statement_recursive(s,i,orig.index,nstage,tmpvar_map)
    }
  end

  length = pipeline[0].length

  ret = Nest.new([])
  ret.body += decl
  #initialize unrolled FORCE var
  tmpvar_map.each{|v|
    iotype = $varhash[v][0]
    type   = $varhash[v][1]
    if iotype == "FORCE"
      type_single = get_single_element_type(type)
      get_vector_elements(type).zip(accum_hash[v]){ |dim,op|
        for i in 0...nstage
          name = v + "_swpl#{i}"
          name = Expression.new([:dot,name,dim,type_single]) if type =~ /vec/
          ret.body += [Duplicate.new([name,get_initial_value(op,type_single),type_single])]
        end
      }
    end
  }
  for i in 0...length
    for j in 0...nstage
      loop_tmp.statements.push(pipeline[j][i])
    end
  end
  ret.body += [loop_tmp]
  # accumulate unrolled FORCE var to original FORCE var
  tmpvar_map.each{|v|
    iotype = $varhash[v][0]
    type   = $varhash[v][1]
    if iotype == "FORCE"
      type_single = get_single_element_type(type)
      get_vector_elements(type).zip(accum_hash[v]){ |dim,op|
        for i in 0...nstage
          src = v + "_swpl#{i}"
          src = Expression.new([:dot,src,dim,type_single]) if type =~ /vec/
          name = v
          name = Expression.new([:dot,name,dim,type_single]) if type =~ /vec/
          if op == "max" || op == "min"
            ret.body += [Statement.new([name,FuncCall.new([op,[name,src],type_single])])]
          elsif op == :plus || op == :minus
            ret.body += [Statement.new([name,Expression.new([:plus,name,src,type_single]),type_single])]
          elsif op == :mult || :div
            ret.body += [Statement.new([name,Expression.new([:mult,name,src,type_single]),type_single])]
          else
            abort "unsupported accumulate operator #{op}"
          end
        end
      }
    end
  }

  ret
end

def unroll_statement_recursive(s,i,index,unroll_stage=$unroll_stage,tmpvar_map)
  ret = Array.new

  s.name.get_type if isStatement(s)
  if s.class == ConditionalBranch
    cb = ConditionalBranch.new([[],[]])
    s.conditions.each{ |c|
      cb.conditions.push(c.copy_for_swpl(i,index,tmpvar_map))
    }
    s.bodies.each{ |b|
      tmp = Array.new
      b.each{ |bs|
        tmp += unroll_statement_recursive(bs,i,index,unroll_stage,tmpvar_map)
      }
      cb.bodies.push(tmp)
    }
    ret.push(cb)
  else
    tmp = s.copy_for_swpl(i,index,tmpvar_map)
    iotype = $varhash[get_name(s)][0]
    $varhash[get_name(tmp)] = $varhash[get_name(s)].dup if $varhash[get_name(tmp)] == nil
    ret.push(tmp) if s.class != Declaration
  end

  ret
end
