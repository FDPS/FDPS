def byte_count(type)
  case type
  when "F64" then
    byte = 64
  when "F32" then
    byte = 32
  when "S64" then
    byte = 64
  when "S32" then
    byte = 32
  when "U64" then
    byte = 64
  when "U32" then
    byte = 32
  when "F64vec2" then
    byte = 128
  when "F32vec2" then
    byte = 64
  when "F64vec3" then
    byte = 192
  when "F32vec3" then
    byte = 96
  when "F64vec4" then
    byte = 256
  when "F32vec4" then
    byte = 128
  when "F64vec" then
    byte = 192
  when "F32vec" then
    byte = 96
  else
    abort "error: unsupported scalar type of #{@type} for A64FX at byte_count"
  end
  byte
end

def find_loop_fission_load_store_vars(ss)
  # build related vars list
  ls = []
  l = []
  ss.each{ |s|
    if s.class == Pragma && s.name == "statement" && s.option == ["loop_fission_point"]
      ls += [l]
      l = []
    elsif isStatement(s)
      l += [[get_name(s),s.expression.get_related_variable]]
    elsif s.class == IfElseState
      l += [[nil,s.expression.get_related_variable]] if s.expression != nil
    elsif s.class == ConditionalBranch
      l += s.get_related_variable
    end
  }
  ls += [l]

  seed = []
  ls.each{ |l|
    tmp = []
    l.each{ |v|
      tmp += v[1]
    }
    seed += [tmp.sort.uniq]
  }
  ls.zip(seed){ |l1,l2|
    l1.each{ |l|
      l2.delete(l[0])
    }
  }

  result = []
  ls.each{|l|
    tmp = []
    l.each{ |v|
      tmp += [v[0]]
    }
    result += [tmp.uniq]
  }

  ret = [[[],seed.shift]]
  count = 0
  seed.zip(result){ |s,r|
    lload  = []
    lstore = []
    seed.each{ |vv|
      vv.each{ |v|
        lload += [v] if s.index(v) != nil
      }
    }
    #seed.shift
    seed.each_with_index{ |vv,i|
      if i >= count
        vv.each{ |v|
          lstore += [v] if r.index(v) != nil
        }
      end
    }
    #result.shift
    ret += [[lstore.uniq,lload.uniq]]
    count += 1
  }

  #abort "find_loop_fission_load_store_vars is testing"
  ret
end

