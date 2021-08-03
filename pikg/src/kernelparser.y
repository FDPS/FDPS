# coding: utf-8
class KernelParser
  prechigh
  left '.'
  nonassoc UMINUS
  left '*' '/'
  left '+' '-'
  left '<' '>' '<=' '>='
  left '&' '|'
  left '==' '!='
  left '&&' '||'
  left '+=' '-='
    
  preclow
  start innerkernel

  rule

  innerkernel: iodeclarations functions statements {result=Kernelprogram.new(val)}
             | iodeclarations statements           {result=Kernelprogram.new(val)}
             | functions statements                {result=Kernelprogram.new(val)}
             | statements                          {result=Kernelprogram.new(val)}
  iodeclarations: iodeclaration
                | iodeclarations iodeclaration {result = val[0]+val[1]}
  iodeclaration: iotype type varname ':' fdpsname EOL          {result = [Iodeclaration.new([val[0],val[1],val[2],val[4],nil])]}
               | iotype "static" type varname ':' fdpsname EOL {result = [Iodeclaration.new([val[0],val[2],val[3],val[5],"static"])]}
               | iotype "local" type varname EOL { result = [Iodeclaration.new([val[0],val[2],val[3],nil,"local"])]}
               | declaration

  functions : function
            | function functions {result = val[0]+val[1]}
  function : funcdeclaration statements ret_state end_state {result = [Function.new([val[0],val[1],val[2]])]}
           | funcdeclaration ret_state end_state            {result = [Function.new([val[0],[],val[1]])]}
  ret_state : "return" expression EOL {result = ReturnState.new(val[1])}
  end_state : "end" EOL
  funcdeclaration : "function" IDENT '(' operands ')' EOL {result = Funcdeclaration.new([val[1],val[3]])}

  operands : operand
           | operands ',' operand {result = val[0] + val[2]}
  operand : varname {result = [val[0]]}

  declaration: type varname EOL {result = [Iodeclaration.new(["MEMBER",val[0],val[1],nil,nil])]}

  iotype : EPI
         | EPJ
         | FORCE
         | TABLE

  type : F64VEC
       | F64
       | F32VEC
       | F32
       | F16VEC
       | F16
       | S64
       | U64
       | S32
       | U32
       | S16
       | U16

  varname : IDENT

  fdpsname : IDENT

  statements : statement
             | statement statements {result = val[0]+val[1]}

  statement : var "="  expression EOL {result = [Statement.new([val[0],val[2]])]}
            | var "="  table EOL      {result = [TableDecl.new([val[0],val[2]])]}
            | var '+=' expression EOL {result = [Statement.new([val[0],val[2],nil,:plus])]}
            | var '-=' expression EOL {result = [Statement.new([val[0],val[2],nil,:minus])]}
            | var '*=' expression EOL {result = [Statement.new([val[0],val[2],nil,:mult])]}
            | var '/=' expression EOL {result = [Statement.new([val[0],val[2],nil,:div])]}
            | pragma EOL {result = [val[0]]}
            | "if" expression EOL    {result = [IfElseState.new([:if,val[1]])]}
            | "else" EOL             {result = [IfElseState.new([:else,nil])]}
            | "elsif" expression EOL {result = [IfElseState.new([:elsif,val[1]])]}
            | "endif" EOL            {result = [IfElseState.new([:endif,nil])]}
  pragma: "#pragma" TEXT {result = Pragma.new([val[1],nil])}
        | "#pragma" TEXT options {result = Pragma.new([val[1],val[2]])}
  options: option { result = val[0]}
         | options option { result = val[0] + val[1] }
  option: TEXT { result = [val[0]]}

  expression: binary

  binary: binary '+'  binary {result = Expression.new([:plus,  val[0], val[2]])}
        | binary '-'  binary {result = Expression.new([:minus, val[0], val[2]])}
        | binary '*'  binary {result = Expression.new([:mult,  val[0], val[2]])}
        | binary '/'  binary {result = Expression.new([:div,   val[0], val[2]])}
        | binary '==' binary {result = Expression.new([:eq,    val[0], val[2]])}
        | binary '!=' binary {result = Expression.new([:neq,   val[0], val[2]])}
        | binary '&'  binary {result = Expression.new([:and,   val[0], val[2]])}
        | binary '|'  binary {result = Expression.new([:or,    val[0], val[2]])}
        | binary '&&' binary {result = Expression.new([:land,   val[0], val[2]])}
        | binary '||' binary {result = Expression.new([:lor,    val[0], val[2]])}
        | binary '>'  binary {result = Expression.new([:gt,    val[0], val[2]])}
        | binary '<'  binary {result = Expression.new([:lt,    val[0], val[2]])}
        | binary '>=' binary {result = Expression.new([:ge,    val[0], val[2]])}
        | binary '<=' binary {result = Expression.new([:le,    val[0], val[2]])}
        | '(' binary ')'     {result = val[1]}
        | funccall           {result = val[0]}
	| '-' binary =UMINUS {result = Expression.new([:uminus,val[1], nil])}
        | var
	| number

  table: '{' args '}' {result = Table.new(val[1])}

  funccall: IDENT '(' args ')' {result = FuncCall.new([val[0], val[2]])}
  args: arg
      | args ',' arg {result = val[0] + val[2]}
  arg: binary {result = [val[0]]}

  var: IDENT             {result = val[0]}
     | var '.' IDENT     {result = Expression.new([:dot,val[0],val[2]])}
     | IDENT '[' number ']' {result = Expression.new([:array,val[0],val[2]])}

  number: DEC             {result = IntegerValue.new(val[0])}
        | DEC 'l'         {result = IntegerValue.new(val[0]+val[1])}
        | DEC 's'         {result = IntegerValue.new(val[0]+val[1])}
        | DEC 'u'         {result = IntegerValue.new(val[0]+val[1])}
        | DEC "ul"        {result = IntegerValue.new(val[0]+val[1])}
        | DEC "us"        {result = IntegerValue.new(val[0]+val[1])}
        | DEC '.' DEC     {result = FloatingPoint.new(val[0]+val[1]+val[2])}
	| DEC '.' DEC "f" {result = FloatingPoint.new(val[0]+val[1]+val[2]+val[3])}
        | DEC '.' DEC "h" {result = FloatingPoint.new(val[0]+val[1]+val[2]+val[3])}



end

---- header
$lines = []
$used_q = []
def get_lineno(token)
  nline = 0
  $used_q.each{ |t|
    nline += 1 if t[0] == :EOL
    break if t[1] == token
  }
  $lines[nline][0]
end

def get_line(token)
  nline = 0
  $used_q.each{ |t|
    nline += 1 if t[0] == :EOL
    break if t[1] == token
  }
  $lines[nline][1]
end

---- inner
def vartype(string)
[[:F64VEC,:F64,:F32VEC,:F32,:F16VEC,:F16,:S64,:U64,:S32,:U32,:S16,:U16][["F64vec","F64","F32vec","F32","F16vec","F16","S64","U64","S32","U32","S16","U16"].index(string)],string]
end
def iotype(string)
[[:EPI,:EPJ,:FORCE,:MEMBER,:ACCUM,:TABLE][["EPI","EPJ","FORCE","MEMBER","ACCUM","TABLE"].index(string)], string]
end
def dim(string)
  [[:x,:y,:z,:xx,:yy,:zz,:xy,:xz,:yz,:yx,:zx,:zy][["x","y","z","xx","yy","zz","xy","xz","yz","yx","zx","zy"].index(string)],string]
end

def parse(filename)
  @q=[]
  open(filename){ |f|
    count = 0
    f.each_line{|str|
      a=str.chomp.split(/(\/\/|\s|:|\=\=|\!\=|\+\=|-\=|\*\=|\&\&|\|\||\&|\||\+|-|\*|\/|\=|\(|\)|,|\.|>\=|<\=|>|<|\[|\]|;|~|\^)/).select{|s| s=~/\S+/}
      symbols = /^(\=\=|\!\=|\&\&|\|\||\&|\||\+|-|\*|\/|\(|\)|,|\.|>\=|<\=|>|<|\[|\]|\{|\})$/
      if a == [] || a[0] == "//"
        next
      end
      count += 1
      $lines.push([count,str])
      if a[0] =~ /(EPI|EPJ|FORCE|MEMBER|ACCUM|TABLE)/ && a[1] != '.'
        io = a[0]
        @q << iotype(a.shift) # iotype
      end
      if a[0] == "static"
        @q << ["static",a.shift]
      elsif a[0] == "local"
        @q << ["local",a.shift]
      end

      if a[0] =~ /(F(64|32|16)|F(64|32|16)vec|(U|S)(64|32|16))/
        @q << vartype(a.shift) # type
        @q << [:IDENT, a.shift] # varname
        if io != nil && a[0] == ":"
          @q << [":",a.shift] # :
          @q << [:IDENT, a.shift] # fdpsname
        end
      elsif a[0] == "function"
        #      print "function decl"
        @q << ["function", a.shift]
        @q << [:IDENT, a.shift]
        @q << ["(", a.shift]
        a.each{|x|
          if x =~ /(F(64|32|16)|F(64|32|16)vec|(U|S)(64|32|16))/
            @q << vartype(x)
          elsif x  == ','
            @q << [',', x]
          elsif x  == ')'
            @q << [')', x]
	  elsif x == "//"
	       break
          else
            @q << [:IDENT, x]
          end
        }
      elsif a[0] == "return"
        @q << ["return", a.shift]
        a.each{ |x|
          if x =~ symbols
            @q << [x, x]
	  elsif x =~ /^\d+(f|h|s|l|u|us|ul)?$/
            @q << [:DEC, x]
	  elsif x == "//"
		break
          else
            @q << [:IDENT, x]
          end
        }
      elsif a[0] == "end"
        @q << ["end", a.shift]
      elsif a[0] == "#pragma"
        #p a
        @q << ["#pragma",a.shift]
        a.each{|x|
	  if x == "//"
	       break
	  else
	    @q << [:TEXT,x]
          end
	       
        }
      elsif a[0] == "if"
        @q << ["if","if"]
        a.shift
        a.each{|x|
	  if x =~ symbols
	    @q << [x,x]
	  elsif x =~ /^\d+(f|h|s|l|u|us|ul)?$/
            @q << [:DEC, x]
	  elsif x == "//"
	    break
          else
            @q << [:IDENT, x]
	  end
        }
      elsif a[0] == "elsif"
        @q << ["elsif","elsif"]
        a.shift
        a.each{|x|
	  if x =~ symbols
	    @q << [x,x]
	  elsif x =~ /^\d+(f|h|s|l|u|us|ul)?$/
	    @q << [:DEC, x]
	  elsif x == "//"
	    break
	  else
            @q << [:IDENT, x]
	  end
        }
      elsif a[0] == "else"
        @q << ["else","else"]
      elsif a[0] == "endif"
        @q << ["endif","endif"]
      elsif a[1] =~ /(\=|\+\=|-\=)/
        #print "statement \n"
        @q << [:IDENT,a.shift]
        #p a
        @q << [a[0], a[0]]
        a.shift
        a.each{|x|
          if x =~ symbols
            @q << [x,x]
	  elsif x =~/^\d+(f|h|s|l|u|us|ul)?$/
            @q << [:DEC,x]
          elsif x == "//"
	    break   
          else
            @q << [:IDENT,x]
          end
	}
      elsif a[3] =~ /(\=|\+\=|-\=)/
        #print "statement \n"
        @q << [:IDENT,a.shift] # var
        @q << [a[0], a[0]]    #.
        a.shift
        @q << [:IDENT,a.shift]
        @q << [a[0], a[0]]    # =|+=|-=|*=|/=
        a.shift
        a.each{|x|
          if x =~ symbols
            @q << [x,x]
	  elsif x =~/^\d+(f|h|s|l|u|us|ul)?$/
            @q << [:DEC,x]
	  elsif x == "break"
	    break
          else
            @q << [:IDENT,x]
          end
        }
      else
        warn "error: unsupported DSL description"
        warn "  line #{f.lineno}: #{str}"
        abort
      end
      @q << [:EOL,:EOL]
    }
    @q << nil
    do_parse
  }
end

def next_token
  tmp = @q.shift
  $used_q << tmp
  tmp
end


def on_error(id,val,stack)
  #lineno = get_lineno(val)
  #line   = get_line(val)
  line0 = []
  $used_q.reverse_each{ |v|
    break if v[1] == :EOL
    line0.push(v[1])
  }
  line0 = "parse error: " + line0.reverse.join('')

  line1 = []
  @q.each{ |v|
    break if v[1] == :EOL
    line1.push(v[1])
  }
  line1 = line1.join('')

  #warning = "parse error :line #{lineno}: #{line}"
  warning = line0 + line1 + "\n"
  for i in 1..line0.length
    warning += " "
  end
  warning += "^"
  abort warning
end
