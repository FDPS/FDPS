#
# add_cominfo_if_fortran.rb
#
# create new  (Version 7)  FDPS_module_blueprint.F90
# from old one (version 6)
#

def process_interface_additions
  print <<EOF
      procedure :: ci_initialize
      procedure :: ci_set_communicator
      procedure :: ci_delete
      procedure :: ci_create
      procedure :: ci_split
EOF
end
def process_privates_additions
  print <<EOF
      private :: ci_initialize
      private :: ci_set_communicator
      private :: ci_delete
      private :: ci_create
      private :: ci_split
EOF
end

def process_f_interface_additions
  print <<EOF
      function fdps_f_ci_initialize(comm) bind (c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_f_ci_initialize
         integer(kind=c_int), value, intent(in) :: comm
      end function fdps_f_ci_initialize

      function fdps_f_ci_set_communicator(ci, comm) bind (c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_f_ci_set_communicator
         integer(kind=c_int), value, intent(in) :: ci
         integer(kind=c_int), value, intent(in) :: comm
      end function fdps_f_ci_set_communicator

      function fdps_f_ci_delete(ci) bind (c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_f_ci_delete
         integer(kind=c_int), value, intent(in) :: ci
      end function fdps_f_ci_delete

      function fdps_f_ci_create(ci,n,rank) bind (c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_f_ci_create
         integer(kind=c_int), value, intent(in) :: ci
         integer(kind=c_int), value, intent(in) :: n
         integer(kind=c_int), dimension(n), intent(in) :: rank
      end function fdps_f_ci_create

      function fdps_f_ci_split(ci,color,key) bind (c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_f_ci_split
         integer(kind=c_int), value, intent(in) :: ci
         integer(kind=c_int), value, intent(in) :: color
         integer(kind=c_int), value, intent(in) :: key
      end function fdps_f_ci_split

EOF
end

def process_f_body_additions
  print <<EOF
      function ci_initialize(this, comm)
         implicit none
         class(FDPS_controller) :: this
         integer(kind=c_int) :: ci_initialize
         integer(kind=c_int), intent(in) :: comm
         ci_initialize = fdps_f_ci_initialize(comm)
      end function ci_initialize

      function ci_set_communicator(this, ci, comm)
         implicit none
         class(FDPS_controller) :: this
         integer(kind=c_int) :: ci_set_communicator
         integer(kind=c_int), value, intent(in) :: ci
         integer(kind=c_int), value, intent(in) :: comm
         ci_set_communicator = fdps_f_ci_set_communicator(ci, comm)
      end function ci_set_communicator

      function ci_delete(this,ci) 
         implicit none
         class(FDPS_controller) :: this
         integer(kind=c_int) :: ci_delete
         integer(kind=c_int), value, intent(in) :: ci
         ci_delete = fdps_f_ci_delete(ci)
      end function ci_delete

      function ci_create(this, ci,n,rank) 
         implicit none
         class(FDPS_controller) :: this
         integer(kind=c_int) :: ci_create
         integer(kind=c_int), value, intent(in) :: ci
         integer(kind=c_int), value, intent(in) :: n
         integer(kind=c_int), dimension(n), intent(in) :: rank
         ci_create = fdps_f_ci_create(ci,n,rank)
      end function ci_create

      function ci_split(this, ci, color, key)
         implicit none
         class(FDPS_controller) :: this
         integer(kind=c_int) :: ci_split
         integer(kind=c_int), value, intent(in) :: ci
         integer(kind=c_int), value, intent(in) :: color
         integer(kind=c_int), value, intent(in) :: key
         ci_split = fdps_f_ci_split(ci, color, key)
      end function ci_split

EOF
end

def process_interface(s)
  while s !=nil
    args=[]
    a=s.chomp.split
    ss=nil
    if a[0] == "procedure"
      header = "      "+a[0]+" "+a[1]+" "
      funcname = a[2]
    elsif a[0] == "procedure,"
      header = "      "+a[0]+" "+a[1]+" "+a[2]+" "
      funcname = a[3]
    elsif    a[0] == "generic"
      header = "      "+a[0]+" "+a[1]+" ci_"+a[2]+ " => "
      funcname = a[4..5].join(" ")
      ss=gets.chomp
      while ss.index("                       ")
        args.push ss
        ss=gets.chomp
      end
      ss += "\n"
    end
    print s
    args.each{|x|       print x,"\n"}
    print header, "ci_"+funcname, "\n"
    args.each{|x|
      a= x.chomp.split
      print "                                  "
      print "ci_"+a.join(" "), "\n"
    }
    return if funcname =~ /barrier/
    if ss
      s=ss
    else
      s=gets
    end
  end
end

def process_privates(s)
  while s !=nil
    a=s.chomp.split
    if a[0] == "private"
      header = "   "+a[0]+" "+a[1]+" "
      funcname = a[2]
    end
    print s
    print header, "ci_"+funcname, "\n"
    return if funcname =~ /barrier/
    s=gets
  end
end


def process_f_interface(s)
  while s !=nil
    ss=nil
    if s =~/function|subroutine/
      args=[]
      a=s.chomp.split
      header = "      "+a[0]+" "
      funcname = a[1..(a.size)].join(" ")
      ss=gets.chomp
      while ss =~/\S/
        args.push ss
        ss=gets.chomp
      end
      print s
      args.each{|x|       print x,"\n"}
      print "\n"
      funcname =~/(\S+)\(/
      fname=$1
      funcname.sub!("fdps_", "fdps_f_ci_")
      unless funcname.gsub!("()", "(ci)")
        funcname.sub!(/\((\w*)/){"(ci, "+$1}
      end

      
      print header, funcname,"\n"
      args[0..(args.size()-2)].each{|x|
        print x.gsub("fdps_", "fdps_f_ci_"),"\n"
      }
      print "         integer(kind=c_int), value, intent(in) :: ci\n"
      print args[args.size()-1].gsub("fdps_", "fdps_f_ci_"),"\n\n"
      return if funcname =~ /barrier/
    else
      print s
    end
    if ss
      s=ss
    else
      s=gets
    end
  end
end


def process_f_body(s)
  while s !=nil
    ss=nil
    if s =~/function|subroutine/
      args=[]
      a=s.chomp.split
      header = "   "+a[0]+" "
      funcname = a[1..(a.size)].join(" ")
      end_reached = false
      until end_reached
        ss=gets.chomp
        args.push ss
        end_reached = ss =~/end (subroutine|function)/
      end
      ss=gets.chomp
      print s
      args.each{|x|       print x,"\n"}
      print "\n"
      funcname =~/(\S+)\(/
      fname=$1
      funcname.sub!(/\(this/, "(this, ci")
      
      print header, "ci_"+funcname,"\n"
      args[0..(args.size()-2)].each{|x|
        x.gsub!(fname, "ci_"+fname)
        x.gsub!("fdps_ci_", "fdps_f_ci_")
        unless x.gsub!(/fdps_f_ci_(\S+)\(\)/){"fdps_f_ci_"+$1+"(ci)"}
          x.gsub!(/fdps_f_ci_(\S+)\(/){"fdps_f_ci_"+$1+"(ci,"}
        end
        
        print x, "\n"
        if x =~/class\(FDPS/
          print "      integer(kind=c_int), intent(IN) :: ci\n"
        end
      }
      a = args[args.size()-1].split
      print "   ", a[0], " ", a[1], " ci_", a[2],"\n\n"
      return if funcname =~ /barrier/
    else
      print s
    end
    if ss
      s=ss
    else
      s=gets
    end
  end
end


in_comm = 0
header_printed = false
while s=gets
  in_comm = 1 if s =~ /procedure :: get_rank /
  in_comm = 2 if s =~ /private :: get_rank /
  in_comm = 3 if s =~ /function fdps_get_rank/
  in_comm = 4 if s =~ /function get_rank/
  if in_comm == 0
    print s.chomp, "\n"
  elsif  in_comm == 1
    process_interface(s)
    process_interface_additions()
    in_comm=0
  elsif  in_comm == 2
    process_privates(s)
    process_privates_additions()
    in_comm=0
  elsif  in_comm == 3
    process_f_interface(s)
    process_f_interface_additions()
    in_comm=0
  elsif  in_comm == 4
    process_f_body(s)
    process_f_body_additions()
    in_comm=0
  end
end


    
