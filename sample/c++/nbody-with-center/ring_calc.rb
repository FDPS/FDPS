#
# ring_calc.rb
#
# create division
#
#usage:
#ruby ring_calc.rb width nx ny
#
#ruby ring_calc.rb 0.01 10 10

#
include Math

def cut_fan_area(r,x)
  if x >= r
    0
  else
    theta=acos(x/r)
    (r*r*theta -r*x*sin(theta))/2
  end
end
def band_area(x0,x1,r0,r1)
  # print "enter band", x0, " ", x1, " ", r0, " ", r1, "\n"
  # p  cut_fan_area(r1,x0)
  # p cut_fan_area(r1,x1)
  # p cut_fan_area(r0,x0)
  # p cut_fan_area(r0,x1)
  cut_fan_area(r1,x0)-cut_fan_area(r1,x1) -
    cut_fan_area(r0,x0)+cut_fan_area(r0,x1)
end

def patch_dx(y,x0,x1,r0,r1)
  yx0r1 = sqrt(r1*r1-x0*x0)
  yx0r0=0
  yx1r1=sqrt(r1*r1-x1*x1)
  yx1r0=0
  yx0r0 = sqrt(r0*r0-x0*x0) if r0>x0
  yx1r0 = sqrt(r0*r0-x1*x1) if r0>x1
  xmin=x0
  xmax=x1
  if y < yx0r0
    xmin = sqrt(r0*r0-y*y)
  end
  if y> yx1r1
    xmax = sqrt(r1*r1-y*y)
  end
  dx = xmax-xmin
  dx=0  if dx<0
  dx
end


def patch_area(y0,y1,x0,x1,r0,r1)
  n=20
  s= (patch_dx(y0,x0,x1,r0,r1)+patch_dx(y1,x0,x1,r0,r1))/2
  (n-1).times{|i|
    s+=patch_dx(y0+(i+1)*(y1-y0)/n,x0,x1,r0,r1)
  }
  s/n*(y1-y0)
end
  
def calc_x0(r0, r1, n, x1)
  s= Math::PI*(r1*r1-r0*r0)/n/4
  #  p s
  a=0
  b=x1
  x=b
  50.times{
    x= (a+b)/2
    # print "x=", x, "\n"
    # p band_area(x,x1,r0,r1)
    if band_area(x,x1,r0,r1)-s > 0
      a=x
    else
      b=x
    end
  }
  x
end

def calc_y0(x0,x1,r0, r1, n, y1)
  s= Math::PI*(r1*r1-r0*r0)/n/4
  #  p s
  yx1r0=0
  yx1r0 = sqrt(r0*r0-x1*x1) if r0>x1
  a=yx1r0
  b=y1
  y=b
  50.times{
    y= (a+b)/2
#    print "y, area=", y," ", patch_area(y,y1,x0,x1,r0,r1),"\n"
    if patch_area(y,y1,x0,x1,r0,r1)-s > 0
      a=y
    else
      b=y
    end
  }
  y
end

#term
#square
dr= ARGV[0].to_f
n= ARGV[1].to_i
m= ARGV[2].to_i
r0=1-dr/2
r1=1+dr/2
x1=r1
#limit -r1 r1 -r1 r1
print "rmax= ",r1,"\n"
#box
#print(x1,"\n")
n.times{
  x0=calc_x0(r0,r1,n,x1)
#  print "newx=", x0,"\n"
  y1 = sqrt(r1*r1-x0*x0)
#  print "ytop=", y1,"\n"
  m.times{
    y0=calc_y0(x0,x1,r0,r1,n*m,y1)
    print "box ", x0," ", x1, " ", y0, " ", y1, "\n"
    # reloc x0 y0
    # draw x1 y0
    # draw x1 y1
    # draw x0 y1
    # draw x0 y0
    y1=y0
  }
  x1=x0
}

  
  
