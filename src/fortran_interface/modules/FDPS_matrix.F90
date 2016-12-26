!==============================
!   MODULE: FDPS matrix
!==============================
module fdps_matrix
   use, intrinsic :: iso_c_binding
   implicit none

   !**** PS::F32mat
   type, public, bind(c) :: fdps_f32mat
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      real(kind=c_float) :: xx,yy,zz,xy,xz,yz 
#else
      real(kind=c_float) :: xx,yy,xy
#endif
   end type fdps_f32mat

   !**** PS::F64mat
   type, public, bind(c) :: fdps_f64mat
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      real(kind=c_double) :: xx,yy,zz,xy,xz,yz 
#else
      real(kind=c_double) :: xx,yy,xy
#endif
   end type fdps_f64mat
   
   !* Private routines
   private :: set_m32_i
   private :: set_m32_f
   private :: set_m32_d
   private :: set_m32_m32
   private :: set_m32_m64
   private :: set_m64_i
   private :: set_m64_f
   private :: set_m64_d
   private :: set_m64_m32
   private :: set_m64_m64

   private :: add_m32_m32
   private :: add_m32_m64
   private :: add_m64_m32
   private :: add_m64_m64
   private :: pos_m32
   private :: pos_m64

   private :: sub_m32_m32
   private :: sub_m32_m64
   private :: sub_m64_m32
   private :: sub_m64_m64
   private :: neg_m32
   private :: neg_m64

   private :: mul_m32_i
   private :: mul_m32_f
   private :: mul_m32_d
   private :: mul_i_m32
   private :: mul_f_m32
   private :: mul_d_m32
   private :: mul_m64_i
   private :: mul_m64_f
   private :: mul_m64_d
   private :: mul_i_m64
   private :: mul_f_m64
   private :: mul_d_m64
   private :: mul_m32_m32
   private :: mul_m64_m32
   private :: mul_m32_m64
   private :: mul_m64_m64

   private :: div_m32_i
   private :: div_m32_f
   private :: div_m32_d
   private :: div_m64_i
   private :: div_m64_f
   private :: div_m64_d


   !* Declare interface operators
   interface assignment (=)
      module procedure set_m32_i
      module procedure set_m32_f
      module procedure set_m32_d
      module procedure set_m32_m32
      module procedure set_m32_m64
      module procedure set_m64_i
      module procedure set_m64_f
      module procedure set_m64_d
      module procedure set_m64_m32
      module procedure set_m64_m64
   end interface

   interface operator (+)
      module procedure add_m32_m32
      module procedure add_m32_m64
      module procedure add_m64_m32
      module procedure add_m64_m64
      module procedure pos_m32
      module procedure pos_m64
   end interface

   interface operator (-)
      module procedure sub_m32_m32
      module procedure sub_m32_m64
      module procedure sub_m64_m32
      module procedure sub_m64_m64
      module procedure neg_m32
      module procedure neg_m64
   end interface 

   interface operator (*)
      !* scalar product
      module procedure mul_m32_i
      module procedure mul_m32_f
      module procedure mul_m32_d
      module procedure mul_i_m32
      module procedure mul_f_m32
      module procedure mul_d_m32
      module procedure mul_m64_i
      module procedure mul_m64_f
      module procedure mul_m64_d
      module procedure mul_i_m64
      module procedure mul_f_m64
      module procedure mul_d_m64
      module procedure mul_m32_m32
      module procedure mul_m64_m32
      module procedure mul_m32_m64
      module procedure mul_m64_m64
   end interface

   interface operator (/)
      module procedure div_m32_i
      module procedure div_m32_f
      module procedure div_m32_d
      module procedure div_m64_i
      module procedure div_m64_f
      module procedure div_m64_d
   end interface

   contains

   !--------------------------------------------------------------------
   elemental pure subroutine set_m32_i(m,i)
      type(fdps_f32mat), intent(inout) :: m
      integer(kind=c_int), intent(in) :: i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      m%xx = real(i,kind=c_float)
      m%yy = real(i,kind=c_float)
      m%zz = real(i,kind=c_float)
      m%xy = real(i,kind=c_float)
      m%xz = real(i,kind=c_float)
      m%yz = real(i,kind=c_float)
#else
      m%xx = real(i,kind=c_float)
      m%yy = real(i,kind=c_float)
      m%xy = real(i,kind=c_float)
#endif
   end subroutine set_m32_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_m32_f(m,f)
      type(fdps_f32mat), intent(inout) :: m
      real(kind=c_float), intent(in) :: f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      m%xx = f
      m%yy = f
      m%zz = f
      m%xy = f
      m%xz = f
      m%yz = f
#else
      m%xx = f
      m%yy = f
      m%xy = f
#endif
   end subroutine set_m32_f

   !--------------------------------------------------------------------
   elemental pure subroutine set_m32_d(m,d)
      type(fdps_f32mat), intent(inout) :: m
      real(kind=c_double), intent(in) :: d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      m%xx = d
      m%yy = d
      m%zz = d
      m%xy = d
      m%xz = d
      m%yz = d
#else
      m%xx = d
      m%yy = d
      m%xy = d
#endif
   end subroutine set_m32_d

   !--------------------------------------------------------------------
   elemental pure subroutine set_m32_m32(a,b)
      type(fdps_f32mat), intent(inout) :: a
      type(fdps_f32mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      a%xx = b%xx
      a%yy = b%yy
      a%zz = b%zz
      a%xy = b%xy
      a%xz = b%xz
      a%yz = b%yz
#else
      a%xx = b%xx
      a%yy = b%yy
      a%xy = b%xy
#endif
   end subroutine set_m32_m32

   !--------------------------------------------------------------------
   elemental pure subroutine set_m32_m64(a,b)
      type(fdps_f32mat), intent(inout) :: a
      type(fdps_f64mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      a%xx = b%xx
      a%yy = b%yy
      a%zz = b%zz
      a%xy = b%xy
      a%xz = b%xz
      a%yz = b%yz
#else
      a%xx = b%xx
      a%yy = b%yy
      a%xy = b%xy
#endif
   end subroutine set_m32_m64

   !--------------------------------------------------------------------
   elemental pure subroutine set_m64_i(m,i)
      type(fdps_f64mat), intent(inout) :: m
      integer(kind=c_int), intent(in) :: i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      m%xx = real(i,kind=c_double)
      m%yy = real(i,kind=c_double)
      m%zz = real(i,kind=c_double)
      m%xy = real(i,kind=c_double)
      m%xz = real(i,kind=c_double)
      m%yz = real(i,kind=c_double)
#else
      m%xx = real(i,kind=c_double)
      m%yy = real(i,kind=c_double)
      m%xy = real(i,kind=c_double)
#endif
   end subroutine set_m64_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_m64_f(m,f)
      type(fdps_f64mat), intent(inout) :: m
      real(kind=c_float), intent(in) :: f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      m%xx = f
      m%yy = f
      m%zz = f
      m%xy = f
      m%xz = f
      m%yz = f
#else
      m%xx = f
      m%yy = f
      m%xy = f
#endif
   end subroutine set_m64_f

   !--------------------------------------------------------------------
   elemental pure subroutine set_m64_d(m,d)
      type(fdps_f64mat), intent(inout) :: m
      real(kind=c_double), intent(in) :: d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      m%xx = d
      m%yy = d
      m%zz = d
      m%xy = d
      m%xz = d
      m%yz = d
#else
      m%xx = d
      m%yy = d
      m%xy = d
#endif
   end subroutine set_m64_d

   !--------------------------------------------------------------------
   elemental pure subroutine set_m64_m32(a,b)
      type(fdps_f64mat), intent(inout) :: a
      type(fdps_f32mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      a%xx = b%xx
      a%yy = b%yy
      a%zz = b%zz
      a%xy = b%xy
      a%xz = b%xz
      a%yz = b%yz
#else
      a%xx = b%xx
      a%yy = b%yy
      a%xy = b%xy
#endif
   end subroutine set_m64_m32

   !--------------------------------------------------------------------
   elemental pure subroutine set_m64_m64(a,b)
      type(fdps_f64mat), intent(inout) :: a
      type(fdps_f64mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      a%xx = b%xx
      a%yy = b%yy
      a%zz = b%zz
      a%xy = b%xy
      a%xz = b%xz
      a%yz = b%yz
#else
      a%xx = b%xx
      a%yy = b%yy
      a%xy = b%xy
#endif
   end subroutine set_m64_m64

   !####################################################################
   elemental pure type(fdps_f32mat) function add_m32_m32(a,b)
      type(fdps_f32mat), intent(in) :: a,b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      add_m32_m32%xx = a%xx + b%xx
      add_m32_m32%yy = a%yy + b%yy
      add_m32_m32%zz = a%zz + b%zz
      add_m32_m32%xy = a%xy + b%xy
      add_m32_m32%xz = a%xz + b%xz
      add_m32_m32%yz = a%yz + b%yz
#else
      add_m32_m32%xx = a%xx + b%xx
      add_m32_m32%yy = a%yy + b%yy
      add_m32_m32%xy = a%xy + b%xy
#endif
   end function add_m32_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function add_m32_m64(a,b)
      type(fdps_f32mat), intent(in) :: a
      type(fdps_f64mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      add_m32_m64%xx = a%xx + b%xx
      add_m32_m64%yy = a%yy + b%yy
      add_m32_m64%zz = a%zz + b%zz
      add_m32_m64%xy = a%xy + b%xy
      add_m32_m64%xz = a%xz + b%xz
      add_m32_m64%yz = a%yz + b%yz
#else
      add_m32_m64%xx = a%xx + b%xx
      add_m32_m64%yy = a%yy + b%yy
      add_m32_m64%xy = a%xy + b%xy
#endif
   end function add_m32_m64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function add_m64_m32(a,b)
      type(fdps_f64mat), intent(in) :: a
      type(fdps_f32mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      add_m64_m32%xx = a%xx + b%xx
      add_m64_m32%yy = a%yy + b%yy
      add_m64_m32%zz = a%zz + b%zz
      add_m64_m32%xy = a%xy + b%xy
      add_m64_m32%xz = a%xz + b%xz
      add_m64_m32%yz = a%yz + b%yz
#else
      add_m64_m32%xx = a%xx + b%xx
      add_m64_m32%yy = a%yy + b%yy
      add_m64_m32%xy = a%xy + b%xy
#endif
   end function add_m64_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function add_m64_m64(a,b)
      type(fdps_f64mat), intent(in) :: a
      type(fdps_f64mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      add_m64_m64%xx = a%xx + b%xx
      add_m64_m64%yy = a%yy + b%yy
      add_m64_m64%zz = a%zz + b%zz
      add_m64_m64%xy = a%xy + b%xy
      add_m64_m64%xz = a%xz + b%xz
      add_m64_m64%yz = a%yz + b%yz
#else
      add_m64_m64%xx = a%xx + b%xx
      add_m64_m64%yy = a%yy + b%yy
      add_m64_m64%xy = a%xy + b%xy
#endif
   end function add_m64_m64

   !====================================================================
   elemental pure type(fdps_f32mat) function pos_m32(m)
      type(fdps_f32mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      pos_m32%xx = m%xx
      pos_m32%yy = m%yy
      pos_m32%zz = m%zz
      pos_m32%xy = m%xy
      pos_m32%xz = m%xz
      pos_m32%yz = m%yz
#else
      pos_m32%xx = m%xx
      pos_m32%yy = m%yy
      pos_m32%xy = m%xy
#endif
   end function pos_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function pos_m64(m)
      type(fdps_f64mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      pos_m64%xx = m%xx
      pos_m64%yy = m%yy
      pos_m64%zz = m%zz
      pos_m64%xy = m%xy
      pos_m64%xz = m%xz
      pos_m64%yz = m%yz
#else
      pos_m64%xx = m%xx
      pos_m64%yy = m%yy
      pos_m64%xy = m%xy
#endif
   end function pos_m64

   !####################################################################
   elemental pure type(fdps_f32mat) function sub_m32_m32(a,b)
      type(fdps_f32mat), intent(in) :: a,b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      sub_m32_m32%xx = a%xx - b%xx
      sub_m32_m32%yy = a%yy - b%yy
      sub_m32_m32%zz = a%zz - b%zz
      sub_m32_m32%xy = a%xy - b%xy
      sub_m32_m32%xz = a%xz - b%xz
      sub_m32_m32%yz = a%yz - b%yz
#else
      sub_m32_m32%xx = a%xx - b%xx
      sub_m32_m32%yy = a%yy - b%yy
      sub_m32_m32%xy = a%xy - b%xy
#endif
   end function sub_m32_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function sub_m32_m64(a,b)
      type(fdps_f32mat), intent(in) :: a
      type(fdps_f64mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      sub_m32_m64%xx = a%xx - b%xx
      sub_m32_m64%yy = a%yy - b%yy
      sub_m32_m64%zz = a%zz - b%zz
      sub_m32_m64%xy = a%xy - b%xy
      sub_m32_m64%xz = a%xz - b%xz
      sub_m32_m64%yz = a%yz - b%yz
#else
      sub_m32_m64%xx = a%xx - b%xx
      sub_m32_m64%yy = a%yy - b%yy
      sub_m32_m64%xy = a%xy - b%xy
#endif
   end function sub_m32_m64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function sub_m64_m32(a,b)
      type(fdps_f64mat), intent(in) :: a
      type(fdps_f32mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      sub_m64_m32%xx = a%xx - b%xx
      sub_m64_m32%yy = a%yy - b%yy
      sub_m64_m32%zz = a%zz - b%zz
      sub_m64_m32%xy = a%xy - b%xy
      sub_m64_m32%xz = a%xz - b%xz
      sub_m64_m32%yz = a%yz - b%yz
#else
      sub_m64_m32%xx = a%xx - b%xx
      sub_m64_m32%yy = a%yy - b%yy
      sub_m64_m32%xy = a%xy - b%xy
#endif
   end function sub_m64_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function sub_m64_m64(a,b)
      type(fdps_f64mat), intent(in) :: a
      type(fdps_f64mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      sub_m64_m64%xx = a%xx - b%xx
      sub_m64_m64%yy = a%yy - b%yy
      sub_m64_m64%zz = a%zz - b%zz
      sub_m64_m64%xy = a%xy - b%xy
      sub_m64_m64%xz = a%xz - b%xz
      sub_m64_m64%yz = a%yz - b%yz
#else
      sub_m64_m64%xx = a%xx - b%xx
      sub_m64_m64%yy = a%yy - b%yy
      sub_m64_m64%xy = a%xy - b%xy
#endif
   end function sub_m64_m64

   !====================================================================
   elemental pure type(fdps_f32mat) function neg_m32(m)
      type(fdps_f32mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      neg_m32%xx = - m%xx
      neg_m32%yy = - m%yy
      neg_m32%zz = - m%zz
      neg_m32%xy = - m%xy
      neg_m32%xz = - m%xz
      neg_m32%yz = - m%yz
#else
      neg_m32%xx = - m%xx
      neg_m32%yy = - m%yy
      neg_m32%xy = - m%xy
#endif
   end function neg_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function neg_m64(m)
      type(fdps_f64mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      neg_m64%xx = - m%xx
      neg_m64%yy = - m%yy
      neg_m64%zz = - m%zz
      neg_m64%xy = - m%xy
      neg_m64%xz = - m%xz
      neg_m64%yz = - m%yz
#else
      neg_m64%xx = - m%xx
      neg_m64%yy = - m%yy
      neg_m64%xy = - m%xy
#endif
   end function neg_m64

   !####################################################################
   elemental pure type(fdps_f32mat) function mul_m32_i(m,i)
      type(fdps_f32mat), intent(in) :: m
      integer(kind=c_int), intent(in) :: i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m32_i%xx = m%xx * i 
      mul_m32_i%yy = m%yy * i
      mul_m32_i%zz = m%zz * i
      mul_m32_i%xy = m%xy * i
      mul_m32_i%xz = m%xz * i
      mul_m32_i%yz = m%yz * i
#else
      mul_m32_i%xx = m%xx * i 
      mul_m32_i%yy = m%yy * i
      mul_m32_i%xy = m%xy * i
#endif
   end function mul_m32_i

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32mat) function mul_m32_f(m,f)
      type(fdps_f32mat), intent(in) :: m
      real(kind=c_float), intent(in) :: f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m32_f%xx = m%xx * f
      mul_m32_f%yy = m%yy * f
      mul_m32_f%zz = m%zz * f
      mul_m32_f%xy = m%xy * f
      mul_m32_f%xz = m%xz * f
      mul_m32_f%yz = m%yz * f
#else
      mul_m32_f%xx = m%xx * f
      mul_m32_f%yy = m%yy * f
      mul_m32_f%xy = m%xy * f
#endif
   end function mul_m32_f

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32mat) function mul_m32_d(m,d)
      type(fdps_f32mat), intent(in) :: m
      real(kind=c_double), intent(in) :: d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m32_d%xx = m%xx * d
      mul_m32_d%yy = m%yy * d
      mul_m32_d%zz = m%zz * d
      mul_m32_d%xy = m%xy * d
      mul_m32_d%xz = m%xz * d
      mul_m32_d%yz = m%yz * d
#else
      mul_m32_d%xx = m%xx * d
      mul_m32_d%yy = m%yy * d
      mul_m32_d%xy = m%xy * d
#endif
   end function mul_m32_d

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32mat) function mul_i_m32(i,m)
      integer(kind=c_int), intent(in) :: i
      type(fdps_f32mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_i_m32%xx = i * m%xx
      mul_i_m32%yy = i * m%yy
      mul_i_m32%zz = i * m%zz
      mul_i_m32%xy = i * m%xy
      mul_i_m32%xz = i * m%xz
      mul_i_m32%yz = i * m%yz
#else
      mul_i_m32%xx = i * m%xx
      mul_i_m32%yy = i * m%yy
      mul_i_m32%xy = i * m%xy
#endif
   end function mul_i_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32mat) function mul_f_m32(f,m)
      real(kind=c_float), intent(in) :: f
      type(fdps_f32mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_f_m32%xx = f * m%xx
      mul_f_m32%yy = f * m%yy
      mul_f_m32%zz = f * m%zz
      mul_f_m32%xy = f * m%xy
      mul_f_m32%xz = f * m%xz
      mul_f_m32%yz = f * m%yz
#else
      mul_f_m32%xx = f * m%xx
      mul_f_m32%yy = f * m%yy
      mul_f_m32%xy = f * m%xy
#endif
   end function mul_f_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32mat) function mul_d_m32(d,m)
      real(kind=c_double), intent(in) :: d
      type(fdps_f32mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_d_m32%xx = d * m%xx
      mul_d_m32%yy = d * m%yy
      mul_d_m32%zz = d * m%zz
      mul_d_m32%xy = d * m%xy
      mul_d_m32%xz = d * m%xz
      mul_d_m32%yz = d * m%yz
#else
      mul_d_m32%xx = d * m%xx
      mul_d_m32%yy = d * m%yy
      mul_d_m32%xy = d * m%xy
#endif
   end function mul_d_m32

   !####################################################################
   elemental pure type(fdps_f64mat) function mul_m64_i(m,i)
      type(fdps_f64mat), intent(in) :: m
      integer(kind=c_int), intent(in) :: i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m64_i%xx = m%xx * i 
      mul_m64_i%yy = m%yy * i
      mul_m64_i%zz = m%zz * i
      mul_m64_i%xy = m%xy * i
      mul_m64_i%xz = m%xz * i
      mul_m64_i%yz = m%yz * i
#else
      mul_m64_i%xx = m%xx * i 
      mul_m64_i%yy = m%yy * i
      mul_m64_i%xy = m%xy * i
#endif
   end function mul_m64_i

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function mul_m64_f(m,f)
      type(fdps_f64mat), intent(in) :: m
      real(kind=c_float), intent(in) :: f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m64_f%xx = m%xx * f
      mul_m64_f%yy = m%yy * f
      mul_m64_f%zz = m%zz * f
      mul_m64_f%xy = m%xy * f
      mul_m64_f%xz = m%xz * f
      mul_m64_f%yz = m%yz * f
#else
      mul_m64_f%xx = m%xx * f
      mul_m64_f%yy = m%yy * f
      mul_m64_f%xy = m%xy * f
#endif
   end function mul_m64_f

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function mul_m64_d(m,d)
      type(fdps_f64mat), intent(in) :: m
      real(kind=c_double), intent(in) :: d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m64_d%xx = m%xx * d
      mul_m64_d%yy = m%yy * d
      mul_m64_d%zz = m%zz * d
      mul_m64_d%xy = m%xy * d
      mul_m64_d%xz = m%xz * d
      mul_m64_d%yz = m%yz * d
#else
      mul_m64_d%xx = m%xx * d
      mul_m64_d%yy = m%yy * d
      mul_m64_d%xy = m%xy * d
#endif
   end function mul_m64_d

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function mul_i_m64(i,m)
      integer(kind=c_int), intent(in) :: i
      type(fdps_f64mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_i_m64%xx = i * m%xx
      mul_i_m64%yy = i * m%yy
      mul_i_m64%zz = i * m%zz
      mul_i_m64%xy = i * m%xy
      mul_i_m64%xz = i * m%xz
      mul_i_m64%yz = i * m%yz
#else
      mul_i_m64%xx = i * m%xx
      mul_i_m64%yy = i * m%yy
      mul_i_m64%xy = i * m%xy
#endif
   end function mul_i_m64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function mul_f_m64(f,m)
      real(kind=c_float), intent(in) :: f
      type(fdps_f64mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_f_m64%xx = f * m%xx
      mul_f_m64%yy = f * m%yy
      mul_f_m64%zz = f * m%zz
      mul_f_m64%xy = f * m%xy
      mul_f_m64%xz = f * m%xz
      mul_f_m64%yz = f * m%yz
#else
      mul_f_m64%xx = f * m%xx
      mul_f_m64%yy = f * m%yy
      mul_f_m64%xy = f * m%xy
#endif
   end function mul_f_m64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function mul_d_m64(d,m)
      real(kind=c_double), intent(in) :: d
      type(fdps_f64mat), intent(in) :: m
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_d_m64%xx = d * m%xx
      mul_d_m64%yy = d * m%yy
      mul_d_m64%zz = d * m%zz
      mul_d_m64%xy = d * m%xy
      mul_d_m64%xz = d * m%xz
      mul_d_m64%yz = d * m%yz
#else
      mul_d_m64%xx = d * m%xx
      mul_d_m64%yy = d * m%yy
      mul_d_m64%xy = d * m%xy
#endif
   end function mul_d_m64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32mat) function mul_m32_m32(a,b)
      type(fdps_f32mat), intent(in) :: a,b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m32_m32%xx = a%xx * b%xx + a%xy * b%xy + a%xz * b%xz
      mul_m32_m32%yy = a%xy * b%xy + a%yy * b%yy + a%yz * b%yz
      mul_m32_m32%zz = a%xz * b%xz + a%yz * b%yz + a%zz * b%zz
      mul_m32_m32%xy = a%xx * b%xy + a%xy * b%yy + a%xz * b%yz
      mul_m32_m32%xz = a%xx * b%xz + a%xy * b%yz + a%xz * b%zz
      mul_m32_m32%yz = a%xy * b%xz + a%yy * b%yz + a%yz * b%zz
#else
      mul_m32_m32%xx = a%xx * b%xx + a%xy * b%xy
      mul_m32_m32%yy = a%xy * b%xy + a%yy * b%yy
      mul_m32_m32%xy = a%xx * b%xy + a%xy * b%yy
#endif
   end function mul_m32_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function mul_m64_m32(a,b)
      type(fdps_f64mat), intent(in) :: a
      type(fdps_f32mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m64_m32%xx = a%xx * b%xx + a%xy * b%xy + a%xz * b%xz
      mul_m64_m32%yy = a%xy * b%xy + a%yy * b%yy + a%yz * b%yz
      mul_m64_m32%zz = a%xz * b%xz + a%yz * b%yz + a%zz * b%zz
      mul_m64_m32%xy = a%xx * b%xy + a%xy * b%yy + a%xz * b%yz
      mul_m64_m32%xz = a%xx * b%xz + a%xy * b%yz + a%xz * b%zz
      mul_m64_m32%yz = a%xy * b%xz + a%yy * b%yz + a%yz * b%zz
#else
      mul_m64_m32%xx = a%xx * b%xx + a%xy * b%xy
      mul_m64_m32%yy = a%xy * b%xy + a%yy * b%yy
      mul_m64_m32%xy = a%xx * b%xy + a%xy * b%yy
#endif
   end function mul_m64_m32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function mul_m32_m64(a,b)
      type(fdps_f32mat), intent(in) :: a
      type(fdps_f64mat), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m32_m64%xx = a%xx * b%xx + a%xy * b%xy + a%xz * b%xz
      mul_m32_m64%yy = a%xy * b%xy + a%yy * b%yy + a%yz * b%yz
      mul_m32_m64%zz = a%xz * b%xz + a%yz * b%yz + a%zz * b%zz
      mul_m32_m64%xy = a%xx * b%xy + a%xy * b%yy + a%xz * b%yz
      mul_m32_m64%xz = a%xx * b%xz + a%xy * b%yz + a%xz * b%zz
      mul_m32_m64%yz = a%xy * b%xz + a%yy * b%yz + a%yz * b%zz
#else
      mul_m32_m64%xx = a%xx * b%xx + a%xy * b%xy
      mul_m32_m64%yy = a%xy * b%xy + a%yy * b%yy
      mul_m32_m64%xy = a%xx * b%xy + a%xy * b%yy
#endif
   end function mul_m32_m64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function mul_m64_m64(a,b)
      type(fdps_f64mat), intent(in) :: a,b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      mul_m64_m64%xx = a%xx * b%xx + a%xy * b%xy + a%xz * b%xz
      mul_m64_m64%yy = a%xy * b%xy + a%yy * b%yy + a%yz * b%yz
      mul_m64_m64%zz = a%xz * b%xz + a%yz * b%yz + a%zz * b%zz
      mul_m64_m64%xy = a%xx * b%xy + a%xy * b%yy + a%xz * b%yz
      mul_m64_m64%xz = a%xx * b%xz + a%xy * b%yz + a%xz * b%zz
      mul_m64_m64%yz = a%xy * b%xz + a%yy * b%yz + a%yz * b%zz
#else
      mul_m64_m64%xx = a%xx * b%xx + a%xy * b%xy
      mul_m64_m64%yy = a%xy * b%xy + a%yy * b%yy
      mul_m64_m64%xy = a%xx * b%xy + a%xy * b%yy
#endif
   end function mul_m64_m64

   !####################################################################
   elemental pure type(fdps_f32mat) function div_m32_i(m,i)
      type(fdps_f32mat), intent(in) :: m
      integer(kind=c_int), intent(in) :: i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      div_m32_i%xx = m%xx / i
      div_m32_i%yy = m%yy / i
      div_m32_i%zz = m%zz / i
      div_m32_i%xy = m%xy / i
      div_m32_i%xz = m%xz / i
      div_m32_i%yz = m%yz / i
#else
      div_m32_i%xx = m%xx / i
      div_m32_i%yy = m%yy / i
      div_m32_i%xy = m%xy / i
#endif
   end function div_m32_i

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32mat) function div_m32_f(m,f)
      type(fdps_f32mat), intent(in) :: m
      real(kind=c_float), intent(in) :: f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      div_m32_f%xx = m%xx / f
      div_m32_f%yy = m%yy / f
      div_m32_f%zz = m%zz / f
      div_m32_f%xy = m%xy / f
      div_m32_f%xz = m%xz / f
      div_m32_f%yz = m%yz / f
#else
      div_m32_f%xx = m%xx / f
      div_m32_f%yy = m%yy / f
      div_m32_f%xy = m%xy / f
#endif
   end function div_m32_f

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32mat) function div_m32_d(m,d)
      type(fdps_f32mat), intent(in) :: m
      real(kind=c_double), intent(in) :: d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      div_m32_d%xx = m%xx / d
      div_m32_d%yy = m%yy / d
      div_m32_d%zz = m%zz / d
      div_m32_d%xy = m%xy / d
      div_m32_d%xz = m%xz / d
      div_m32_d%yz = m%yz / d
#else
      div_m32_d%xx = m%xx / d
      div_m32_d%yy = m%yy / d
      div_m32_d%xy = m%xy / d
#endif
   end function div_m32_d

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function div_m64_i(m,i)
      type(fdps_f64mat), intent(in) :: m
      integer(kind=c_int), intent(in) :: i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      div_m64_i%xx = m%xx / i
      div_m64_i%yy = m%yy / i
      div_m64_i%zz = m%zz / i
      div_m64_i%xy = m%xy / i
      div_m64_i%xz = m%xz / i
      div_m64_i%yz = m%yz / i
#else
      div_m64_i%xx = m%xx / i
      div_m64_i%yy = m%yy / i
      div_m64_i%xy = m%xy / i
#endif
   end function div_m64_i

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function div_m64_f(m,f)
      type(fdps_f64mat), intent(in) :: m
      real(kind=c_float), intent(in) :: f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      div_m64_f%xx = m%xx / f
      div_m64_f%yy = m%yy / f
      div_m64_f%zz = m%zz / f
      div_m64_f%xy = m%xy / f
      div_m64_f%xz = m%xz / f
      div_m64_f%yz = m%yz / f
#else
      div_m64_f%xx = m%xx / f
      div_m64_f%yy = m%yy / f
      div_m64_f%xy = m%xy / f
#endif
   end function div_m64_f

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64mat) function div_m64_d(m,d)
      type(fdps_f64mat), intent(in) :: m
      real(kind=c_double), intent(in) :: d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      div_m64_d%xx = m%xx / d
      div_m64_d%yy = m%yy / d
      div_m64_d%zz = m%zz / d
      div_m64_d%xy = m%xy / d
      div_m64_d%xz = m%xz / d
      div_m64_d%yz = m%yz / d
#else
      div_m64_d%xx = m%xx / d
      div_m64_d%yy = m%yy / d
      div_m64_d%xy = m%xy / d
#endif
   end function div_m64_d

end module fdps_matrix

