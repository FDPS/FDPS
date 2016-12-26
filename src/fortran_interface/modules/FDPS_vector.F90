!==============================
!   MODULE: FDPS vector
!==============================
module fdps_vector
   use, intrinsic :: iso_c_binding
   implicit none

   !**** PS::F32vec
   type, public, bind(c) :: fdps_f32vec
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION 
      real(kind=c_float) :: x,y
#else
      real(kind=c_float) :: x,y,z
#endif
   end type fdps_f32vec

   !**** PS::F64vec
   type, public, bind(c) :: fdps_f64vec
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION 
      real(kind=c_double) :: x,y
#else
      real(kind=c_double) :: x,y,z
#endif
   end type fdps_f64vec

   !* Private routines
   private :: set_v32_i   ! i := c_int
   private :: set_v32_f   ! f := c_float
   private :: set_v32_fa  ! a := array
   private :: set_v32_d   ! d := c_double
   private :: set_v32_da
   private :: set_v32_v32
   private :: set_v32_v64
   private :: set_v64_i
   private :: set_v64_f
   private :: set_v64_fa
   private :: set_v64_d
   private :: set_v64_da
   private :: set_v64_v32
   private :: set_v64_v64

   private :: add_v32_fa
   private :: add_fa_v32
   private :: add_v32_v32
   private :: add_v32_da
   private :: add_da_v32
   private :: add_v64_da
   private :: add_da_v64
   private :: add_v32_v64
   private :: add_v64_v32
   private :: add_v64_v64
   private :: pos_v32
   private :: pos_v64

   private :: sub_v32_fa
   private :: sub_fa_v32
   private :: sub_v32_v32
   private :: sub_v32_da
   private :: sub_da_v32
   private :: sub_v64_da
   private :: sub_da_v64
   private :: sub_v32_v64
   private :: sub_v64_v32
   private :: sub_v64_v64
   private :: neg_v32
   private :: neg_v64

   private :: mul_v32_i
   private :: mul_v32_f
   private :: mul_v32_d
   private :: mul_i_v32
   private :: mul_f_v32
   private :: mul_d_v32
   private :: mul_v64_i
   private :: mul_v64_f
   private :: mul_v64_d
   private :: mul_i_v64
   private :: mul_f_v64
   private :: mul_d_v64
   private :: mul_v32_fa
   private :: mul_fa_v32
   private :: mul_v32_v32
   private :: mul_v32_da
   private :: mul_da_v32
   private :: mul_v64_da
   private :: mul_da_v64
   private :: mul_v64_v32
   private :: mul_v32_v64
   private :: mul_v64_v64

   private :: div_v32_i
   private :: div_v32_f
   private :: div_v32_d

   private :: div_v64_i
   private :: div_v64_f
   private :: div_v64_d

   !* Declare interface operators
   interface assignment (=)
      module procedure set_v32_i
      module procedure set_v32_f
      module procedure set_v32_fa
      module procedure set_v32_d
      module procedure set_v32_da
      module procedure set_v32_v32
      module procedure set_v32_v64

      module procedure set_v64_i
      module procedure set_v64_f
      module procedure set_v64_fa
      module procedure set_v64_d
      module procedure set_v64_da
      module procedure set_v64_v32
      module procedure set_v64_v64
   end interface 

   interface operator (+)
      module procedure add_v32_fa
      module procedure add_fa_v32
      module procedure add_v32_v32
      module procedure add_v32_da
      module procedure add_da_v32
      module procedure add_v64_da
      module procedure add_da_v64
      module procedure add_v32_v64
      module procedure add_v64_v32
      module procedure add_v64_v64
      module procedure pos_v32
      module procedure pos_v64
   end interface 

   interface operator (-)
      module procedure sub_v32_fa
      module procedure sub_fa_v32
      module procedure sub_v32_v32
      module procedure sub_v32_da
      module procedure sub_da_v32
      module procedure sub_v64_da
      module procedure sub_da_v64
      module procedure sub_v32_v64
      module procedure sub_v64_v32
      module procedure sub_v64_v64
      module procedure neg_v32
      module procedure neg_v64
   end interface

   interface operator (*)
      module procedure mul_v32_i
      module procedure mul_v32_f
      module procedure mul_v32_d
      module procedure mul_i_v32
      module procedure mul_f_v32
      module procedure mul_d_v32
      module procedure mul_v64_i
      module procedure mul_v64_f
      module procedure mul_v64_d
      module procedure mul_i_v64
      module procedure mul_f_v64
      module procedure mul_d_v64

      module procedure mul_v32_fa
      module procedure mul_fa_v32
      module procedure mul_v32_v32
      module procedure mul_v32_da
      module procedure mul_da_v32
      module procedure mul_v64_da
      module procedure mul_da_v64
      module procedure mul_v64_v32
      module procedure mul_v32_v64
      module procedure mul_v64_v64
   end interface

   interface operator (/)
      module procedure div_v32_i
      module procedure div_v32_f
      module procedure div_v32_d

      module procedure div_v64_i
      module procedure div_v64_f
      module procedure div_v64_d
   end interface

   contains

   !--------------------------------------------------------------------
   elemental pure subroutine set_v32_i(v,i)
      type(fdps_f32vec), intent(out) :: v
      integer(kind=c_int), intent(in) :: i
      v%x = real(i,kind=c_float)
      v%y = real(i,kind=c_float)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = real(i,kind=c_float)
#endif
   end subroutine set_v32_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_v32_f(v,f)
      type(fdps_f32vec), intent(out) :: v
      real(kind=c_float), intent(in) :: f
      v%x = f
      v%y = f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = f
#endif
   end subroutine set_v32_f

   !--------------------------------------------------------------------
   pure subroutine set_v32_fa(v,f)
      type(fdps_f32vec), intent(out) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      v%x = f(1)
      v%y = f(2)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = f(3)
#endif
   end subroutine set_v32_fa

   !--------------------------------------------------------------------
   elemental pure subroutine set_v32_d(v,d)
      type(fdps_f32vec), intent(out) :: v
      real(kind=c_double), intent(in) :: d
      v%x = real(d,kind=c_float)
      v%y = real(d,kind=c_float)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = real(d,kind=c_float)
#endif
   end subroutine set_v32_d

   !--------------------------------------------------------------------
   pure subroutine set_v32_da(v,d)
      type(fdps_f32vec), intent(out) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      v%x = real(d(1),kind=c_float)
      v%y = real(d(2),kind=c_float)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = real(d(3),kind=c_float)
#endif
   end subroutine set_v32_da

   !--------------------------------------------------------------------
   elemental pure subroutine set_v32_v32(a,b)
      type(fdps_f32vec), intent(out) :: a
      type(fdps_f32vec), intent(in) :: b
      a%x = b%x
      a%y = b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      a%z = b%z
#endif
   end subroutine set_v32_v32

   !--------------------------------------------------------------------
   elemental pure subroutine set_v32_v64(a,b)
      type(fdps_f32vec), intent(out) :: a
      type(fdps_f64vec), intent(in) :: b
      a%x = real(b%x,kind=c_float)
      a%y = real(b%y,kind=c_float)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      a%z = real(b%z,kind=c_float)
#endif
   end subroutine set_v32_v64

   !--------------------------------------------------------------------
   elemental pure subroutine set_v64_i(v,i)
      type(fdps_f64vec), intent(out) :: v
      integer(kind=c_int), intent(in) :: i
      v%x = real(i,kind=c_double)
      v%y = real(i,kind=c_double)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = real(i,kind=c_double)
#endif
   end subroutine set_v64_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_v64_f(v,f)
      type(fdps_f64vec), intent(out) :: v
      real(kind=c_float), intent(in) :: f
      v%x = real(f,kind=c_double)
      v%y = real(f,kind=c_double)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = real(f,kind=c_double)
#endif
   end subroutine set_v64_f

   !--------------------------------------------------------------------
   pure subroutine set_v64_fa(v,f)
      type(fdps_f64vec), intent(out) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      v%x = real(f(1),kind=c_double)
      v%y = real(f(2),kind=c_double)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = real(f(3),kind=c_double)
#endif
   end subroutine set_v64_fa

   !--------------------------------------------------------------------
   elemental pure subroutine set_v64_d(v,d)
      type(fdps_f64vec), intent(out) :: v
      real(kind=c_double), intent(in) :: d
      v%x = d
      v%y = d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = d
#endif
   end subroutine set_v64_d

   !--------------------------------------------------------------------
   pure subroutine set_v64_da(v,d)
      type(fdps_f64vec), intent(out) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      v%x = d(1)
      v%y = d(2)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      v%z = d(3)
#endif
   end subroutine set_v64_da

   !--------------------------------------------------------------------
   elemental pure subroutine set_v64_v32(a,b)
      type(fdps_f64vec), intent(out) :: a
      type(fdps_f32vec), intent(in) :: b
      a%x = real(b%x,kind=c_double)
      a%y = real(b%y,kind=c_double)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      a%z = real(b%z,kind=c_double)
#endif
   end subroutine set_v64_v32

   !--------------------------------------------------------------------
   elemental pure subroutine set_v64_v64(a,b)
      type(fdps_f64vec), intent(out) :: a
      type(fdps_f64vec), intent(in) :: b
      a%x = b%x
      a%y = b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      a%z = b%z
#endif
   end subroutine set_v64_v64

   !####################################################################
   pure type(fdps_f32vec) function add_v32_fa(v,f)
     type(fdps_f32vec), intent(in) :: v
     real(kind=c_float), dimension(:), intent(IN) :: f
     add_v32_fa%x = v%x + f(1)
     add_v32_fa%y = v%y + f(2)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
     add_v32_fa%z = v%z + f(3)
#endif
   end function add_v32_fa

   !--------------------------------------------------------------------
   pure type(fdps_f32vec) function add_fa_v32(f,v)
     real(kind=c_float), dimension(:), intent(IN) :: f
     type(fdps_f32vec), intent(in) :: v
     add_fa_v32%x = f(1) + v%x
     add_fa_v32%y = f(2) + v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
     add_fa_v32%z = f(3) + v%z
#endif
   end function add_fa_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function add_v32_v32(a,b)
     type(fdps_f32vec), intent(in) :: a,b
     add_v32_v32%x = a%x + b%x
     add_v32_v32%y = a%y + b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
     add_v32_v32%z = a%z + b%z
#endif
   end function add_v32_v32

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function add_v32_da(v,d)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      add_v32_da%x = v%x + d(1)
      add_v32_da%y = v%y + d(2)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      add_v32_da%z = v%z + d(3)
#endif
   end function add_v32_da

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function add_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(fdps_f32vec), intent(in) :: v
      add_da_v32%x = d(1) + v%x 
      add_da_v32%y = d(2) + v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      add_da_v32%z = d(3) + v%z
#endif
   end function add_da_v32

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function add_v64_da(v,d)
      type(fdps_f64vec), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      add_v64_da%x = v%x + d(1)
      add_v64_da%y = v%y + d(2)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      add_v64_da%z = v%z + d(3)
#endif
   end function add_v64_da

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function add_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(fdps_f64vec), intent(in) :: v
      add_da_v64%x = d(1) + v%x 
      add_da_v64%y = d(2) + v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      add_da_v64%z = d(3) + v%z
#endif
   end function add_da_v64

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function add_v32_v64(a,b)
      type(fdps_f32vec), intent(in) :: a
      type(fdps_f64vec), intent(in) :: b
      add_v32_v64%x = a%x + b%x
      add_v32_v64%y = a%y + b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      add_v32_v64%z = a%z + b%z
#endif
   end function add_v32_v64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function add_v64_v32(a,b)
      type(fdps_f64vec), intent(in) :: a
      type(fdps_f32vec), intent(in) :: b
      add_v64_v32%x = a%x + b%x
      add_v64_v32%y = a%y + b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      add_v64_v32%z = a%z + b%z
#endif
   end function add_v64_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function add_v64_v64(a,b)
      type(fdps_f64vec), intent(in) :: a,b
      add_v64_v64%x = a%x + b%x
      add_v64_v64%y = a%y + b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      add_v64_v64%z = a%z + b%z
#endif
   end function add_v64_v64

   !====================================================================
   elemental pure type(fdps_f32vec) function pos_v32(v)
      type(fdps_f32vec), intent(in) :: v
      pos_v32%x = v%x
      pos_v32%y = v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      pos_v32%z = v%z
#endif
   end function pos_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function pos_v64(v)
      type(fdps_f64vec), intent(in) :: v
      pos_v64%x = v%x
      pos_v64%y = v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      pos_v64%z = v%z
#endif
   end function pos_v64

   !####################################################################
   pure type(fdps_f32vec) function sub_v32_fa(v,f)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      sub_v32_fa%x = v%x - f(1)
      sub_v32_fa%y = v%y - f(2)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_v32_fa%z = v%z - f(3)
#endif
   end function sub_v32_fa

   !--------------------------------------------------------------------
   pure type(fdps_f32vec) function sub_fa_v32(f,v)
      real(kind=c_float), dimension(:), intent(in) :: f
      type(fdps_f32vec), intent(in) :: v
      sub_fa_v32%x = f(1) - v%x
      sub_fa_v32%y = f(2) - v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_fa_v32%z = f(3) - v%z
#endif
   end function sub_fa_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function sub_v32_v32(a,b)
      type(fdps_f32vec), intent(in) :: a,b
      sub_v32_v32%x = a%x - b%x
      sub_v32_v32%y = a%y - b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_v32_v32%z = a%z - b%z
#endif
   end function sub_v32_v32

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function sub_v32_da(v,d)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      sub_v32_da%x = v%x - d(1)
      sub_v32_da%y = v%y - d(2)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_v32_da%z = v%z - d(3)
#endif
   end function sub_v32_da

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function sub_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(fdps_f32vec), intent(in) :: v
      sub_da_v32%x = d(1) - v%x
      sub_da_v32%y = d(2) - v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_da_v32%z = d(3) - v%z
#endif
   end function sub_da_v32

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function sub_v64_da(v,d)
      type(fdps_f64vec), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      sub_v64_da%x = v%x - d(1)
      sub_v64_da%y = v%y - d(2)
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_v64_da%z = v%z - d(3)
#endif
   end function sub_v64_da

   !--------------------------------------------------------------------
   pure type(fdps_f64vec) function sub_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(fdps_f64vec), intent(in) :: v
      sub_da_v64%x = d(1) - v%x
      sub_da_v64%y = d(2) - v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_da_v64%z = d(3) - v%z
#endif
   end function sub_da_v64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function sub_v32_v64(a,b)
      type(fdps_f32vec), intent(in) :: a
      type(fdps_f64vec), intent(in) :: b
      sub_v32_v64%x = a%x - b%x
      sub_v32_v64%y = a%y - b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_v32_v64%z = a%z - b%z
#endif
   end function sub_v32_v64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function sub_v64_v32(a,b)
      type(fdps_f64vec), intent(in) :: a
      type(fdps_f32vec), intent(in) :: b
      sub_v64_v32%x = a%x - b%x
      sub_v64_v32%y = a%y - b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_v64_v32%z = a%z - b%z
#endif
   end function sub_v64_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function sub_v64_v64(a,b)
      type(fdps_f64vec), intent(in) :: a,b
      sub_v64_v64%x = a%x - b%x
      sub_v64_v64%y = a%y - b%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      sub_v64_v64%z = a%z - b%z
#endif
   end function sub_v64_v64

   !====================================================================
   elemental pure type(fdps_f32vec) function neg_v32(v)
      type(fdps_f32vec), intent(in) :: v
      neg_v32%x = - v%x
      neg_v32%y = - v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      neg_v32%z = - v%z
#endif
   end function neg_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function neg_v64(v)
      type(fdps_f64vec), intent(in) :: v
      neg_v64%x = - v%x
      neg_v64%y = - v%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      neg_v64%z = - v%z
#endif
   end function neg_v64

   !####################################################################
   elemental pure type(fdps_f32vec) function mul_v32_i(v,i)
      type(fdps_f32vec), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      mul_v32_i%x = v%x * i
      mul_v32_i%y = v%y * i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v32_i%z = v%z * i
#endif
   end function mul_v32_i

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function mul_v32_f(v,f)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      mul_v32_f%x = v%x * f
      mul_v32_f%y = v%y * f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v32_f%z = v%z * f
#endif
   end function mul_v32_f

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function mul_v32_d(v,d)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      mul_v32_d%x = v%x * d
      mul_v32_d%y = v%y * d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v32_d%z = v%z * d
#endif
   end function mul_v32_d

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function mul_i_v32(i,v)
      integer(kind=c_int), intent(in) :: i
      type(fdps_f32vec), intent(in) :: v
      mul_i_v32%x = i * v%x 
      mul_i_v32%y = i * v%y 
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_i_v32%z = i * v%z 
#endif
   end function mul_i_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function mul_f_v32(f,v)
      real(kind=c_float), intent(in) :: f
      type(fdps_f32vec), intent(in) :: v
      mul_f_v32%x = f * v%x 
      mul_f_v32%y = f * v%y 
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_f_v32%z = f * v%z 
#endif
   end function mul_f_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function mul_d_v32(d,v)
      real(kind=c_double), intent(in) :: d
      type(fdps_f32vec), intent(in) :: v
      mul_d_v32%x = d * v%x 
      mul_d_v32%y = d * v%y 
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_d_v32%z = d * v%z 
#endif
   end function mul_d_v32

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function mul_v64_i(v,i)
      type(fdps_f64vec), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      mul_v64_i%x = v%x * i
      mul_v64_i%y = v%y * i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v64_i%z = v%z * i
#endif
   end function mul_v64_i

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function mul_v64_f(v,f)
      type(fdps_f64vec), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      mul_v64_f%x = v%x * f
      mul_v64_f%y = v%y * f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v64_f%z = v%z * f
#endif
   end function mul_v64_f

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function mul_v64_d(v,d)
      type(fdps_f64vec), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      mul_v64_d%x = v%x * d
      mul_v64_d%y = v%y * d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v64_d%z = v%z * d
#endif
   end function mul_v64_d

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function mul_i_v64(i,v)
      integer(kind=c_int), intent(in) :: i
      type(fdps_f64vec), intent(in) :: v
      mul_i_v64%x = i * v%x 
      mul_i_v64%y = i * v%y 
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_i_v64%z = i * v%z 
#endif
   end function mul_i_v64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function mul_f_v64(f,v)
      real(kind=c_float), intent(in) :: f
      type(fdps_f64vec), intent(in) :: v
      mul_f_v64%x = f * v%x 
      mul_f_v64%y = f * v%y 
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_f_v64%z = f * v%z 
#endif
   end function mul_f_v64

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function mul_d_v64(d,v)
      real(kind=c_double), intent(in) :: d
      type(fdps_f64vec), intent(in) :: v
      mul_d_v64%x = d * v%x 
      mul_d_v64%y = d * v%y 
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_d_v64%z = d * v%z 
#endif
   end function mul_d_v64

   !--------------------------------------------------------------------
   pure real(kind=c_float) function mul_v32_fa(v,f)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v32_fa = v%x * f(1) + v%y * f(2) + v%z * f(3)
#else
      mul_v32_fa = v%x * f(1) + v%y * f(2)
#endif
   end function mul_v32_fa

   !--------------------------------------------------------------------
   pure real(kind=c_float) function mul_fa_v32(f,v)
      real(kind=c_float), dimension(:), intent(in) :: f
      type(fdps_f32vec), intent(in) :: v
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_fa_v32 = f(1) * v%x  + f(2) * v%y + f(3) * v%z
#else
      mul_fa_v32 = f(1) * v%x + f(2) * v%y
#endif
   end function mul_fa_v32

   !--------------------------------------------------------------------
   elemental pure real(kind=c_float) function mul_v32_v32(a,b)
      type(fdps_f32vec), intent(in) :: a,b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v32_v32 = a%x * b%x + a%y * b%y + a%z * b%z 
#else
      mul_v32_v32 = a%x * b%x + a%y * b%y
#endif
   end function mul_v32_v32

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_v32_da(v,d)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v32_da = v%x * d(1) + v%y * d(2) + v%z * d(3)
#else
      mul_v32_da = v%x * d(1) + v%y * d(2)
#endif
   end function mul_v32_da

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(fdps_f32vec), intent(in) :: v
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_da_v32 = d(1) * v%x + d(2) * v%y + d(3) * v%z
#else
      mul_da_v32 = d(1) * v%x + d(2) * v%y
#endif
   end function mul_da_v32

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_v64_da(v,d)
      type(fdps_f64vec), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v64_da = v%x * d(1) + v%y * d(2) + v%z * d(3)
#else
      mul_v64_da = v%x * d(1) + v%y * d(2)
#endif
   end function mul_v64_da

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(fdps_f64vec), intent(in) :: v
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_da_v64 = d(1) * v%x + d(2) * v%y + d(3) * v%z
#else
      mul_da_v64 = d(1) * v%x + d(2) * v%y
#endif
   end function mul_da_v64

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_v64_v32(a,b)
      type(fdps_f64vec), intent(in) :: a
      type(fdps_f32vec), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v64_v32 = a%x * b%x + a%y * b%y + a%z * b%z 
#else
      mul_v64_v32 = a%x * b%x + a%y * b%y 
#endif
   end function mul_v64_v32

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_v32_v64(a,b)
      type(fdps_f32vec), intent(in) :: a
      type(fdps_f64vec), intent(in) :: b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v32_v64 = a%x * b%x + a%y * b%y + a%z * b%z 
#else
      mul_v32_v64 = a%x * b%x + a%y * b%y
#endif
   end function mul_v32_v64

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_v64_v64(a,b)
      type(fdps_f64vec), intent(in) :: a,b
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      mul_v64_v64 = a%x * b%x + a%y * b%y + a%z * b%z 
#else
      mul_v64_v64 = a%x * b%x + a%y * b%y 
#endif
   end function mul_v64_v64

   !####################################################################
   elemental pure type(fdps_f32vec) function div_v32_i(v,i)
      type(fdps_f32vec), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      div_v32_i%x = v%x / i
      div_v32_i%y = v%y / i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      div_v32_i%z = v%z / i
#endif
   end function div_v32_i

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function div_v32_f(v,f)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      div_v32_f%x = v%x / f
      div_v32_f%y = v%y / f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      div_v32_f%z = v%z / f
#endif
   end function div_v32_f

   !--------------------------------------------------------------------
   elemental pure type(fdps_f32vec) function div_v32_d(v,d)
      type(fdps_f32vec), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      div_v32_d%x = v%x / d
      div_v32_d%y = v%y / d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      div_v32_d%z = v%z / d
#endif
   end function div_v32_d

   !====================================================================
   elemental pure type(fdps_f64vec) function div_v64_i(v,i)
      type(fdps_f64vec), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      div_v64_i%x = v%x / i
      div_v64_i%y = v%y / i
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      div_v64_i%z = v%z / i
#endif
   end function div_v64_i

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function div_v64_f(v,f)
      type(fdps_f64vec), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      div_v64_f%x = v%x / f
      div_v64_f%y = v%y / f
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      div_v64_f%z = v%z / f
#endif
   end function div_v64_f

   !--------------------------------------------------------------------
   elemental pure type(fdps_f64vec) function div_v64_d(v,d)
      type(fdps_f64vec), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      div_v64_d%x = v%x / d
      div_v64_d%y = v%y / d
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION 
      div_v64_d%z = v%z / d
#endif
   end function div_v64_d

end module fdps_vector
