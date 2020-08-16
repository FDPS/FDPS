!==============================
!   MODULE: PIKG vector
!==============================
module pikg_vector
   use, intrinsic :: iso_c_binding
   implicit none

   !**** PIKG::F32vec2
   type, public, bind(c) :: pikg_f32vec2
      real(kind=c_float) :: x,y
   end type pikg_f32vec2

   !**** PIKG::F32vec3
   type, public, bind(c) :: pikg_f32vec3
      real(kind=c_float) :: x,y,z
   end type pikg_f32vec3

   !**** PIKG::F32vec4
   type, public, bind(c) :: pikg_f32vec4
      real(kind=c_float) :: x,y,z,w
   end type pikg_f32vec4

   !**** PIKG::F64vec2
   type, public, bind(c) :: pikg_f64vec2
      real(kind=c_double) :: x,y
   end type pikg_f64vec2

   !**** PIKG::F64vec3
   type, public, bind(c) :: pikg_f64vec3
      real(kind=c_double) :: x,y,z
   end type pikg_f64vec3

   !**** PIKG::F64vec4
   type, public, bind(c) :: pikg_f64vec4
      real(kind=c_double) :: x,y,z,w
   end type pikg_f64vec4

   !* Private routines
   private :: set_2d_v32_i   ! i := c_int
   private :: set_2d_v32_f   ! f := c_float
   private :: set_2d_v32_fa  ! a := array
   private :: set_2d_v32_d   ! d := c_double
   private :: set_2d_v32_da
   private :: set_2d_v32_v32
   private :: set_2d_v32_v64
   private :: set_2d_v64_i
   private :: set_2d_v64_f
   private :: set_2d_v64_fa
   private :: set_2d_v64_d
   private :: set_2d_v64_da
   private :: set_2d_v64_v32
   private :: set_2d_v64_v64

   private :: set_3d_v32_i   
   private :: set_3d_v32_f   
   private :: set_3d_v32_fa  
   private :: set_3d_v32_d   
   private :: set_3d_v32_da
   private :: set_3d_v32_v32
   private :: set_3d_v32_v64
   private :: set_3d_v64_i
   private :: set_3d_v64_f
   private :: set_3d_v64_fa
   private :: set_3d_v64_d
   private :: set_3d_v64_da
   private :: set_3d_v64_v32
   private :: set_3d_v64_v64

   private :: set_4d_v32_i   
   private :: set_4d_v32_f   
   private :: set_4d_v32_fa  
   private :: set_4d_v32_d   
   private :: set_4d_v32_da
   private :: set_4d_v32_v32
   private :: set_4d_v32_v64
   private :: set_4d_v64_i
   private :: set_4d_v64_f
   private :: set_4d_v64_fa
   private :: set_4d_v64_d
   private :: set_4d_v64_da
   private :: set_4d_v64_v32
   private :: set_4d_v64_v64

   private :: add_2d_v32_fa
   private :: add_2d_fa_v32
   private :: add_2d_v32_v32
   private :: add_2d_v32_da
   private :: add_2d_da_v32
   private :: add_2d_v64_da
   private :: add_2d_da_v64
   private :: add_2d_v32_v64
   private :: add_2d_v64_v32
   private :: add_2d_v64_v64
   private :: pos_2d_v32
   private :: pos_2d_v64

   private :: add_3d_v32_fa
   private :: add_3d_fa_v32
   private :: add_3d_v32_v32
   private :: add_3d_v32_da
   private :: add_3d_da_v32
   private :: add_3d_v64_da
   private :: add_3d_da_v64
   private :: add_3d_v32_v64
   private :: add_3d_v64_v32
   private :: add_3d_v64_v64
   private :: pos_3d_v32
   private :: pos_3d_v64

   private :: add_4d_v32_fa
   private :: add_4d_fa_v32
   private :: add_4d_v32_v32
   private :: add_4d_v32_da
   private :: add_4d_da_v32
   private :: add_4d_v64_da
   private :: add_4d_da_v64
   private :: add_4d_v32_v64
   private :: add_4d_v64_v32
   private :: add_4d_v64_v64
   private :: pos_4d_v32
   private :: pos_4d_v64

   private :: sub_2d_v32_fa
   private :: sub_2d_fa_v32
   private :: sub_2d_v32_v32
   private :: sub_2d_v32_da
   private :: sub_2d_da_v32
   private :: sub_2d_v64_da
   private :: sub_2d_da_v64
   private :: sub_2d_v32_v64
   private :: sub_2d_v64_v32
   private :: sub_2d_v64_v64
   private :: neg_2d_v32
   private :: neg_2d_v64

   private :: sub_3d_v32_fa
   private :: sub_3d_fa_v32
   private :: sub_3d_v32_v32
   private :: sub_3d_v32_da
   private :: sub_3d_da_v32
   private :: sub_3d_v64_da
   private :: sub_3d_da_v64
   private :: sub_3d_v32_v64
   private :: sub_3d_v64_v32
   private :: sub_3d_v64_v64
   private :: neg_3d_v32
   private :: neg_3d_v64

   private :: sub_4d_v32_fa
   private :: sub_4d_fa_v32
   private :: sub_4d_v32_v32
   private :: sub_4d_v32_da
   private :: sub_4d_da_v32
   private :: sub_4d_v64_da
   private :: sub_4d_da_v64
   private :: sub_4d_v32_v64
   private :: sub_4d_v64_v32
   private :: sub_4d_v64_v64
   private :: neg_4d_v32
   private :: neg_4d_v64

   private :: mul_2d_v32_i
   private :: mul_2d_v32_f
   private :: mul_2d_v32_d
   private :: mul_2d_i_v32
   private :: mul_2d_f_v32
   private :: mul_2d_d_v32
   private :: mul_2d_v64_i
   private :: mul_2d_v64_f
   private :: mul_2d_v64_d
   private :: mul_2d_i_v64
   private :: mul_2d_f_v64
   private :: mul_2d_d_v64
   private :: mul_2d_v32_fa
   private :: mul_2d_fa_v32
   private :: mul_2d_v32_v32
   private :: mul_2d_v32_da
   private :: mul_2d_da_v32
   private :: mul_2d_v64_da
   private :: mul_2d_da_v64
   private :: mul_2d_v64_v32
   private :: mul_2d_v32_v64
   private :: mul_2d_v64_v64

   private :: mul_3d_v32_i
   private :: mul_3d_v32_f
   private :: mul_3d_v32_d
   private :: mul_3d_i_v32
   private :: mul_3d_f_v32
   private :: mul_3d_d_v32
   private :: mul_3d_v64_i
   private :: mul_3d_v64_f
   private :: mul_3d_v64_d
   private :: mul_3d_i_v64
   private :: mul_3d_f_v64
   private :: mul_3d_d_v64
   private :: mul_3d_v32_fa
   private :: mul_3d_fa_v32
   private :: mul_3d_v32_v32
   private :: mul_3d_v32_da
   private :: mul_3d_da_v32
   private :: mul_3d_v64_da
   private :: mul_3d_da_v64
   private :: mul_3d_v64_v32
   private :: mul_3d_v32_v64
   private :: mul_3d_v64_v64

   private :: mul_4d_v32_i
   private :: mul_4d_v32_f
   private :: mul_4d_v32_d
   private :: mul_4d_i_v32
   private :: mul_4d_f_v32
   private :: mul_4d_d_v32
   private :: mul_4d_v64_i
   private :: mul_4d_v64_f
   private :: mul_4d_v64_d
   private :: mul_4d_i_v64
   private :: mul_4d_f_v64
   private :: mul_4d_d_v64
   private :: mul_4d_v32_fa
   private :: mul_4d_fa_v32
   private :: mul_4d_v32_v32
   private :: mul_4d_v32_da
   private :: mul_4d_da_v32
   private :: mul_4d_v64_da
   private :: mul_4d_da_v64
   private :: mul_4d_v64_v32
   private :: mul_4d_v32_v64
   private :: mul_4d_v64_v64

   private :: div_2d_v32_i
   private :: div_2d_v32_f
   private :: div_2d_v32_d
   private :: div_2d_v64_i
   private :: div_2d_v64_f
   private :: div_2d_v64_d

   private :: div_3d_v32_i
   private :: div_3d_v32_f
   private :: div_3d_v32_d
   private :: div_3d_v64_i
   private :: div_3d_v64_f
   private :: div_3d_v64_d

   private :: div_4d_v32_i
   private :: div_4d_v32_f
   private :: div_4d_v32_d
   private :: div_4d_v64_i
   private :: div_4d_v64_f
   private :: div_4d_v64_d

   !* Declare interface operators
   interface assignment (=)
      module procedure set_2d_v32_i
      module procedure set_2d_v32_f
      module procedure set_2d_v32_fa
      module procedure set_2d_v32_d
      module procedure set_2d_v32_da
      module procedure set_2d_v32_v32
      module procedure set_2d_v32_v64
      module procedure set_2d_v64_i
      module procedure set_2d_v64_f
      module procedure set_2d_v64_fa
      module procedure set_2d_v64_d
      module procedure set_2d_v64_da
      module procedure set_2d_v64_v32
      module procedure set_2d_v64_v64

      module procedure set_3d_v32_i
      module procedure set_3d_v32_f
      module procedure set_3d_v32_fa
      module procedure set_3d_v32_d
      module procedure set_3d_v32_da
      module procedure set_3d_v32_v32
      module procedure set_3d_v32_v64
      module procedure set_3d_v64_i
      module procedure set_3d_v64_f
      module procedure set_3d_v64_fa
      module procedure set_3d_v64_d
      module procedure set_3d_v64_da
      module procedure set_3d_v64_v32
      module procedure set_3d_v64_v64

      module procedure set_4d_v32_i
      module procedure set_4d_v32_f
      module procedure set_4d_v32_fa
      module procedure set_4d_v32_d
      module procedure set_4d_v32_da
      module procedure set_4d_v32_v32
      module procedure set_4d_v32_v64
      module procedure set_4d_v64_i
      module procedure set_4d_v64_f
      module procedure set_4d_v64_fa
      module procedure set_4d_v64_d
      module procedure set_4d_v64_da
      module procedure set_4d_v64_v32
      module procedure set_4d_v64_v64
   end interface 

   interface operator (+)
      module procedure add_2d_v32_fa
      module procedure add_2d_fa_v32
      module procedure add_2d_v32_v32
      module procedure add_2d_v32_da
      module procedure add_2d_da_v32
      module procedure add_2d_v64_da
      module procedure add_2d_da_v64
      module procedure add_2d_v32_v64
      module procedure add_2d_v64_v32
      module procedure add_2d_v64_v64
      module procedure pos_2d_v32
      module procedure pos_2d_v64

      module procedure add_3d_v32_fa
      module procedure add_3d_fa_v32
      module procedure add_3d_v32_v32
      module procedure add_3d_v32_da
      module procedure add_3d_da_v32
      module procedure add_3d_v64_da
      module procedure add_3d_da_v64
      module procedure add_3d_v32_v64
      module procedure add_3d_v64_v32
      module procedure add_3d_v64_v64
      module procedure pos_3d_v32
      module procedure pos_3d_v64

      module procedure add_4d_v32_fa
      module procedure add_4d_fa_v32
      module procedure add_4d_v32_v32
      module procedure add_4d_v32_da
      module procedure add_4d_da_v32
      module procedure add_4d_v64_da
      module procedure add_4d_da_v64
      module procedure add_4d_v32_v64
      module procedure add_4d_v64_v32
      module procedure add_4d_v64_v64
      module procedure pos_4d_v32
      module procedure pos_4d_v64
   end interface 

   interface operator (-)
      module procedure sub_2d_v32_fa
      module procedure sub_2d_fa_v32
      module procedure sub_2d_v32_v32
      module procedure sub_2d_v32_da
      module procedure sub_2d_da_v32
      module procedure sub_2d_v64_da
      module procedure sub_2d_da_v64
      module procedure sub_2d_v32_v64
      module procedure sub_2d_v64_v32
      module procedure sub_2d_v64_v64
      module procedure neg_2d_v32
      module procedure neg_2d_v64

      module procedure sub_3d_v32_fa
      module procedure sub_3d_fa_v32
      module procedure sub_3d_v32_v32
      module procedure sub_3d_v32_da
      module procedure sub_3d_da_v32
      module procedure sub_3d_v64_da
      module procedure sub_3d_da_v64
      module procedure sub_3d_v32_v64
      module procedure sub_3d_v64_v32
      module procedure sub_3d_v64_v64
      module procedure neg_3d_v32
      module procedure neg_3d_v64

      module procedure sub_4d_v32_fa
      module procedure sub_4d_fa_v32
      module procedure sub_4d_v32_v32
      module procedure sub_4d_v32_da
      module procedure sub_4d_da_v32
      module procedure sub_4d_v64_da
      module procedure sub_4d_da_v64
      module procedure sub_4d_v32_v64
      module procedure sub_4d_v64_v32
      module procedure sub_4d_v64_v64
      module procedure neg_4d_v32
      module procedure neg_4d_v64
   end interface

   interface operator (*)
      module procedure mul_2d_v32_i
      module procedure mul_2d_v32_f
      module procedure mul_2d_v32_d
      module procedure mul_2d_i_v32
      module procedure mul_2d_f_v32
      module procedure mul_2d_d_v32
      module procedure mul_2d_v64_i
      module procedure mul_2d_v64_f
      module procedure mul_2d_v64_d
      module procedure mul_2d_i_v64
      module procedure mul_2d_f_v64
      module procedure mul_2d_d_v64
      module procedure mul_2d_v32_fa
      module procedure mul_2d_fa_v32
      module procedure mul_2d_v32_v32
      module procedure mul_2d_v32_da
      module procedure mul_2d_da_v32
      module procedure mul_2d_v64_da
      module procedure mul_2d_da_v64
      module procedure mul_2d_v64_v32
      module procedure mul_2d_v32_v64
      module procedure mul_2d_v64_v64

      module procedure mul_3d_v32_i
      module procedure mul_3d_v32_f
      module procedure mul_3d_v32_d
      module procedure mul_3d_i_v32
      module procedure mul_3d_f_v32
      module procedure mul_3d_d_v32
      module procedure mul_3d_v64_i
      module procedure mul_3d_v64_f
      module procedure mul_3d_v64_d
      module procedure mul_3d_i_v64
      module procedure mul_3d_f_v64
      module procedure mul_3d_d_v64
      module procedure mul_3d_v32_fa
      module procedure mul_3d_fa_v32
      module procedure mul_3d_v32_v32
      module procedure mul_3d_v32_da
      module procedure mul_3d_da_v32
      module procedure mul_3d_v64_da
      module procedure mul_3d_da_v64
      module procedure mul_3d_v64_v32
      module procedure mul_3d_v32_v64
      module procedure mul_3d_v64_v64

      module procedure mul_4d_v32_i
      module procedure mul_4d_v32_f
      module procedure mul_4d_v32_d
      module procedure mul_4d_i_v32
      module procedure mul_4d_f_v32
      module procedure mul_4d_d_v32
      module procedure mul_4d_v64_i
      module procedure mul_4d_v64_f
      module procedure mul_4d_v64_d
      module procedure mul_4d_i_v64
      module procedure mul_4d_f_v64
      module procedure mul_4d_d_v64
      module procedure mul_4d_v32_fa
      module procedure mul_4d_fa_v32
      module procedure mul_4d_v32_v32
      module procedure mul_4d_v32_da
      module procedure mul_4d_da_v32
      module procedure mul_4d_v64_da
      module procedure mul_4d_da_v64
      module procedure mul_4d_v64_v32
      module procedure mul_4d_v32_v64
      module procedure mul_4d_v64_v64
   end interface

   interface operator (/)
      module procedure div_2d_v32_i
      module procedure div_2d_v32_f
      module procedure div_2d_v32_d
      module procedure div_2d_v64_i
      module procedure div_2d_v64_f
      module procedure div_2d_v64_d

      module procedure div_3d_v32_i
      module procedure div_3d_v32_f
      module procedure div_3d_v32_d
      module procedure div_3d_v64_i
      module procedure div_3d_v64_f
      module procedure div_3d_v64_d

      module procedure div_4d_v32_i
      module procedure div_4d_v32_f
      module procedure div_4d_v32_d
      module procedure div_4d_v64_i
      module procedure div_4d_v64_f
      module procedure div_4d_v64_d
   end interface

   contains

   !####################################################################
   !========
   !   2D
   !========
   elemental pure subroutine set_2d_v32_i(v,i)
      type(pikg_f32vec2), intent(out) :: v
      integer(kind=c_int), intent(in) :: i
      v%x = real(i,kind=c_float)
      v%y = real(i,kind=c_float)
   end subroutine set_2d_v32_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v32_f(v,f)
      type(pikg_f32vec2), intent(out) :: v
      real(kind=c_float), intent(in) :: f
      v%x = f
      v%y = f
   end subroutine set_2d_v32_f

   !--------------------------------------------------------------------
   pure subroutine set_2d_v32_fa(v,f)
      type(pikg_f32vec2), intent(out) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      v%x = f(1)
      v%y = f(2)
   end subroutine set_2d_v32_fa

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v32_d(v,d)
      type(pikg_f32vec2), intent(out) :: v
      real(kind=c_double), intent(in) :: d
      v%x = real(d,kind=c_float)
      v%y = real(d,kind=c_float)
   end subroutine set_2d_v32_d

   !--------------------------------------------------------------------
   pure subroutine set_2d_v32_da(v,d)
      type(pikg_f32vec2), intent(out) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      v%x = real(d(1),kind=c_float)
      v%y = real(d(2),kind=c_float)
   end subroutine set_2d_v32_da

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v32_v32(a,b)
      type(pikg_f32vec2), intent(out) :: a
      type(pikg_f32vec2), intent(in) :: b
      a%x = b%x
      a%y = b%y
   end subroutine set_2d_v32_v32

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v32_v64(a,b)
      type(pikg_f32vec2), intent(out) :: a
      type(pikg_f64vec2), intent(in) :: b
      a%x = real(b%x,kind=c_float)
      a%y = real(b%y,kind=c_float)
   end subroutine set_2d_v32_v64

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v64_i(v,i)
      type(pikg_f64vec2), intent(out) :: v
      integer(kind=c_int), intent(in) :: i
      v%x = real(i,kind=c_double)
      v%y = real(i,kind=c_double)
   end subroutine set_2d_v64_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v64_f(v,f)
      type(pikg_f64vec2), intent(out) :: v
      real(kind=c_float), intent(in) :: f
      v%x = real(f,kind=c_double)
      v%y = real(f,kind=c_double)
   end subroutine set_2d_v64_f

   !--------------------------------------------------------------------
   pure subroutine set_2d_v64_fa(v,f)
      type(pikg_f64vec2), intent(out) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      v%x = real(f(1),kind=c_double)
      v%y = real(f(2),kind=c_double)
   end subroutine set_2d_v64_fa

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v64_d(v,d)
      type(pikg_f64vec2), intent(out) :: v
      real(kind=c_double), intent(in) :: d
      v%x = d
      v%y = d
   end subroutine set_2d_v64_d

   !--------------------------------------------------------------------
   pure subroutine set_2d_v64_da(v,d)
      type(pikg_f64vec2), intent(out) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      v%x = d(1)
      v%y = d(2)
   end subroutine set_2d_v64_da

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v64_v32(a,b)
      type(pikg_f64vec2), intent(out) :: a
      type(pikg_f32vec2), intent(in) :: b
      a%x = real(b%x,kind=c_double)
      a%y = real(b%y,kind=c_double)
   end subroutine set_2d_v64_v32

   !--------------------------------------------------------------------
   elemental pure subroutine set_2d_v64_v64(a,b)
      type(pikg_f64vec2), intent(out) :: a
      type(pikg_f64vec2), intent(in) :: b
      a%x = b%x
      a%y = b%y
   end subroutine set_2d_v64_v64
   !========
   !   3D
   !========
   elemental pure subroutine set_3d_v32_i(v,i)
      type(pikg_f32vec3), intent(out) :: v
      integer(kind=c_int), intent(in) :: i
      v%x = real(i,kind=c_float)
      v%y = real(i,kind=c_float)
      v%z = real(i,kind=c_float)
   end subroutine set_3d_v32_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v32_f(v,f)
      type(pikg_f32vec3), intent(out) :: v
      real(kind=c_float), intent(in) :: f
      v%x = f
      v%y = f
      v%z = f
   end subroutine set_3d_v32_f

   !--------------------------------------------------------------------
   pure subroutine set_3d_v32_fa(v,f)
      type(pikg_f32vec3), intent(out) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      v%x = f(1)
      v%y = f(2)
      v%z = f(3)
   end subroutine set_3d_v32_fa

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v32_d(v,d)
      type(pikg_f32vec3), intent(out) :: v
      real(kind=c_double), intent(in) :: d
      v%x = real(d,kind=c_float)
      v%y = real(d,kind=c_float)
      v%z = real(d,kind=c_float)
   end subroutine set_3d_v32_d

   !--------------------------------------------------------------------
   pure subroutine set_3d_v32_da(v,d)
      type(pikg_f32vec3), intent(out) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      v%x = real(d(1),kind=c_float)
      v%y = real(d(2),kind=c_float)
      v%z = real(d(3),kind=c_float)
   end subroutine set_3d_v32_da

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v32_v32(a,b)
      type(pikg_f32vec3), intent(out) :: a
      type(pikg_f32vec3), intent(in) :: b
      a%x = b%x
      a%y = b%y
      a%z = b%z
   end subroutine set_3d_v32_v32

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v32_v64(a,b)
      type(pikg_f32vec3), intent(out) :: a
      type(pikg_f64vec3), intent(in) :: b
      a%x = real(b%x,kind=c_float)
      a%y = real(b%y,kind=c_float)
      a%z = real(b%z,kind=c_float)
   end subroutine set_3d_v32_v64

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v64_i(v,i)
      type(pikg_f64vec3), intent(out) :: v
      integer(kind=c_int), intent(in) :: i
      v%x = real(i,kind=c_double)
      v%y = real(i,kind=c_double)
      v%z = real(i,kind=c_double)
   end subroutine set_3d_v64_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v64_f(v,f)
      type(pikg_f64vec3), intent(out) :: v
      real(kind=c_float), intent(in) :: f
      v%x = real(f,kind=c_double)
      v%y = real(f,kind=c_double)
      v%z = real(f,kind=c_double)
   end subroutine set_3d_v64_f

   !--------------------------------------------------------------------
   pure subroutine set_3d_v64_fa(v,f)
      type(pikg_f64vec3), intent(out) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      v%x = real(f(1),kind=c_double)
      v%y = real(f(2),kind=c_double)
      v%z = real(f(3),kind=c_double)
   end subroutine set_3d_v64_fa

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v64_d(v,d)
      type(pikg_f64vec3), intent(out) :: v
      real(kind=c_double), intent(in) :: d
      v%x = d
      v%y = d
      v%z = d
   end subroutine set_3d_v64_d

   !--------------------------------------------------------------------
   pure subroutine set_3d_v64_da(v,d)
      type(pikg_f64vec3), intent(out) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      v%x = d(1)
      v%y = d(2)
      v%z = d(3)
   end subroutine set_3d_v64_da

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v64_v32(a,b)
      type(pikg_f64vec3), intent(out) :: a
      type(pikg_f32vec3), intent(in) :: b
      a%x = real(b%x,kind=c_double)
      a%y = real(b%y,kind=c_double)
      a%z = real(b%z,kind=c_double)
   end subroutine set_3d_v64_v32

   !--------------------------------------------------------------------
   elemental pure subroutine set_3d_v64_v64(a,b)
      type(pikg_f64vec3), intent(out) :: a
      type(pikg_f64vec3), intent(in) :: b
      a%x = b%x
      a%y = b%y
      a%z = b%z
   end subroutine set_3d_v64_v64
   !========
   !   4D
   !========
   elemental pure subroutine set_4d_v32_i(v,i)
      type(pikg_f32vec4), intent(out) :: v
      integer(kind=c_int), intent(in) :: i
      v%x = real(i,kind=c_float)
      v%y = real(i,kind=c_float)
      v%z = real(i,kind=c_float)
      v%w = real(i,kind=c_float)
   end subroutine set_4d_v32_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v32_f(v,f)
      type(pikg_f32vec4), intent(out) :: v
      real(kind=c_float), intent(in) :: f
      v%x = f
      v%y = f
      v%z = f
      v%w = f
   end subroutine set_4d_v32_f

   !--------------------------------------------------------------------
   pure subroutine set_4d_v32_fa(v,f)
      type(pikg_f32vec4), intent(out) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      v%x = f(1)
      v%y = f(2)
      v%z = f(3)
      v%w = f(4)
   end subroutine set_4d_v32_fa

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v32_d(v,d)
      type(pikg_f32vec4), intent(out) :: v
      real(kind=c_double), intent(in) :: d
      v%x = real(d,kind=c_float)
      v%y = real(d,kind=c_float)
      v%z = real(d,kind=c_float)
      v%w = real(d,kind=c_float)
   end subroutine set_4d_v32_d

   !--------------------------------------------------------------------
   pure subroutine set_4d_v32_da(v,d)
      type(pikg_f32vec4), intent(out) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      v%x = real(d(1),kind=c_float)
      v%y = real(d(2),kind=c_float)
      v%z = real(d(3),kind=c_float)
      v%w = real(d(4),kind=c_float)
   end subroutine set_4d_v32_da

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v32_v32(a,b)
      type(pikg_f32vec4), intent(out) :: a
      type(pikg_f32vec4), intent(in) :: b
      a%x = b%x
      a%y = b%y
      a%z = b%z
      a%w = b%w
   end subroutine set_4d_v32_v32

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v32_v64(a,b)
      type(pikg_f32vec4), intent(out) :: a
      type(pikg_f64vec4), intent(in) :: b
      a%x = real(b%x,kind=c_float)
      a%y = real(b%y,kind=c_float)
      a%z = real(b%z,kind=c_float)
      a%w = real(b%w,kind=c_float)
   end subroutine set_4d_v32_v64

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v64_i(v,i)
      type(pikg_f64vec4), intent(out) :: v
      integer(kind=c_int), intent(in) :: i
      v%x = real(i,kind=c_double)
      v%y = real(i,kind=c_double)
      v%z = real(i,kind=c_double)
      v%z = real(i,kind=c_double)
   end subroutine set_4d_v64_i

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v64_f(v,f)
      type(pikg_f64vec4), intent(out) :: v
      real(kind=c_float), intent(in) :: f
      v%x = real(f,kind=c_double)
      v%y = real(f,kind=c_double)
      v%z = real(f,kind=c_double)
      v%w = real(f,kind=c_double)
   end subroutine set_4d_v64_f

   !--------------------------------------------------------------------
   pure subroutine set_4d_v64_fa(v,f)
      type(pikg_f64vec4), intent(out) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      v%x = real(f(1),kind=c_double)
      v%y = real(f(2),kind=c_double)
      v%z = real(f(3),kind=c_double)
      v%w = real(f(4),kind=c_double)
   end subroutine set_4d_v64_fa

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v64_d(v,d)
      type(pikg_f64vec4), intent(out) :: v
      real(kind=c_double), intent(in) :: d
      v%x = d
      v%y = d
      v%z = d
      v%w = d
   end subroutine set_4d_v64_d

   !--------------------------------------------------------------------
   pure subroutine set_4d_v64_da(v,d)
      type(pikg_f64vec4), intent(out) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      v%x = d(1)
      v%y = d(2)
      v%z = d(3)
      v%w = d(4)
   end subroutine set_4d_v64_da

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v64_v32(a,b)
      type(pikg_f64vec4), intent(out) :: a
      type(pikg_f32vec4), intent(in) :: b
      a%x = real(b%x,kind=c_double)
      a%y = real(b%y,kind=c_double)
      a%z = real(b%z,kind=c_double)
      a%w = real(b%w,kind=c_double)
   end subroutine set_4d_v64_v32

   !--------------------------------------------------------------------
   elemental pure subroutine set_4d_v64_v64(a,b)
      type(pikg_f64vec4), intent(out) :: a
      type(pikg_f64vec4), intent(in) :: b
      a%x = b%x
      a%y = b%y
      a%z = b%z
      a%w = b%w
   end subroutine set_4d_v64_v64

   !####################################################################
   !========
   !   2D
   !========
   pure type(pikg_f32vec2) function add_2d_v32_fa(v,f)
     type(pikg_f32vec2), intent(in) :: v
     real(kind=c_float), dimension(:), intent(IN) :: f
     add_2d_v32_fa%x = v%x + f(1)
     add_2d_v32_fa%y = v%y + f(2)
   end function add_2d_v32_fa

   !--------------------------------------------------------------------
   pure type(pikg_f32vec2) function add_2d_fa_v32(f,v)
     real(kind=c_float), dimension(:), intent(IN) :: f
     type(pikg_f32vec2), intent(in) :: v
     add_2d_fa_v32%x = f(1) + v%x
     add_2d_fa_v32%y = f(2) + v%y
   end function add_2d_fa_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function add_2d_v32_v32(a,b)
     type(pikg_f32vec2), intent(in) :: a,b
     add_2d_v32_v32%x = a%x + b%x
     add_2d_v32_v32%y = a%y + b%y
   end function add_2d_v32_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function add_2d_v32_da(v,d)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      add_2d_v32_da%x = v%x + d(1)
      add_2d_v32_da%y = v%y + d(2)
   end function add_2d_v32_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function add_2d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec2), intent(in) :: v
      add_2d_da_v32%x = d(1) + v%x 
      add_2d_da_v32%y = d(2) + v%y
   end function add_2d_da_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function add_2d_v64_da(v,d)
      type(pikg_f64vec2), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      add_2d_v64_da%x = v%x + d(1)
      add_2d_v64_da%y = v%y + d(2)
   end function add_2d_v64_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function add_2d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec2), intent(in) :: v
      add_2d_da_v64%x = d(1) + v%x 
      add_2d_da_v64%y = d(2) + v%y
   end function add_2d_da_v64

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function add_2d_v32_v64(a,b)
      type(pikg_f32vec2), intent(in) :: a
      type(pikg_f64vec2), intent(in) :: b
      add_2d_v32_v64%x = a%x + b%x
      add_2d_v32_v64%y = a%y + b%y
   end function add_2d_v32_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function add_2d_v64_v32(a,b)
      type(pikg_f64vec2), intent(in) :: a
      type(pikg_f32vec2), intent(in) :: b
      add_2d_v64_v32%x = a%x + b%x
      add_2d_v64_v32%y = a%y + b%y
   end function add_2d_v64_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function add_2d_v64_v64(a,b)
      type(pikg_f64vec2), intent(in) :: a,b
      add_2d_v64_v64%x = a%x + b%x
      add_2d_v64_v64%y = a%y + b%y
   end function add_2d_v64_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function pos_2d_v32(v)
      type(pikg_f32vec2), intent(in) :: v
      pos_2d_v32%x = v%x
      pos_2d_v32%y = v%y
   end function pos_2d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function pos_2d_v64(v)
      type(pikg_f64vec2), intent(in) :: v
      pos_2d_v64%x = v%x
      pos_2d_v64%y = v%y
   end function pos_2d_v64
   !========
   !   3D
   !========
   pure type(pikg_f32vec3) function add_3d_v32_fa(v,f)
     type(pikg_f32vec3), intent(in) :: v
     real(kind=c_float), dimension(:), intent(IN) :: f
     add_3d_v32_fa%x = v%x + f(1)
     add_3d_v32_fa%y = v%y + f(2)
     add_3d_v32_fa%z = v%z + f(3)
   end function add_3d_v32_fa

   !--------------------------------------------------------------------
   pure type(pikg_f32vec3) function add_3d_fa_v32(f,v)
     real(kind=c_float), dimension(:), intent(IN) :: f
     type(pikg_f32vec3), intent(in) :: v
     add_3d_fa_v32%x = f(1) + v%x
     add_3d_fa_v32%y = f(2) + v%y
     add_3d_fa_v32%z = f(3) + v%z
   end function add_3d_fa_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function add_3d_v32_v32(a,b)
     type(pikg_f32vec3), intent(in) :: a,b
     add_3d_v32_v32%x = a%x + b%x
     add_3d_v32_v32%y = a%y + b%y
     add_3d_v32_v32%z = a%z + b%z
   end function add_3d_v32_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function add_3d_v32_da(v,d)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      add_3d_v32_da%x = v%x + d(1)
      add_3d_v32_da%y = v%y + d(2)
      add_3d_v32_da%z = v%z + d(3)
   end function add_3d_v32_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function add_3d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec3), intent(in) :: v
      add_3d_da_v32%x = d(1) + v%x 
      add_3d_da_v32%y = d(2) + v%y
      add_3d_da_v32%z = d(3) + v%z
   end function add_3d_da_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function add_3d_v64_da(v,d)
      type(pikg_f64vec3), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      add_3d_v64_da%x = v%x + d(1)
      add_3d_v64_da%y = v%y + d(2)
      add_3d_v64_da%z = v%z + d(3)
   end function add_3d_v64_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function add_3d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec3), intent(in) :: v
      add_3d_da_v64%x = d(1) + v%x 
      add_3d_da_v64%y = d(2) + v%y
      add_3d_da_v64%z = d(3) + v%z
   end function add_3d_da_v64

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function add_3d_v32_v64(a,b)
      type(pikg_f32vec3), intent(in) :: a
      type(pikg_f64vec3), intent(in) :: b
      add_3d_v32_v64%x = a%x + b%x
      add_3d_v32_v64%y = a%y + b%y
      add_3d_v32_v64%z = a%z + b%z
   end function add_3d_v32_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function add_3d_v64_v32(a,b)
      type(pikg_f64vec3), intent(in) :: a
      type(pikg_f32vec3), intent(in) :: b
      add_3d_v64_v32%x = a%x + b%x
      add_3d_v64_v32%y = a%y + b%y
      add_3d_v64_v32%z = a%z + b%z
   end function add_3d_v64_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function add_3d_v64_v64(a,b)
      type(pikg_f64vec3), intent(in) :: a,b
      add_3d_v64_v64%x = a%x + b%x
      add_3d_v64_v64%y = a%y + b%y
      add_3d_v64_v64%z = a%z + b%z
   end function add_3d_v64_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function pos_3d_v32(v)
      type(pikg_f32vec3), intent(in) :: v
      pos_3d_v32%x = v%x
      pos_3d_v32%y = v%y
      pos_3d_v32%z = v%z
   end function pos_3d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function pos_3d_v64(v)
      type(pikg_f64vec3), intent(in) :: v
      pos_3d_v64%x = v%x
      pos_3d_v64%y = v%y
      pos_3d_v64%z = v%z
   end function pos_3d_v64
   !========
   !   4D
   !========
   pure type(pikg_f32vec4) function add_4d_v32_fa(v,f)
     type(pikg_f32vec4), intent(in) :: v
     real(kind=c_float), dimension(:), intent(IN) :: f
     add_4d_v32_fa%x = v%x + f(1)
     add_4d_v32_fa%y = v%y + f(2)
     add_4d_v32_fa%z = v%z + f(3)
     add_4d_v32_fa%w = v%w + f(4)
   end function add_4d_v32_fa

   !--------------------------------------------------------------------
   pure type(pikg_f32vec4) function add_4d_fa_v32(f,v)
     real(kind=c_float), dimension(:), intent(IN) :: f
     type(pikg_f32vec4), intent(in) :: v
     add_4d_fa_v32%x = f(1) + v%x
     add_4d_fa_v32%y = f(2) + v%y
     add_4d_fa_v32%z = f(3) + v%z
     add_4d_fa_v32%w = f(4) + v%w
   end function add_4d_fa_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function add_4d_v32_v32(a,b)
     type(pikg_f32vec4), intent(in) :: a,b
     add_4d_v32_v32%x = a%x + b%x
     add_4d_v32_v32%y = a%y + b%y
     add_4d_v32_v32%z = a%z + b%z
     add_4d_v32_v32%w = a%w + b%w
   end function add_4d_v32_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function add_4d_v32_da(v,d)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      add_4d_v32_da%x = v%x + d(1)
      add_4d_v32_da%y = v%y + d(2)
      add_4d_v32_da%z = v%z + d(3)
      add_4d_v32_da%w = v%w + d(4)
   end function add_4d_v32_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function add_4d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec4), intent(in) :: v
      add_4d_da_v32%x = d(1) + v%x 
      add_4d_da_v32%y = d(2) + v%y
      add_4d_da_v32%z = d(3) + v%z
      add_4d_da_v32%w = d(4) + v%w
   end function add_4d_da_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function add_4d_v64_da(v,d)
      type(pikg_f64vec4), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      add_4d_v64_da%x = v%x + d(1)
      add_4d_v64_da%y = v%y + d(2)
      add_4d_v64_da%z = v%z + d(3)
      add_4d_v64_da%w = v%w + d(4)
   end function add_4d_v64_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function add_4d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec4), intent(in) :: v
      add_4d_da_v64%x = d(1) + v%x 
      add_4d_da_v64%y = d(2) + v%y
      add_4d_da_v64%z = d(3) + v%z
      add_4d_da_v64%w = d(4) + v%w
   end function add_4d_da_v64

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function add_4d_v32_v64(a,b)
      type(pikg_f32vec4), intent(in) :: a
      type(pikg_f64vec4), intent(in) :: b
      add_4d_v32_v64%x = a%x + b%x
      add_4d_v32_v64%y = a%y + b%y
      add_4d_v32_v64%z = a%z + b%z
      add_4d_v32_v64%w = a%w + b%w
   end function add_4d_v32_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function add_4d_v64_v32(a,b)
      type(pikg_f64vec4), intent(in) :: a
      type(pikg_f32vec4), intent(in) :: b
      add_4d_v64_v32%x = a%x + b%x
      add_4d_v64_v32%y = a%y + b%y
      add_4d_v64_v32%z = a%z + b%z
      add_4d_v64_v32%w = a%w + b%w
   end function add_4d_v64_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function add_4d_v64_v64(a,b)
      type(pikg_f64vec4), intent(in) :: a,b
      add_4d_v64_v64%x = a%x + b%x
      add_4d_v64_v64%y = a%y + b%y
      add_4d_v64_v64%z = a%z + b%z
      add_4d_v64_v64%w = a%w + b%w
   end function add_4d_v64_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function pos_4d_v32(v)
      type(pikg_f32vec4), intent(in) :: v
      pos_4d_v32%x = v%x
      pos_4d_v32%y = v%y
      pos_4d_v32%z = v%z
      pos_4d_v32%w = v%w
   end function pos_4d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function pos_4d_v64(v)
      type(pikg_f64vec4), intent(in) :: v
      pos_4d_v64%x = v%x
      pos_4d_v64%y = v%y
      pos_4d_v64%z = v%z
      pos_4d_v64%w = v%w
   end function pos_4d_v64

   !####################################################################
   !========
   !   2D
   !========
   pure type(pikg_f32vec2) function sub_2d_v32_fa(v,f)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      sub_2d_v32_fa%x = v%x - f(1)
      sub_2d_v32_fa%y = v%y - f(2)
   end function sub_2d_v32_fa

   !--------------------------------------------------------------------
   pure type(pikg_f32vec2) function sub_2d_fa_v32(f,v)
      real(kind=c_float), dimension(:), intent(in) :: f
      type(pikg_f32vec2), intent(in) :: v
      sub_2d_fa_v32%x = f(1) - v%x
      sub_2d_fa_v32%y = f(2) - v%y
   end function sub_2d_fa_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function sub_2d_v32_v32(a,b)
      type(pikg_f32vec2), intent(in) :: a,b
      sub_2d_v32_v32%x = a%x - b%x
      sub_2d_v32_v32%y = a%y - b%y
   end function sub_2d_v32_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function sub_2d_v32_da(v,d)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      sub_2d_v32_da%x = v%x - d(1)
      sub_2d_v32_da%y = v%y - d(2)
   end function sub_2d_v32_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function sub_2d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec2), intent(in) :: v
      sub_2d_da_v32%x = d(1) - v%x
      sub_2d_da_v32%y = d(2) - v%y
   end function sub_2d_da_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function sub_2d_v64_da(v,d)
      type(pikg_f64vec2), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      sub_2d_v64_da%x = v%x - d(1)
      sub_2d_v64_da%y = v%y - d(2)
   end function sub_2d_v64_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec2) function sub_2d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec2), intent(in) :: v
      sub_2d_da_v64%x = d(1) - v%x
      sub_2d_da_v64%y = d(2) - v%y
   end function sub_2d_da_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function sub_2d_v32_v64(a,b)
      type(pikg_f32vec2), intent(in) :: a
      type(pikg_f64vec2), intent(in) :: b
      sub_2d_v32_v64%x = a%x - b%x
      sub_2d_v32_v64%y = a%y - b%y
   end function sub_2d_v32_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function sub_2d_v64_v32(a,b)
      type(pikg_f64vec2), intent(in) :: a
      type(pikg_f32vec2), intent(in) :: b
      sub_2d_v64_v32%x = a%x - b%x
      sub_2d_v64_v32%y = a%y - b%y
   end function sub_2d_v64_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function sub_2d_v64_v64(a,b)
      type(pikg_f64vec2), intent(in) :: a,b
      sub_2d_v64_v64%x = a%x - b%x
      sub_2d_v64_v64%y = a%y - b%y
   end function sub_2d_v64_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function neg_2d_v32(v)
      type(pikg_f32vec2), intent(in) :: v
      neg_2d_v32%x = - v%x
      neg_2d_v32%y = - v%y
   end function neg_2d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function neg_2d_v64(v)
      type(pikg_f64vec2), intent(in) :: v
      neg_2d_v64%x = - v%x
      neg_2d_v64%y = - v%y
   end function neg_2d_v64
   !========
   !   3D
   !========
   pure type(pikg_f32vec3) function sub_3d_v32_fa(v,f)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      sub_3d_v32_fa%x = v%x - f(1)
      sub_3d_v32_fa%y = v%y - f(2)
      sub_3d_v32_fa%z = v%z - f(3)
   end function sub_3d_v32_fa

   !--------------------------------------------------------------------
   pure type(pikg_f32vec3) function sub_3d_fa_v32(f,v)
      real(kind=c_float), dimension(:), intent(in) :: f
      type(pikg_f32vec3), intent(in) :: v
      sub_3d_fa_v32%x = f(1) - v%x
      sub_3d_fa_v32%y = f(2) - v%y
      sub_3d_fa_v32%z = f(3) - v%z
   end function sub_3d_fa_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function sub_3d_v32_v32(a,b)
      type(pikg_f32vec3), intent(in) :: a,b
      sub_3d_v32_v32%x = a%x - b%x
      sub_3d_v32_v32%y = a%y - b%y
      sub_3d_v32_v32%z = a%z - b%z
   end function sub_3d_v32_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function sub_3d_v32_da(v,d)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      sub_3d_v32_da%x = v%x - d(1)
      sub_3d_v32_da%y = v%y - d(2)
      sub_3d_v32_da%z = v%z - d(3)
   end function sub_3d_v32_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function sub_3d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec3), intent(in) :: v
      sub_3d_da_v32%x = d(1) - v%x
      sub_3d_da_v32%y = d(2) - v%y
      sub_3d_da_v32%z = d(3) - v%z
   end function sub_3d_da_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function sub_3d_v64_da(v,d)
      type(pikg_f64vec3), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      sub_3d_v64_da%x = v%x - d(1)
      sub_3d_v64_da%y = v%y - d(2)
      sub_3d_v64_da%z = v%z - d(3)
   end function sub_3d_v64_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec3) function sub_3d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec3), intent(in) :: v
      sub_3d_da_v64%x = d(1) - v%x
      sub_3d_da_v64%y = d(2) - v%y
      sub_3d_da_v64%z = d(3) - v%z
   end function sub_3d_da_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function sub_3d_v32_v64(a,b)
      type(pikg_f32vec3), intent(in) :: a
      type(pikg_f64vec3), intent(in) :: b
      sub_3d_v32_v64%x = a%x - b%x
      sub_3d_v32_v64%y = a%y - b%y
      sub_3d_v32_v64%z = a%z - b%z
   end function sub_3d_v32_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function sub_3d_v64_v32(a,b)
      type(pikg_f64vec3), intent(in) :: a
      type(pikg_f32vec3), intent(in) :: b
      sub_3d_v64_v32%x = a%x - b%x
      sub_3d_v64_v32%y = a%y - b%y
      sub_3d_v64_v32%z = a%z - b%z
   end function sub_3d_v64_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function sub_3d_v64_v64(a,b)
      type(pikg_f64vec3), intent(in) :: a,b
      sub_3d_v64_v64%x = a%x - b%x
      sub_3d_v64_v64%y = a%y - b%y
      sub_3d_v64_v64%z = a%z - b%z
   end function sub_3d_v64_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function neg_3d_v32(v)
      type(pikg_f32vec3), intent(in) :: v
      neg_3d_v32%x = - v%x
      neg_3d_v32%y = - v%y
      neg_3d_v32%z = - v%z
   end function neg_3d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function neg_3d_v64(v)
      type(pikg_f64vec3), intent(in) :: v
      neg_3d_v64%x = - v%x
      neg_3d_v64%y = - v%y
      neg_3d_v64%z = - v%z
   end function neg_3d_v64
   !========
   !   4D
   !========
   pure type(pikg_f32vec4) function sub_4d_v32_fa(v,f)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      sub_4d_v32_fa%x = v%x - f(1)
      sub_4d_v32_fa%y = v%y - f(2)
      sub_4d_v32_fa%z = v%z - f(3)
      sub_4d_v32_fa%w = v%w - f(4)
   end function sub_4d_v32_fa

   !--------------------------------------------------------------------
   pure type(pikg_f32vec4) function sub_4d_fa_v32(f,v)
      real(kind=c_float), dimension(:), intent(in) :: f
      type(pikg_f32vec4), intent(in) :: v
      sub_4d_fa_v32%x = f(1) - v%x
      sub_4d_fa_v32%y = f(2) - v%y
      sub_4d_fa_v32%z = f(3) - v%z
      sub_4d_fa_v32%w = f(4) - v%w
   end function sub_4d_fa_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function sub_4d_v32_v32(a,b)
      type(pikg_f32vec4), intent(in) :: a,b
      sub_4d_v32_v32%x = a%x - b%x
      sub_4d_v32_v32%y = a%y - b%y
      sub_4d_v32_v32%z = a%z - b%z
      sub_4d_v32_v32%w = a%w - b%w
   end function sub_4d_v32_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function sub_4d_v32_da(v,d)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      sub_4d_v32_da%x = v%x - d(1)
      sub_4d_v32_da%y = v%y - d(2)
      sub_4d_v32_da%z = v%z - d(3)
      sub_4d_v32_da%w = v%w - d(4)
   end function sub_4d_v32_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function sub_4d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec4), intent(in) :: v
      sub_4d_da_v32%x = d(1) - v%x
      sub_4d_da_v32%y = d(2) - v%y
      sub_4d_da_v32%z = d(3) - v%z
      sub_4d_da_v32%w = d(4) - v%w
   end function sub_4d_da_v32

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function sub_4d_v64_da(v,d)
      type(pikg_f64vec4), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      sub_4d_v64_da%x = v%x - d(1)
      sub_4d_v64_da%y = v%y - d(2)
      sub_4d_v64_da%z = v%z - d(3)
      sub_4d_v64_da%w = v%w - d(4)
   end function sub_4d_v64_da

   !--------------------------------------------------------------------
   pure type(pikg_f64vec4) function sub_4d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec4), intent(in) :: v
      sub_4d_da_v64%x = d(1) - v%x
      sub_4d_da_v64%y = d(2) - v%y
      sub_4d_da_v64%z = d(3) - v%z
      sub_4d_da_v64%w = d(4) - v%w
   end function sub_4d_da_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function sub_4d_v32_v64(a,b)
      type(pikg_f32vec4), intent(in) :: a
      type(pikg_f64vec4), intent(in) :: b
      sub_4d_v32_v64%x = a%x - b%x
      sub_4d_v32_v64%y = a%y - b%y
      sub_4d_v32_v64%z = a%z - b%z
      sub_4d_v32_v64%w = a%w - b%w
   end function sub_4d_v32_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function sub_4d_v64_v32(a,b)
      type(pikg_f64vec4), intent(in) :: a
      type(pikg_f32vec4), intent(in) :: b
      sub_4d_v64_v32%x = a%x - b%x
      sub_4d_v64_v32%y = a%y - b%y
      sub_4d_v64_v32%z = a%z - b%z
      sub_4d_v64_v32%w = a%w - b%w
   end function sub_4d_v64_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function sub_4d_v64_v64(a,b)
      type(pikg_f64vec4), intent(in) :: a,b
      sub_4d_v64_v64%x = a%x - b%x
      sub_4d_v64_v64%y = a%y - b%y
      sub_4d_v64_v64%z = a%z - b%z
      sub_4d_v64_v64%w = a%w - b%w
   end function sub_4d_v64_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function neg_4d_v32(v)
      type(pikg_f32vec4), intent(in) :: v
      neg_4d_v32%x = - v%x
      neg_4d_v32%y = - v%y
      neg_4d_v32%z = - v%z
      neg_4d_v32%w = - v%w
   end function neg_4d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function neg_4d_v64(v)
      type(pikg_f64vec4), intent(in) :: v
      neg_4d_v64%x = - v%x
      neg_4d_v64%y = - v%y
      neg_4d_v64%z = - v%z
      neg_4d_v64%w = - v%w
   end function neg_4d_v64

   !####################################################################
   !========
   !   2D
   !========
   elemental pure type(pikg_f32vec2) function mul_2d_v32_i(v,i)
      type(pikg_f32vec2), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      mul_2d_v32_i%x = v%x * i
      mul_2d_v32_i%y = v%y * i
   end function mul_2d_v32_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function mul_2d_v32_f(v,f)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      mul_2d_v32_f%x = v%x * f
      mul_2d_v32_f%y = v%y * f
   end function mul_2d_v32_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function mul_2d_v32_d(v,d)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      mul_2d_v32_d%x = v%x * d
      mul_2d_v32_d%y = v%y * d
   end function mul_2d_v32_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function mul_2d_i_v32(i,v)
      integer(kind=c_int), intent(in) :: i
      type(pikg_f32vec2), intent(in) :: v
      mul_2d_i_v32%x = i * v%x 
      mul_2d_i_v32%y = i * v%y 
   end function mul_2d_i_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function mul_2d_f_v32(f,v)
      real(kind=c_float), intent(in) :: f
      type(pikg_f32vec2), intent(in) :: v
      mul_2d_f_v32%x = f * v%x 
      mul_2d_f_v32%y = f * v%y 
   end function mul_2d_f_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function mul_2d_d_v32(d,v)
      real(kind=c_double), intent(in) :: d
      type(pikg_f32vec2), intent(in) :: v
      mul_2d_d_v32%x = d * v%x 
      mul_2d_d_v32%y = d * v%y 
   end function mul_2d_d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function mul_2d_v64_i(v,i)
      type(pikg_f64vec2), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      mul_2d_v64_i%x = v%x * i
      mul_2d_v64_i%y = v%y * i
   end function mul_2d_v64_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function mul_2d_v64_f(v,f)
      type(pikg_f64vec2), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      mul_2d_v64_f%x = v%x * f
      mul_2d_v64_f%y = v%y * f
   end function mul_2d_v64_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function mul_2d_v64_d(v,d)
      type(pikg_f64vec2), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      mul_2d_v64_d%x = v%x * d
      mul_2d_v64_d%y = v%y * d
   end function mul_2d_v64_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function mul_2d_i_v64(i,v)
      integer(kind=c_int), intent(in) :: i
      type(pikg_f64vec2), intent(in) :: v
      mul_2d_i_v64%x = i * v%x 
      mul_2d_i_v64%y = i * v%y 
   end function mul_2d_i_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function mul_2d_f_v64(f,v)
      real(kind=c_float), intent(in) :: f
      type(pikg_f64vec2), intent(in) :: v
      mul_2d_f_v64%x = f * v%x 
      mul_2d_f_v64%y = f * v%y 
   end function mul_2d_f_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function mul_2d_d_v64(d,v)
      real(kind=c_double), intent(in) :: d
      type(pikg_f64vec2), intent(in) :: v
      mul_2d_d_v64%x = d * v%x 
      mul_2d_d_v64%y = d * v%y 
   end function mul_2d_d_v64

   !--------------------------------------------------------------------
   pure real(kind=c_float) function mul_2d_v32_fa(v,f)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      mul_2d_v32_fa = v%x * f(1) + v%y * f(2)
   end function mul_2d_v32_fa

   !--------------------------------------------------------------------
   pure real(kind=c_float) function mul_2d_fa_v32(f,v)
      real(kind=c_float), dimension(:), intent(in) :: f
      type(pikg_f32vec2), intent(in) :: v
      mul_2d_fa_v32 = f(1) * v%x + f(2) * v%y
   end function mul_2d_fa_v32

   !--------------------------------------------------------------------
   elemental pure real(kind=c_float) function mul_2d_v32_v32(a,b)
      type(pikg_f32vec2), intent(in) :: a,b
      mul_2d_v32_v32 = a%x * b%x + a%y * b%y
   end function mul_2d_v32_v32

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_2d_v32_da(v,d)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      mul_2d_v32_da = v%x * d(1) + v%y * d(2)
   end function mul_2d_v32_da

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_2d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec2), intent(in) :: v
      mul_2d_da_v32 = d(1) * v%x + d(2) * v%y
   end function mul_2d_da_v32

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_2d_v64_da(v,d)
      type(pikg_f64vec2), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      mul_2d_v64_da = v%x * d(1) + v%y * d(2)
   end function mul_2d_v64_da

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_2d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec2), intent(in) :: v
      mul_2d_da_v64 = d(1) * v%x + d(2) * v%y
   end function mul_2d_da_v64

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_2d_v64_v32(a,b)
      type(pikg_f64vec2), intent(in) :: a
      type(pikg_f32vec2), intent(in) :: b
      mul_2d_v64_v32 = a%x * b%x + a%y * b%y 
   end function mul_2d_v64_v32

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_2d_v32_v64(a,b)
      type(pikg_f32vec2), intent(in) :: a
      type(pikg_f64vec2), intent(in) :: b
      mul_2d_v32_v64 = a%x * b%x + a%y * b%y
   end function mul_2d_v32_v64

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_2d_v64_v64(a,b)
      type(pikg_f64vec2), intent(in) :: a,b
      mul_2d_v64_v64 = a%x * b%x + a%y * b%y 
   end function mul_2d_v64_v64
   !========
   !   3D
   !========
   elemental pure type(pikg_f32vec3) function mul_3d_v32_i(v,i)
      type(pikg_f32vec3), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      mul_3d_v32_i%x = v%x * i
      mul_3d_v32_i%y = v%y * i
      mul_3d_v32_i%z = v%z * i
   end function mul_3d_v32_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function mul_3d_v32_f(v,f)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      mul_3d_v32_f%x = v%x * f
      mul_3d_v32_f%y = v%y * f
      mul_3d_v32_f%z = v%z * f
   end function mul_3d_v32_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function mul_3d_v32_d(v,d)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      mul_3d_v32_d%x = v%x * d
      mul_3d_v32_d%y = v%y * d
      mul_3d_v32_d%z = v%z * d
   end function mul_3d_v32_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function mul_3d_i_v32(i,v)
      integer(kind=c_int), intent(in) :: i
      type(pikg_f32vec3), intent(in) :: v
      mul_3d_i_v32%x = i * v%x 
      mul_3d_i_v32%y = i * v%y 
      mul_3d_i_v32%z = i * v%z 
   end function mul_3d_i_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function mul_3d_f_v32(f,v)
      real(kind=c_float), intent(in) :: f
      type(pikg_f32vec3), intent(in) :: v
      mul_3d_f_v32%x = f * v%x 
      mul_3d_f_v32%y = f * v%y 
      mul_3d_f_v32%z = f * v%z 
   end function mul_3d_f_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function mul_3d_d_v32(d,v)
      real(kind=c_double), intent(in) :: d
      type(pikg_f32vec3), intent(in) :: v
      mul_3d_d_v32%x = d * v%x 
      mul_3d_d_v32%y = d * v%y 
      mul_3d_d_v32%z = d * v%z 
   end function mul_3d_d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function mul_3d_v64_i(v,i)
      type(pikg_f64vec3), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      mul_3d_v64_i%x = v%x * i
      mul_3d_v64_i%y = v%y * i
      mul_3d_v64_i%z = v%z * i
   end function mul_3d_v64_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function mul_3d_v64_f(v,f)
      type(pikg_f64vec3), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      mul_3d_v64_f%x = v%x * f
      mul_3d_v64_f%y = v%y * f
      mul_3d_v64_f%z = v%z * f
   end function mul_3d_v64_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function mul_3d_v64_d(v,d)
      type(pikg_f64vec3), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      mul_3d_v64_d%x = v%x * d
      mul_3d_v64_d%y = v%y * d
      mul_3d_v64_d%z = v%z * d
   end function mul_3d_v64_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function mul_3d_i_v64(i,v)
      integer(kind=c_int), intent(in) :: i
      type(pikg_f64vec3), intent(in) :: v
      mul_3d_i_v64%x = i * v%x 
      mul_3d_i_v64%y = i * v%y 
      mul_3d_i_v64%z = i * v%z 
   end function mul_3d_i_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function mul_3d_f_v64(f,v)
      real(kind=c_float), intent(in) :: f
      type(pikg_f64vec3), intent(in) :: v
      mul_3d_f_v64%x = f * v%x 
      mul_3d_f_v64%y = f * v%y 
      mul_3d_f_v64%z = f * v%z 
   end function mul_3d_f_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function mul_3d_d_v64(d,v)
      real(kind=c_double), intent(in) :: d
      type(pikg_f64vec3), intent(in) :: v
      mul_3d_d_v64%x = d * v%x 
      mul_3d_d_v64%y = d * v%y 
      mul_3d_d_v64%z = d * v%z 
   end function mul_3d_d_v64

   !--------------------------------------------------------------------
   pure real(kind=c_float) function mul_3d_v32_fa(v,f)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      mul_3d_v32_fa = v%x * f(1) + v%y * f(2) + v%z * f(3)
   end function mul_3d_v32_fa

   !--------------------------------------------------------------------
   pure real(kind=c_float) function mul_3d_fa_v32(f,v)
      real(kind=c_float), dimension(:), intent(in) :: f
      type(pikg_f32vec3), intent(in) :: v
      mul_3d_fa_v32 = f(1) * v%x  + f(2) * v%y + f(3) * v%z
   end function mul_3d_fa_v32

   !--------------------------------------------------------------------
   elemental pure real(kind=c_float) function mul_3d_v32_v32(a,b)
      type(pikg_f32vec3), intent(in) :: a,b
      mul_3d_v32_v32 = a%x * b%x + a%y * b%y + a%z * b%z 
   end function mul_3d_v32_v32

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_3d_v32_da(v,d)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      mul_3d_v32_da = v%x * d(1) + v%y * d(2) + v%z * d(3)
   end function mul_3d_v32_da

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_3d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec3), intent(in) :: v
      mul_3d_da_v32 = d(1) * v%x + d(2) * v%y + d(3) * v%z
   end function mul_3d_da_v32

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_3d_v64_da(v,d)
      type(pikg_f64vec3), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      mul_3d_v64_da = v%x * d(1) + v%y * d(2) + v%z * d(3)
   end function mul_3d_v64_da

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_3d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec3), intent(in) :: v
      mul_3d_da_v64 = d(1) * v%x + d(2) * v%y + d(3) * v%z
   end function mul_3d_da_v64

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_3d_v64_v32(a,b)
      type(pikg_f64vec3), intent(in) :: a
      type(pikg_f32vec3), intent(in) :: b
      mul_3d_v64_v32 = a%x * b%x + a%y * b%y + a%z * b%z 
   end function mul_3d_v64_v32

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_3d_v32_v64(a,b)
      type(pikg_f32vec3), intent(in) :: a
      type(pikg_f64vec3), intent(in) :: b
      mul_3d_v32_v64 = a%x * b%x + a%y * b%y + a%z * b%z 
   end function mul_3d_v32_v64

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_3d_v64_v64(a,b)
      type(pikg_f64vec3), intent(in) :: a,b
      mul_3d_v64_v64 = a%x * b%x + a%y * b%y + a%z * b%z 
   end function mul_3d_v64_v64
   !========
   !   4D
   !========
   elemental pure type(pikg_f32vec4) function mul_4d_v32_i(v,i)
      type(pikg_f32vec4), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      mul_4d_v32_i%x = v%x * i
      mul_4d_v32_i%y = v%y * i
      mul_4d_v32_i%z = v%z * i
      mul_4d_v32_i%w = v%w * i
   end function mul_4d_v32_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function mul_4d_v32_f(v,f)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      mul_4d_v32_f%x = v%x * f
      mul_4d_v32_f%y = v%y * f
      mul_4d_v32_f%z = v%z * f
      mul_4d_v32_f%w = v%w * f
   end function mul_4d_v32_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function mul_4d_v32_d(v,d)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      mul_4d_v32_d%x = v%x * d
      mul_4d_v32_d%y = v%y * d
      mul_4d_v32_d%z = v%z * d
      mul_4d_v32_d%w = v%w * d
   end function mul_4d_v32_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function mul_4d_i_v32(i,v)
      integer(kind=c_int), intent(in) :: i
      type(pikg_f32vec4), intent(in) :: v
      mul_4d_i_v32%x = i * v%x 
      mul_4d_i_v32%y = i * v%y 
      mul_4d_i_v32%z = i * v%z 
      mul_4d_i_v32%w = i * v%w
   end function mul_4d_i_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function mul_4d_f_v32(f,v)
      real(kind=c_float), intent(in) :: f
      type(pikg_f32vec4), intent(in) :: v
      mul_4d_f_v32%x = f * v%x 
      mul_4d_f_v32%y = f * v%y 
      mul_4d_f_v32%z = f * v%z 
      mul_4d_f_v32%w = f * v%w
   end function mul_4d_f_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function mul_4d_d_v32(d,v)
      real(kind=c_double), intent(in) :: d
      type(pikg_f32vec4), intent(in) :: v
      mul_4d_d_v32%x = d * v%x 
      mul_4d_d_v32%y = d * v%y 
      mul_4d_d_v32%z = d * v%z 
      mul_4d_d_v32%w = d * v%w
   end function mul_4d_d_v32

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function mul_4d_v64_i(v,i)
      type(pikg_f64vec4), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      mul_4d_v64_i%x = v%x * i
      mul_4d_v64_i%y = v%y * i
      mul_4d_v64_i%z = v%z * i
      mul_4d_v64_i%w = v%w * i
   end function mul_4d_v64_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function mul_4d_v64_f(v,f)
      type(pikg_f64vec4), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      mul_4d_v64_f%x = v%x * f
      mul_4d_v64_f%y = v%y * f
      mul_4d_v64_f%z = v%z * f
      mul_4d_v64_f%w = v%w * f
   end function mul_4d_v64_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function mul_4d_v64_d(v,d)
      type(pikg_f64vec4), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      mul_4d_v64_d%x = v%x * d
      mul_4d_v64_d%y = v%y * d
      mul_4d_v64_d%z = v%z * d
      mul_4d_v64_d%w = v%w * d
   end function mul_4d_v64_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function mul_4d_i_v64(i,v)
      integer(kind=c_int), intent(in) :: i
      type(pikg_f64vec4), intent(in) :: v
      mul_4d_i_v64%x = i * v%x 
      mul_4d_i_v64%y = i * v%y 
      mul_4d_i_v64%z = i * v%z 
      mul_4d_i_v64%w = i * v%w
   end function mul_4d_i_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function mul_4d_f_v64(f,v)
      real(kind=c_float), intent(in) :: f
      type(pikg_f64vec4), intent(in) :: v
      mul_4d_f_v64%x = f * v%x 
      mul_4d_f_v64%y = f * v%y 
      mul_4d_f_v64%z = f * v%z 
      mul_4d_f_v64%w = f * v%w
   end function mul_4d_f_v64

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function mul_4d_d_v64(d,v)
      real(kind=c_double), intent(in) :: d
      type(pikg_f64vec4), intent(in) :: v
      mul_4d_d_v64%x = d * v%x 
      mul_4d_d_v64%y = d * v%y 
      mul_4d_d_v64%z = d * v%z 
      mul_4d_d_v64%w = d * v%w
   end function mul_4d_d_v64

   !--------------------------------------------------------------------
   pure real(kind=c_float) function mul_4d_v32_fa(v,f)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_float), dimension(:), intent(in) :: f
      mul_4d_v32_fa = v%x * f(1) + v%y * f(2) + v%z * f(3) + v%w * f(4)
   end function mul_4d_v32_fa

   !--------------------------------------------------------------------
   pure real(kind=c_float) function mul_4d_fa_v32(f,v)
      real(kind=c_float), dimension(:), intent(in) :: f
      type(pikg_f32vec4), intent(in) :: v
      mul_4d_fa_v32 = f(1) * v%x  + f(2) * v%y + f(3) * v%z + f(4) * v%w
   end function mul_4d_fa_v32

   !--------------------------------------------------------------------
   elemental pure real(kind=c_float) function mul_4d_v32_v32(a,b)
      type(pikg_f32vec4), intent(in) :: a,b
      mul_4d_v32_v32 = a%x * b%x + a%y * b%y + a%z * b%z + a%w * b%w
   end function mul_4d_v32_v32

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_4d_v32_da(v,d)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      mul_4d_v32_da = v%x * d(1) + v%y * d(2) + v%z * d(3) + v%w * d(4)
   end function mul_4d_v32_da

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_4d_da_v32(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f32vec4), intent(in) :: v
      mul_4d_da_v32 = d(1) * v%x + d(2) * v%y + d(3) * v%z + d(4) * v%w
   end function mul_4d_da_v32

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_4d_v64_da(v,d)
      type(pikg_f64vec4), intent(in) :: v
      real(kind=c_double), dimension(:), intent(in) :: d
      mul_4d_v64_da = v%x * d(1) + v%y * d(2) + v%z * d(3) + v%w * d(4)
   end function mul_4d_v64_da

   !--------------------------------------------------------------------
   pure real(kind=c_double) function mul_4d_da_v64(d,v)
      real(kind=c_double), dimension(:), intent(in) :: d
      type(pikg_f64vec4), intent(in) :: v
      mul_4d_da_v64 = d(1) * v%x + d(2) * v%y + d(3) * v%z + d(4) * v%w
   end function mul_4d_da_v64

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_4d_v64_v32(a,b)
      type(pikg_f64vec4), intent(in) :: a
      type(pikg_f32vec4), intent(in) :: b
      mul_4d_v64_v32 = a%x * b%x + a%y * b%y + a%z * b%z + a%w * b%w
   end function mul_4d_v64_v32

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_4d_v32_v64(a,b)
      type(pikg_f32vec4), intent(in) :: a
      type(pikg_f64vec4), intent(in) :: b
      mul_4d_v32_v64 = a%x * b%x + a%y * b%y + a%z * b%z + a%w * b%w
   end function mul_4d_v32_v64

   !--------------------------------------------------------------------
   elemental pure real(kind=c_double) function mul_4d_v64_v64(a,b)
      type(pikg_f64vec4), intent(in) :: a,b
      mul_4d_v64_v64 = a%x * b%x + a%y * b%y + a%z * b%z + a%w * b%w
   end function mul_4d_v64_v64

   !####################################################################
   !========
   !   2D
   !========
   elemental pure type(pikg_f32vec2) function div_2d_v32_i(v,i)
      type(pikg_f32vec2), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      div_2d_v32_i%x = v%x / i
      div_2d_v32_i%y = v%y / i
   end function div_2d_v32_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function div_2d_v32_f(v,f)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      div_2d_v32_f%x = v%x / f
      div_2d_v32_f%y = v%y / f
   end function div_2d_v32_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec2) function div_2d_v32_d(v,d)
      type(pikg_f32vec2), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      div_2d_v32_d%x = v%x / d
      div_2d_v32_d%y = v%y / d
   end function div_2d_v32_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function div_2d_v64_i(v,i)
      type(pikg_f64vec2), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      div_2d_v64_i%x = v%x / i
      div_2d_v64_i%y = v%y / i
   end function div_2d_v64_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function div_2d_v64_f(v,f)
      type(pikg_f64vec2), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      div_2d_v64_f%x = v%x / f
      div_2d_v64_f%y = v%y / f
   end function div_2d_v64_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec2) function div_2d_v64_d(v,d)
      type(pikg_f64vec2), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      div_2d_v64_d%x = v%x / d
      div_2d_v64_d%y = v%y / d
   end function div_2d_v64_d
   !========
   !   3D
   !========
   elemental pure type(pikg_f32vec3) function div_3d_v32_i(v,i)
      type(pikg_f32vec3), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      div_3d_v32_i%x = v%x / i
      div_3d_v32_i%y = v%y / i
      div_3d_v32_i%z = v%z / i
   end function div_3d_v32_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function div_3d_v32_f(v,f)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      div_3d_v32_f%x = v%x / f
      div_3d_v32_f%y = v%y / f
      div_3d_v32_f%z = v%z / f
   end function div_3d_v32_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec3) function div_3d_v32_d(v,d)
      type(pikg_f32vec3), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      div_3d_v32_d%x = v%x / d
      div_3d_v32_d%y = v%y / d
      div_3d_v32_d%z = v%z / d
   end function div_3d_v32_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function div_3d_v64_i(v,i)
      type(pikg_f64vec3), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      div_3d_v64_i%x = v%x / i
      div_3d_v64_i%y = v%y / i
      div_3d_v64_i%z = v%z / i
   end function div_3d_v64_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function div_3d_v64_f(v,f)
      type(pikg_f64vec3), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      div_3d_v64_f%x = v%x / f
      div_3d_v64_f%y = v%y / f
      div_3d_v64_f%z = v%z / f
   end function div_3d_v64_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec3) function div_3d_v64_d(v,d)
      type(pikg_f64vec3), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      div_3d_v64_d%x = v%x / d
      div_3d_v64_d%y = v%y / d
      div_3d_v64_d%z = v%z / d
   end function div_3d_v64_d
   !========
   !   4D
   !========
   elemental pure type(pikg_f32vec4) function div_4d_v32_i(v,i)
      type(pikg_f32vec4), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      div_4d_v32_i%x = v%x / i
      div_4d_v32_i%y = v%y / i
      div_4d_v32_i%z = v%z / i
      div_4d_v32_i%w = v%w / i
   end function div_4d_v32_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function div_4d_v32_f(v,f)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      div_4d_v32_f%x = v%x / f
      div_4d_v32_f%y = v%y / f
      div_4d_v32_f%z = v%z / f
      div_4d_v32_f%w = v%w / f
   end function div_4d_v32_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f32vec4) function div_4d_v32_d(v,d)
      type(pikg_f32vec4), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      div_4d_v32_d%x = v%x / d
      div_4d_v32_d%y = v%y / d
      div_4d_v32_d%z = v%z / d
      div_4d_v32_d%w = v%w / d
   end function div_4d_v32_d

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function div_4d_v64_i(v,i)
      type(pikg_f64vec4), intent(in) :: v
      integer(kind=c_int), intent(in) :: i
      div_4d_v64_i%x = v%x / i
      div_4d_v64_i%y = v%y / i
      div_4d_v64_i%z = v%z / i
      div_4d_v64_i%w = v%w / i
   end function div_4d_v64_i

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function div_4d_v64_f(v,f)
      type(pikg_f64vec4), intent(in) :: v
      real(kind=c_float), intent(in) :: f
      div_4d_v64_f%x = v%x / f
      div_4d_v64_f%y = v%y / f
      div_4d_v64_f%z = v%z / f
      div_4d_v64_f%w = v%w / f
   end function div_4d_v64_f

   !--------------------------------------------------------------------
   elemental pure type(pikg_f64vec4) function div_4d_v64_d(v,d)
      type(pikg_f64vec4), intent(in) :: v
      real(kind=c_double), intent(in) :: d
      div_4d_v64_d%x = v%x / d
      div_4d_v64_d%y = v%y / d
      div_4d_v64_d%z = v%z / d
      div_4d_v64_d%w = v%w / d
   end function div_4d_v64_d

end module pikg_vector
