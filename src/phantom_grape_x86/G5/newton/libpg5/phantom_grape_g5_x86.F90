!================================
! MODULE: phantom_grape_g5_x86
!================================
module phantom_grape_g5_x86
   use, intrinsic :: iso_c_binding

   !* Definitions of interfaces to C functions
   interface

      function g5_get_number_of_pipelines() &
            bind(c,name='g5_get_number_of_pipelines')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: g5_get_number_of_pipelines
      end function g5_get_number_of_pipelines

      function g5_get_jmemsize() bind(c,name='g5_get_jmemsize')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: g5_get_jmemsize
      end function g5_get_jmemsize

      subroutine g5_open() bind(c,name='g5_open')
         use, intrinsic :: iso_c_binding
         implicit none
      end subroutine g5_open

      subroutine g5_close() bind(c,name='g5_close')
         use, intrinsic :: iso_c_binding
         implicit none
      end subroutine g5_close

      subroutine g5_set_eps_to_all(eps) bind(c,name='g5_set_eps_to_all')
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double), value, intent(in) :: eps
      end subroutine g5_set_eps_to_all

      subroutine g5_set_range(xmin,xmax,mmin) bind(c,name='g5_set_range')
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double), value, intent(in) :: xmin,xmax,mmin
      end subroutine g5_set_range

      subroutine g5_set_n(n) bind(c,name='g5_set_n')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: n
      end subroutine g5_set_n

      subroutine g5_set_nMC(devid,n) bind(c,name='g5_set_nMC')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: devid,n
      end subroutine g5_set_nMC

      subroutine g5_set_xi(ni,xi) bind(c,name='g5_set_xi')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: ni
         real(kind=c_double), dimension(3,ni) :: xi
      end subroutine g5_set_xi

      subroutine g5_set_xiMC(devid,ni,xi) bind(c,name='g5_set_xiMC')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: devid,ni
         real(kind=c_double), dimension(3,ni) :: xi
      end subroutine g5_set_xiMC

      subroutine g5_set_xmj(adr,nj,xj,mj) bind(c,name='g5_set_xmj')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: adr,nj
         real(kind=c_double), dimension(3,nj) :: xj
         real(kind=c_double), dimension(nj) :: mj
      end subroutine g5_set_xmj

      subroutine g5_set_xmjMC(devid,adr,nj,xj,mj) &
            bind(c,name='g5_set_xmjMC')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: devid,adr,nj
         real(kind=c_double), dimension(3,nj) :: xj
         real(kind=c_double), dimension(nj) :: mj
      end subroutine g5_set_xmjMC

      subroutine g5_run() bind(c,name='g5_run')
         use, intrinsic :: iso_c_binding
         implicit none
      end subroutine g5_run

      subroutine g5_runMC(devid) bind(c,name='g5_runMC')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: devid 
      end subroutine g5_runMC

      subroutine g5_get_force(ni,a,p) bind(c,name='g5_get_force')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: ni
         real(kind=c_double), dimension(3,ni) :: a
         real(kind=c_double), dimension(ni) :: p
      end subroutine g5_get_force

      subroutine g5_get_forceMC(devid,ni,a,p) &
            bind(c,name='g5_get_forceMC')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: devid,ni
         real(kind=c_double), dimension(3,ni) :: a
         real(kind=c_double), dimension(ni) :: p
      end subroutine g5_get_forceMC

      subroutine g5_calculate_force_on_x(x,a,p,ni) &
            bind(c,name='g5_calculate_force_on_x')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: ni
         real(kind=c_double), dimension(3,ni) :: x,a 
         real(kind=c_double), dimension(ni) :: p
      end subroutine g5_calculate_force_on_x

      subroutine g5_calculate_force_on_xMC(devid,x,a,p,ni) &
            bind(c,name='g5_calculate_force_on_xMC')
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: devid
         integer(kind=c_int), value, intent(in) :: ni
         real(kind=c_double), dimension(3,ni) :: x,a 
         real(kind=c_double), dimension(ni) :: p
      end subroutine g5_calculate_force_on_xMC

   end interface

end module phantom_grape_g5_x86


