!===============================
!   MODULE: User defined types
!===============================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   use fdps_super_particle
   implicit none
  
   !**** Full particle type
   type, public, bind(c) :: fplj !$fdps FP,EPI,EPJ,Force
       !$fdps copyFromForce fplj (pot,pot) (acc,acc)
       !$fdps copyFromFP fplj (id,id) (mass,mass) (pos,pos) (search_radius,search_radius)
       !$fdps clear id=keep, mass=keep, pos=keep, vel=keep, search_radius=keep
       integer(kind=c_long_long) :: id
       real(kind=c_double) :: mass !$fdps charge
       type(fdps_f64vec) :: pos !$fdps position
       type(fdps_f64vec) :: vel !$fdps velocity
       type(fdps_f64vec) :: acc
       real(kind=c_double) :: pot
       real(kind=c_double) :: search_radius !$fdps rsearch
   end type fplj

   contains

   !**** Interaction function 
   subroutine calc_force_fpfp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(fplj), dimension(n_ip), intent(in) :: ep_i
      type(fplj), dimension(n_jp), intent(in) :: ep_j
      type(fplj), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      integer(c_long_long) :: idi
      real(c_double) :: r0,r0sq,r0inv,r0invp6,r0invp7,foffset,poffset
      real(c_double) :: poti,r2,r,r_inv,r2_inv,r6_inv,r12_inv,coef
      type(fdps_f64vec) :: xi,ai,rij

      r0 = 3.0d0
      r0sq = r0*r0
      r0inv = 1.0d0/r0
      r0invp6 = 1.0d0/(r0sq*r0sq*r0sq)
      r0invp7 = r0invp6*r0inv
      foffset = -12.0d0*r0invp6*r0invp7 + 6.0d0*r0invp7
      poffset = -13.0d0*r0invp6*r0invp6 + 7.0d0*r0invp6

      do i=1,n_ip
         idi = ep_i(i)%id
         xi = ep_i(i)%pos
         ai = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            if (idi == ep_j(j)%id) cycle
            rij%x = xi%x - ep_j(j)%pos%x
            rij%y = xi%y - ep_j(j)%pos%y
            rij%z = xi%z - ep_j(j)%pos%z
            r2 = rij%x * rij%x &
               + rij%y * rij%y &
               + rij%z * rij%z
            if (r2 < r0sq) then
               r_inv = 1.0d0/sqrt(r2)
               r = r2*r_inv
               r2_inv = r_inv * r_inv
               r6_inv = r2_inv * r2_inv * r2_inv
               r12_inv = r6_inv * r6_inv
               poti = poti + (r12_inv - r6_inv - foffset*r + poffset)
               coef = (12.0d0*r12_inv*r2_inv - 6.0d0*r6_inv*r2_inv + foffset*r_inv)
               ai%x = ai%x + coef * rij%x
               ai%y = ai%y + coef * rij%y
               ai%z = ai%z + coef * rij%z
            end if
         end do 
         f(i)%acc = f(i)%acc + ai
         f(i)%pot = f(i)%pot + poti
      end do

   end subroutine calc_force_fpfp

end module user_defined_types
