C> @file eqnsolver_cmt.f Routines for entire terms on RHS. Mostly volume integrals

C> \ingroup diffhvol
C> @{
C> Volume integral for diffusive terms. Compute \f$\mathbf{H}^d\f$
C> and store it for one element. Store faces of \f$\mathbf{H}^d\f$ for IGU. 
      subroutine viscous_cmt(e,eq)
      include  'SIZE'
      include  'CMTDATA'
      include  'DG'
      include  'INPUT'

      common /CMTSURFLX/ fatface(heresize),graduf(hdsize)
      real fatface,graduf

      integer e,eq

      if (eq .lt. toteq) then ! not energy
         if (eq .gt. ldim+1) return ! not if3d
      endif

      nxyz=lx1*ly1*lz1
      nfq=lx1*lz1*2*ldim*nelt
      nstate = nqq
! where different things live
      iqm =1
      iqp =iqm+nstate*nfq
      iuj =iqp+nstate*nfq

      call rzero(diffh,3*nxyz)

      call agradu(diffh,gradu,e,eq)

! FIND OUT IF YOU NEED TO COMPUTE {{Hd.n}}
!     call diffh2gradu(e,eq,graduf)
! OR IF YOU NEED TO COMPUTE {{Hd}}.n
      call diffh2face(e,eq,graduf)  ! on faces {{Hd}} after element loop

! volume integral involving "DG-ish" stiffness operator K
      call half_iku_cmt(res1(1,1,1,e,eq),diffh,e)

C> @}
      return
      end

!-----------------------------------------------------------------------

C> \ingroup vsurf
C> @{
C> \f$G^T U\f$
      subroutine igtu_cmt(qminus,ummcu,hface)

!     Vol integral [[u]].{{gradv}}. adapted from Lu's dgf3.f;
! interior penalty stuff could go
! here too. Debug subroutine ihu in heat.usr and try to fold it in here.

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'DG'      ! iface
      include 'CMTDATA'
      include 'SOLN' ! for vz. goes away when agradu_ns works

! arguments
      real qminus(lx1*lz1,2*ldim,nelt,*)    ! intent(in)
      real ummcu(lx1*lz1*2*ldim,nelt,toteq) ! intent(in)
      real hface(lx1*lz1*2*ldim*nelt,toteq,3) ! intent(out) scratch

! commons and scratch
      common /scrns/ superhugeh(lx1*ly1*lz1*lelt,3) ! like totalh, but super-huge
      common /scruz/ gradm1_t_overwrites(lx1*ly1*lz1*lelt) ! sigh
      real superhugeh,gradm1_t_overwrites
!     common /ctmp0/ viscscr(lx1,ly1,lz1)
!     real viscscr
!     parameter (lfq=lx1*lz1*2*ldim*lelt)
!     common /ctmp0/ ftmp1(lfq),ftmp2(lfq)
!     real ftmp1,ftmp2

      integer eijk3(3,3,3)
!     data eijk2 / 0, -1, 1, 0/
      data eijk3
     >/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/

      integer e, eq, n, npl, nf, f, i, k, eq2

      nxz    = lx1*lz1
      nfaces = 2*ldim
      nf     = nxz*nfaces ! 1 element's face points
      nfq    = nf*nelt ! all points in a pile of faces
!     if (ifsip) then
!        const=-1.0 ! SIP
!     else
         const=1.0 ! Baumann-Oden
!     endif

! compute (U-{{U}})_i * n_k
! OK FOLKS GIANT BUG UMMCU IS BAD AT INFLOW
      l=1
      do e=1,nelt
         do eq=1,toteq
            call col3(hface(l,eq,3),ummcu(1,e,eq),area(1,1,1,e),nf)
            call col3(hface(l,eq,1),hface(l,eq,3),unx(1,1,1,e), nf)
            call col3(hface(l,eq,2),hface(l,eq,3),uny(1,1,1,e), nf)
            if(if3d) call col2(hface(l,eq,3),unz(1,1,1,e),nf)
         enddo
         l=l+nf
      enddo

      nxyz  =lx1*ly1*lz1
      nvol  =nxyz*nelt
      ngradu=nxyz*toteq*3
      do eq=1,toteq
         call rzero(superhugeh,3*lx1*ly1*lz1*lelt)
         if (eq .eq. 4 .and. .not. if3d) goto 133
         l=1
         m=1
         do e=1,nelt
            call rzero(diffh,3*nxyz)
            call rzero(gradu,ngradu) ! this too goes away when gradu is global
            do j=1,ldim
               do eq2=1,toteq ! sigh
                  call add_face2full_cmt(1,lx1,ly1,lz1,iface_flux(1,e),
     >                                gradu(1,eq2,j),hface(l,eq2,j))
               enddo
            enddo

            l=l+nf ! index for hface, which is global. this all goes away
                ! once you get rid of that execrable "element loop" in
                ! compute_rhs_and_dt
!           call fluxj_ns(superhugeh,... THIS will be correctly strided as well
! JH110716 AND someday it will work
!!            do j=1,ldim    ! flux direction
!!               do k=1,ldim ! dU   direction
!!                  ieijk=0
!!                  if (eq .lt. toteq) ieijk=eijk3(eq-1,j,k) ! does this work in 2D?
!!                  if (ieijk .eq. 0) then
!!                     call agradu_ns(superhugeh(m,j),gradu(1,1,k),viscscr
!!     >                             ,e,eq,j,k) ! the worst stride ever
!!                  endif
!!               enddo
!!            enddo
! JH110716 but not today. for now, let maxima do my thinking in the fluxj* routines

            call agradu(diffh,gradu,e,eq)

            do j=1,ldim
               call copy(superhugeh(m,j),diffh(1,j),nxyz)
            enddo

            m=m+nxyz

         enddo ! element loop

! gradm1_t uses /ctmp1/
         call gradm1_t(gradm1_t_overwrites,superhugeh(1,1),
     >                        superhugeh(1,2),superhugeh(1,3))
         call cmult(gradm1_t_overwrites,const,nvol)
         call add2(res1(1,1,1,1,eq),gradm1_t_overwrites,nvol)
133      continue
      enddo ! equation loop

C> @}

      return
      end

!-----------------------------------------------------------------------

C> \ingroup convhvol
C> @{
C> Convective volume terms formed and differentiated^T here
      subroutine convective_cmt(e)
! JH081916 convective flux divergence integrated in weak form and
!          placed in res1.
! JH052418 convective flux divergence integrated in strong/ultraweak form and
!          placed in res1.
! JH060518 interchanged with equation loop (now inside here!)
      include 'SIZE'
      include 'CMTDATA'
      include 'GEOM'
      integer e,eq

      n=3*lxd*lyd*lzd*toteq

      call rzero(convh,n)

      if (lxd.gt.lx1) then

!Comment by BAD Feb 14. 2020: Conv_h storage has changed thus
!evaluate_dealiased_conv_h must also change. For now simple eq loop
!added. However probs should build it like Jason did in
!evaulate_aliased_conv_h.
         call evaluate_dealiased_conv_h(e,eq)
         call copy(totalh,convh,n)
         call fluxdiv_dealiased_weak_chain(e,eq)
      else

       call evaluate_aliased_conv_h(e)
       do eq=1,toteq
       call contravariant_flux(totalh(1,1,eq),convh(1,1,eq),rx(1,1,e),1)
       enddo

!     call fluxdiv_2point_scr(res1,totalh,e,rx(1,1,e))
!     call fluxdiv_2point_noscr(res1,totalh,e,rx(1,1,e))
! two-point, KEP/EC
       call fluxdiv_2point_noscr(res1,totalh,e,rx(1,1,e))
      endif  
      return
      end
C> @}
C> Convective pointwise flux function \f$\mathbf{H}^c\f$ on fine grid.
      subroutine evaluate_dealiased_conv_h(e,eq)
! computed as products between primitive variables and conserved variables.
! if you want to write rho u_i u_j as (rho u_i) (rho u_j) (rho^{-1}), this
! is the place to do it

!Comment by BAD Feb 6th 2020, This routine still assumes toteq is 5. ie 
!the energy equation is the 5th Will need to change. Alot of indicis 
!for convh are explicitly set with this in mind. This will also have to 
!change. The variable eq can be used to make the code more versitile 
!but the "if" statments will have to change and the add2col2 call that
!uses eq twice as an indici
      include  'SIZE'
      include  'SOLN'
      include  'DEALIAS'
      include  'CMTDATA'
      include  'INPUT'
      integer  e,eq

      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ JIu1(ldd),JIu2(ldd)!,ur(ldd),us(ldd),ud(ldd),tu(ldd)
      real JIu1,JIu2

      n=lxd*lyd*lzd


c To be consistent with momentum equation, for mass balance flux vector is 
c computed by multiplying rho by u_j
         call intp_rstd(JIu1,phig(1,1,1,e),lx1,lxd,if3d,0)
         call intp_rstd(JIu2,pr(1,1,1,e),lx1,lxd,if3d,0)
!Put velocity here!
         call intp_rstd(vxd(1,1,1,e),vx(1,1,1,e),lx1,lxd,if3d,0)
         call intp_rstd(vyd(1,1,1,e),vy(1,1,1,e),lx1,lxd,if3d,0)
         if (if3d) 
     >    call intp_rstd(vzd(1,1,1,e),vz(1,1,1,e),lx1,lxd,if3d,0)
!Mass
      do j=1,ldim  
         call intp_rstd(convh(1,j,1),u(1,1,1,j+1,e),lx1,lxd,if3d,0)
      enddo  
!MOM
      do eq=2,ldim+1
         call intp_rstd(convh(1,1,eq),u(1,1,1,eq,e),lx1,lxd,if3d,0)
         do j=2,ldim
            call copy(convh(1,j,eq),convh(1,1,eq),n)
         enddo
         call col2(convh(1,1,eq),vxd(1,1,1,e),n)
         call col2(convh(1,2,eq),vyd(1,1,1,e),n)
         if (if3d) call col2(convh(1,3,eq),vzd(1,1,1,e),n)
         call add2col2(convh(1,eq-1,eq),JIu1,JIu2,n)  
      enddo  
!ENE
      eq=toteq  
      call intp_rstd(convh(1,1,eq),u(1,1,1,eq,e),lx1,lxd,if3d,0)
      call add2col2(convh(1,1,eq),JIu1,JIu2,n)
      do j=2,ldim
         call copy(convh(1,j,eq),convh(1,1,eq),n)
      enddo
      call col2(convh(1,1,eq),vxd(1,1,1,e),n)
      call col2(convh(1,2,eq),vyd(1,1,1,e),n)
      if (if3d) call col2(convh(1,3,eq),vzd(1,1,1,e),n)

      return 
      end    

!Comment by BAD Feb. 19 2020. Copy of old code for ref, will delete later
!      if (eq .eq. 1) then ! convective flux of mass=rho u_j=U_{j+1}
!
!         do j=1,ldim
!            call intp_rstd(convh(1,j,eq),u(1,1,1,eq+j,e),lx1,lxd,if3d,0)
!         enddo
!
!      else
!
!c To be consistent with momentum equation, for mass balance flux vector is 
!c computed by multiplying rho by u_j
!         call intp_rstd(JIu1,phig(1,1,1,e),lx1,lxd,if3d,0)
!         call intp_rstd(JIu2,pr(1,1,1,e),lx1,lxd,if3d,0)
!
!         if (eq .lt. 5) then ! self-advection of rho u_i by rho u_i u_j
!
!            call intp_rstd(convh(1,1,eq),u(1,1,1,eq,e),lx1,lxd,if3d,0)
!            do j=2,ldim
!               call copy(convh(1,j,eq),convh(1,1,eq),n)
!            enddo
!            call col2(convh(1,1,eq),vxd(1,1,1,e),n)
!            call col2(convh(1,2,eq),vyd(1,1,1,e),n)
!            if (if3d) call col2(convh(1,3,eq),vzd(1,1,1,e),n)
!            call add2col2(convh(1,eq-1,eq),JIu1,JIu2,n)
!
!         elseif (eq .eq. 5) then
!
!            call intp_rstd(convh(1,1,eq),u(1,1,1,eq,e),lx1,lxd,if3d,0)
!            call add2col2(convh(1,1,eq),JIu1,JIu2,n)
!            do j=2,ldim
!               call copy(convh(1,j,eq),convh(1,1,eq),n)
!            enddo
!            call col2(convh(1,1,eq),vxd(1,1,1,e),n)
!            call col2(convh(1,2,eq),vyd(1,1,1,e),n)
!            call col2(convh(1,3,eq),vzd(1,1,1,e),n)
!
!         else
!            if(nio.eq.0) write(6,*) 'eq=',eq,'really must be <= 5'
!            if(nio.eq.0) write(6,*) 'aborting in evaluate_conv_h'
!            call exitt
!         endif
!
!      endif
!     
!      return
!      end
!


      subroutine fluxdiv_2point_slow(res,e,ja)
      include 'SIZE'
      include 'DG'
      include 'GEOM' ! diagnostic. conflicts with ja
      include 'SOLN'
      include 'CMTDATA'
! JH062018 two-point energy-preserving/SBP flux divergence volume integral, slow and naive
      integer e,eq
      real res(lx1,ly1,lz1,lelt,toteq) ! CMTDATA lurks
      real ja(lx1,ly1,lz1,ldim,ldim)   ! rst outermost
! scratch element for extra variables (hardcoded) and conserved variables U
! transposed to quantity-innermost. unvectorizable?
      common /scrns/ zaux (nparm,lx1,ly1,lz1),ut(toteq,lx1,ly1,lz1),
     >               zauxt(lx1*ly1*lz1,nparm),jat(3,3,lx1,ly1,lz1)
      real zaux,ut,zauxt,jat
      real flx(5)

      nxyz=lx1*ly1*lz1
      call cmt_usrz(zaux,zauxt,ut,e,nparm)

      call rzero(jat,9*lx1*ly1*lz1)
      do j=1,ldim
         do i=1,ldim
            do iz=1,lz1
            do iy=1,ly1
            do ix=1,lx1
               jat(i,j,ix,iy,iz)=ja(ix,iy,iz,i,j)
            enddo
            enddo
            enddo
         enddo
      enddo

      do iz=1,lz1
      do iy=1,ly1
      do ix=1,lx1
! r-direction
         do l=1,lx1
            call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,l,iy,iz),
     >                         zaux(1,ix,iy,iz),zaux(1,l,iy,iz),
     >                         jat(1,1,ix,iy,iz),jat(1,1,l,iy,iz))
            do eq=1,toteq
               res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+
     >                   dstrong(ix,l)*flx(eq)
            enddo
         enddo
! s-direction
         do l=1,ly1
            call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,ix,l,iz),
     >                         zaux(1,ix,iy,iz),zaux(1,ix,l,iz),
     >                         jat(1,2,ix,iy,iz),jat(1,2,ix,l,iz))
            do eq=1,toteq
               res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+
     >                   dstrong(iy,l)*flx(eq)
            enddo
         enddo
      enddo
      enddo
      enddo

      if (lz1.gt.1) then
      do iz=1,lz1
      do iy=1,ly1
      do ix=1,lx1
! t-direction
         do l=1,lz1
            call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,ix,iy,l),
     >                         zaux(1,ix,iy,iz),zaux(1,ix,iy,l),
     >                         jat(1,3,ix,iy,iz),jat(1,3,ix,iy,l))
            do eq=1,toteq
               res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+
     >                   dstrong(iz,l)*flx(eq)
            enddo
         enddo
      enddo
      enddo
      enddo
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine fluxdiv_2point_scr(res,fcons,e,ja)
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
! JH060418 set up two-point energy-preserving/SBP flux divergence volume integral
!          in the contravariant frame, call fluxdiv_2point, and increment res
!          for e^th element.
!          Metric terms Ja probably shouldn't be an argument but whatever.
!          cmt_usr2pt passes arguments to F# in the usr file
      integer e,eq
      real res(lx1,ly1,lz1,lelt,toteq) ! CMTDATA lurks
      real ja(lx1,ly1,lz1,ldim,ldim)   ! rst outermost
      real fcons(lx1,ly1,lz1,3,toteq)   ! consistent ``1-point'' flux
! scratch element for extra variables (hardcoded) and conserved variables U
! transposed to quantity-innermost. unvectorizable?
      common /scrns/ zaux (nparm,lx1,ly1,lz1),ut(toteq,lx1,ly1,lz1),
     >               zauxt(lx1*ly1*lz1,nparm),jat(3,3,lx1,ly1,lz1)
     >              ,rhsscr(lx1,toteq)
      real zaux,ut,zauxt,jat,rhsscr
      real flx(5)

      nxyz=lx1*ly1*lz1
      call cmt_usrz(zaux,zauxt,ut,e,nparm)

      call rzero(jat,9*lx1*ly1*lz1)
      do j=1,ldim
         do i=1,ldim
            do iz=1,lz1
            do iy=1,ly1
            do ix=1,lx1
               jat(i,j,ix,iy,iz)=ja(ix,iy,iz,i,j)
            enddo
            enddo
            enddo
         enddo
      enddo

! r-direction
      do iz=1,lz1
      do iy=1,ly1
         call rzero(rhsscr,toteq*lx1)
         do ix=1,lx1
            do l=ix+1,lx1
               call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,l,iy,iz),
     >                            zaux(1,ix,iy,iz),zaux(1,l,iy,iz),
     >                            jat(1,1,ix,iy,iz),jat(1,1,l,iy,iz))
               do eq=1,toteq
                  rhsscr(ix,eq)=rhsscr(ix,eq)+dstrong(ix,l)*flx(eq)
                  rhsscr(l,eq)=rhsscr(l,eq)+dstrong(l,ix)*flx(eq)
               enddo
            enddo
            do eq=1,toteq
               rhsscr(ix,eq)=rhsscr(ix,eq)+
     >                       dstrong(ix,ix)*fcons(ix,iy,iz,1,eq)
            enddo
         enddo
         do eq=1,toteq
            call add2(res(1,iy,iz,e,eq),rhsscr(1,eq),lx1)
         enddo
      enddo ! iy
      enddo ! iz


! consider repacking ut and zaux with iy in second place
! diagnostic kg is consistent in r and s
!            call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,ix,iy,iz),
!     >                            zaux(1,ix,iy,iz),zaux(1,ix,iy,iz),
!     >                            jat(1,2,ix,iy,iz),jat(1,2,ix,iy,iz))
!               write(60+eq,'(5e15.7)')
!     >xm1(ix,iy,iz,e),ym1(ix,iy,iz,e),flx(eq),fcons(ix,iy,iz,2,eq),
!     >abs(flx(eq)-fcons(ix,iy,iz,2,eq))
! diagnostic
! s-direction
      do iz=1,lz1
      do ix=1,lx1
         call rzero(rhsscr,toteq*ly1)
         do iy=1,ly1
            do l=iy+1,ly1
               call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,ix,l,iz),
     >                            zaux(1,ix,iy,iz),zaux(1,ix,l,iz),
     >                            jat(1,2,ix,iy,iz),jat(1,2,ix,l,iz))
               do eq=1,toteq
                  rhsscr(iy,eq)=rhsscr(iy,eq)+dstrong(iy,l)*flx(eq)
                  rhsscr(l,eq)=rhsscr(l,eq)+dstrong(l,iy)*flx(eq)
               enddo
            enddo
            do eq=1,toteq
               rhsscr(iy,eq)=rhsscr(iy,eq)+
     >                       dstrong(iy,iy)*fcons(ix,iy,iz,2,eq)
            enddo
         enddo
         do eq=1,toteq
            do iy=1,ly1
               res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+
     >                             rhsscr(iy,eq)
            enddo
         enddo
      enddo ! ix
      enddo ! iz

      if (lz1.gt.1) then
      write(6,*) 'duh sir t-direction'
! consider repacking ut and zaux with iz in second place
! t-direction
      do iy=1,ly1
      do ix=1,lx1
         call rzero(rhsscr,toteq*lz1)
         do iz=1,lz1
            do l=iz+1,lz1
               call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,ix,iy,l),
     >                            zaux(1,ix,iy,iz),zaux(1,ix,iy,l),
     >                            jat(1,3,ix,iy,iz),jat(1,3,ix,iy,l))
               do eq=1,toteq
                  rhsscr(iz,eq)=rhsscr(iz,eq)+dstrong(iz,l)*flx(eq)
                  rhsscr(l,eq)=rhsscr(l,eq)+dstrong(l,iz)*flx(eq)
               enddo
            enddo
            do eq=1,toteq
               rhsscr(iz,eq)=rhsscr(iz,eq)+
     >                       dstrong(iz,iz)*fcons(ix,iy,iz,3,eq)
            enddo
         enddo
         do eq=1,toteq
            do iz=1,lz1
               res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+
     >                             rhsscr(iz,eq)
            enddo
         enddo
      enddo ! ix
      enddo ! iy
      endif ! if3d

      return
      end

!-----------------------------------------------------------------------

      subroutine fluxdiv_2point_noscr(res,fcons,e,ja)
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
! JH061118 Access 3D res array directly without regard to stride. compare
!          performance to existing scr
! JH060418 set up two-point energy-preserving/SBP flux divergence volume integral
!          in the contravariant frame, call fluxdiv_2point, and increment res
!          for e^th element.
!          Metric terms Ja probably shouldn't be an argument but whatever.
!          cmt_usr2pt passes arguments to F# in the usr file
!          and FLUXO (github.com/project-fluxo/fluxo)
      integer e,eq
      real res(lx1,ly1,lz1,lelt,toteq) ! res1 lurks in CMTDATA
      real ja(lx1,ly1,lz1,ldim,ldim)   ! rst outermost
      real fcons(lx1,ly1,lz1,3,toteq)   ! consistent ``1-point'' flux
! scratch element for extra variables (hardcoded) and conserved variables U
! transposed to quantity-innermost. unvectorizable?
      common /scrns/ zaux (nparm,lx1,ly1,lz1),ut(toteq,lx1,ly1,lz1),
     >               zauxt(lx1*ly1*lz1,nparm),jat(3,3,lx1,ly1,lz1)
     >              ,rhsscr(lx1,toteq)
      real zaux,ut,zauxt,jat,rhsscr
      real flx(5)

      call cmt_usrz(zaux,zauxt,ut,e,nparm)

      call rzero(jat,9*lx1*ly1*lz1)
      do j=1,ldim
         do i=1,ldim
            do iz=1,lz1
            do iy=1,ly1
            do ix=1,lx1
               jat(i,j,ix,iy,iz)=ja(ix,iy,iz,i,j)
            enddo
            enddo
            enddo
         enddo
      enddo

      do iz=1,lz1
      do iy=1,ly1
      do ix=1,lx1
! r-direction
         do l=ix+1,lx1
            call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,l,iy,iz),
     >                            zaux(1,ix,iy,iz),zaux(1,l,iy,iz),
     >                            jat(1,1,ix,iy,iz),jat(1,1,l,iy,iz))
            do eq=1,toteq
            res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+dstrong(ix,l)*flx(eq)
            res(l,iy,iz,e,eq)=res(l,iy,iz,e,eq)+dstrong(l,ix)*flx(eq)
            enddo
         enddo

         do eq=1,toteq
            res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+
     >                       dstrong(ix,ix)*fcons(ix,iy,iz,1,eq)
         enddo

! s-direction
         do l=iy+1,ly1
            call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,ix,l,iz),
     >                            zaux(1,ix,iy,iz),zaux(1,ix,l,iz),
     >                            jat(1,2,ix,iy,iz),jat(1,2,ix,l,iz))
            do eq=1,toteq
            res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+dstrong(iy,l)*flx(eq)
            res(ix,l,iz,e,eq)=res(ix,l,iz,e,eq)+dstrong(l,iy)*flx(eq)
            enddo
         enddo
         do eq=1,toteq
            res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+
     >                       dstrong(iy,iy)*fcons(ix,iy,iz,2,eq)
         enddo

! t-direction
         do l=iz+1,lz1
            call cmt_usr2pt(flx,ut(1,ix,iy,iz),ut(1,ix,iy,l),
     >                            zaux(1,ix,iy,iz),zaux(1,ix,iy,l),
     >                            jat(1,3,ix,iy,iz),jat(1,3,ix,iy,l))
            do eq=1,toteq
            res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+dstrong(iz,l)*flx(eq)
            res(ix,iy,l,e,eq)=res(ix,iy,l,e,eq)+dstrong(l,iz)*flx(eq)
            enddo
         enddo
         do eq=1,toteq
            res(ix,iy,iz,e,eq)=res(ix,iy,iz,e,eq)+
     >                       dstrong(iz,iz)*fcons(ix,iy,iz,3,eq)
         enddo

      enddo ! ix
      enddo ! iy
      enddo ! iz

      return
      end

      subroutine fluxdiv_strong_contra(e)
! JH052818. Evaluate flux divergence of totalh (in contravariant basis)
!           in strong/ultraweak form
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'CMTDATA'

      integer  e,eq
      parameter (ldd=lx1*ly1*lz1)
      parameter (ldg=lx1**3,lwkd=2*ldg)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      common /dgradl/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
      character*32 cname

      nrstd=ldd
      nxyz=lx1*ly1*lz1
      call get_dgll_ptr(ip,lx1,lx1) ! fills dg, dgt
      mxm1=lx1-1

      do eq=1,toteq
      call rzero(ud,nrstd)

      if (if3d) then ! swapping d and dt should do the trick
                     ! in local_grad#_t
         call local_grad3_t(ud,totalh(1,1,eq),totalh(1,2,eq),
     >                         totalh(1,3,eq),mxm1,
     >                      1,dt(ip),d(ip),wkd)
      else
         call local_grad2_t(ud,totalh(1,1,eq),totalh(1,2,eq),mxm1,
     >                      1,dt(ip),d(ip),wkd)
      endif

      call col2   (ud,bm1(1,1,1,e),nxyz)   ! contravariant rx does not
      call invcol2(ud,jacm1(1,1,1,e),nxyz) ! have quadrature weights
      call add2(res1(1,1,1,e,eq),ud,nxyz)
      enddo

      return
      end

      subroutine evaluate_aliased_conv_h(e)
! JH082418 Unstable. not long for this world
! computed as products between primitive variables and conserved variables.
! if you want to write rho u_i u_j as (rho u_i) (rho u_j) (rho^{-1}), this
! is the place to do it
      include  'SIZE'
      include  'SOLN'
      include  'DEALIAS'
      include  'CMTDATA'
      include  'INPUT'

      parameter (ldd=lx1*ly1*lz1)
      common /ctmp1/ ju(ldd),jv(ldd)!,ur(ldd),us(ldd),ud(ldd),tu(ldd)
      real ju,jv
      integer  e,eq

      n=lx1*ly1*lz1

      call copy(ju,phig(1,1,1,e),n)
      call copy(jv,pr(1,1,1,e),n)
      do j=1,ldim
         call copy(convh(1,j,1),u(1,1,1,j+1,e),n)
      enddo

      do eq=2,ldim+1
         call copy(convh(1,1,eq),u(1,1,1,eq,e),n)
         do j=2,ldim
            call copy(convh(1,j,eq),convh(1,1,eq),n)
         enddo
         call col2(convh(1,1,eq),vx(1,1,1,e),n)
         call col2(convh(1,2,eq),vy(1,1,1,e),n)
         if (if3d) call col2(convh(1,3,eq),vz(1,1,1,e),n)
         call add2col2(convh(1,eq-1,eq),ju,jv,n)
      enddo

      eq=toteq
      call copy(convh(1,1,eq),u(1,1,1,eq,e),n)
      call add2col2(convh(1,1,eq),ju,jv,n)
      do j=2,ldim
         call copy(convh(1,j,eq),convh(1,1,eq),n)
      enddo
      call col2(convh(1,1,eq),vx(1,1,1,e),n)
      call col2(convh(1,2,eq),vy(1,1,1,e),n)
      if (if3d) call col2(convh(1,3,eq),vz(1,1,1,e),n)

      return
      end
C> @}

C> \ingroup convhvol
C> @{
C> \f$(\nabla v)\cdot \mathbf{H}^c=\mathcal{I}^{\intercal}\mathbf{D}^{\intercal}\cdots\f$ for equation eq, element e
      subroutine fluxdiv_dealiased_weak_chain(e)
! JH082418. Formerly flux_div_integral_dealiased. Keep this for simpler
!           conservation laws on particularly simple meshes.
!           totalh: correctly formed flux on the Gauss-Legendre (GL) points
!           rx:     from set_delias_rx.
!                   chain-rule nx1 metrics intp'ed to GL point x GL weights
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'CMTDATA'

      integer  e,eq
      parameter (ldd=lxd*lyd*lzd)
      parameter (ldg=lxd**3,lwkd=2*ldg)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      nrstd=ldd
      nxyz=lx1*ly1*lz1
      call get_dgl_ptr(ip,lxd,lxd) ! fills dg, dgt
      mdm1=lxd-1

      do eq=1,toteq
      call rzero(ur,nrstd)
      call rzero(us,nrstd)
      call rzero(ut,nrstd)
      call rzero(ud,nrstd)
      call rzero(tu,nrstd)

      j0=0
      do j=1,ldim
         j0=j0+1
         call add2col2(ur,totalh(1,j,eq),rx(1,j0,e),nrstd)
      enddo
      do j=1,ldim
         j0=j0+1
         call add2col2(us,totalh(1,j,eq),rx(1,j0,e),nrstd)
      enddo
      if (if3d) then
         do j=1,ldim
            j0=j0+1
            call add2col2(ut,totalh(1,j,eq),rx(1,j0,e),nrstd)
         enddo
         call local_grad3_t(ud,ur,us,ut,mdm1,1,dg(ip),dgt(ip),wkd)
      else
         call local_grad2_t(ud,ur,us,   mdm1,1,dg(ip),dgt(ip),wkd)
      endif

      call intp_rstd(tu,ud,lx1,lxd,if3d,1)

! multiply the residue by mass matrix. Non diagonal B should still be
! one block per element
!     call col2(ud,bm1(1,1,1,e),nxyz)  ! res = B*res  --- B=mass matrix
!     call add2(res1(1,1,1,e,eq),tu,nxyz)
! weak?
      call sub2(res1(1,1,1,e,eq),tu,nxyz)
      enddo ! toteq

      return
      end

C> @}

!-----------------------------------------------------------------------

      subroutine fluxdiv_weak_chain(e)
! JH082418. Formerly flux_div_integral_aliased. Keep this for...linear
!           problems on Cartesian meshes? Not long for this world; scalar
!           DG in nek5000 is probably 100 times better anyway
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'CMTDATA'

      integer  e,eq
      parameter (ldd=lx1*ly1*lz1)
      parameter (ldg=lx1**3,lwkd=2*ldg)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      common /dgradl/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
      character*32 cname

      nrstd=ldd
      nxyz=lx1*ly1*lz1
      call get_dgll_ptr(ip,lx1,lx1) ! fills dg, dgt
      mdm1=lx1-1

      do eq=1,toteq
      call rzero(ur,nrstd)
      call rzero(us,nrstd)
      call rzero(ut,nrstd)
      call rzero(ud,nrstd)
      call rzero(tu,nrstd)

      j0=0
      do j=1,ldim
         j0=j0+1
         call add2col2(ur,totalh(1,j,eq),rx(1,j0,e),nrstd)
      enddo
      do j=1,ldim
         j0=j0+1
         call add2col2(us,totalh(1,j,eq),rx(1,j0,e),nrstd)
      enddo
      if (if3d) then
         do j=1,ldim
            j0=j0+1
            call add2col2(ut,totalh(1,j,eq),rx(1,j0,e),nrstd)
         enddo
         call local_grad3_t(ud,ur,us,ut,mdm1,1,d(ip),dt(ip),wkd)
      else
         call local_grad2_t(ud,ur,us,   mdm1,1,d(ip),dt(ip),wkd)
      endif

      call copy(tu,ud,nxyz)

! needs fleg or removal altogether. not good modularity
      call sub2(res1(1,1,1,e,eq),tu,nxyz)
      enddo ! eq

      return
      end

!-----------------------------------------------------------------------

      subroutine contravariant_flux(frst,fxyz,ja,nel)
      include 'SIZE'
      integer e
      real ja(nx1*ny1*nz1,ldim*ldim,nel)
      real frst(nx1*ny1*nz1,ldim,nel)
      real fxyz(nx1*ny1*nz1,ldim,nel)
      parameter (l11=lx1*ly1*lz1)
      common /ctmp1/ ftmp(l11,ldim)
! JH082418 transform to reference element via contravariant metrics

      nxyz=lx1*ly1*lz1
      nxyz3=ldim*nxyz
      do e=1,nel
         call rzero(ftmp,nxyz3)
         j0=0
         do j=1,ldim    ! rst
            do i=1,ldim ! xyz
               j0=j0+1
               call addcol3(ftmp(1,j),ja(1,j0,e),fxyz(1,i,e),nxyz)
            enddo
         enddo
         call copy(frst(1,1,e),ftmp,nxyz3)
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine compute_forcing(e,eq_num)
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'PARALLEL'
      include  'MASS'
      include  'SOLN'
      include  'CMTDATA'
      include  'DEALIAS'
      
      integer e,eq_num
      parameter (ldd=lxd*lyd*lzd)

      common /lpm_fix/ phigdum,phigvdum
      real phigdum(lx1,ly1,lz1,lelt,3),phigvdum(lx1,ly1,lz1,lelt)

      nxyz=lx1*ly1*lz1
      if(eq_num.ne.1.and.eq_num.ne.5)then


        if (eq_num.eq.4.and.ldim.eq.2)then

        else
#ifdef LPM
           call subcol3(res1(1,1,1,e,eq_num),phigdum(1,1,1,e,eq_num-1)
     >                  ,bm1(1,1,1,e),nxyz)
#endif
           call subcol3(res1(1,1,1,e,eq_num),usrf(1,1,1,eq_num)
     $                  ,bm1(1,1,1,e),nxyz) 
        endif
      elseif(eq_num.eq.5)then

#ifdef LPM
           call subcol3(res1(1,1,1,e,eq_num),phigvdum(1,1,1,e)
     >                  ,bm1(1,1,1,e),nxyz)
#endif
c          call subcol3(res1(1,1,1,e,eq_num),usrf(1,1,1,eq_num)
c    $                  ,bm1(1,1,1,e),nxyz) 

      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine cmtusrf(e)
      include 'SIZE'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'TSTEP'
      include 'PARALLEL'

      integer e,eg

      if(istep.eq.1)then
        n = lx1*ly1*lz1*5
        call rzero(usrf,n)
      endif
      eg = lglel(e)
      do k=1,lz1
         do j=1,ly1
            do i=1,lx1
               call NEKASGN(i,j,k,e)
               call userf(i,j,k,eg)
               rdum4 = 0.
#ifdef LPM
               call lpm_userf(I,J,K,e,rdum1,rdum2,rdum3,rdum4)
               FFX  = FFX + rdum1
               FFY  = FFY + rdum2
               FFZ  = FFZ + rdum3
#endif
               ! note fx,fy,fz multiply by density to stay 
               ! consistent with nek5000 units. Same for phig (cancels)
               usrf(i,j,k,2) = FFX*u(i,j,k,1,e)*phig(i,j,k,e)
               usrf(i,j,k,3) = FFY*u(i,j,k,1,e)*phig(i,j,k,e)
               usrf(i,j,k,4) = FFZ*u(i,j,k,1,e)*phig(i,j,k,e)
               usrf(i,j,k,5) = 0.0
c              usrf(i,j,k,5) = (U(i,j,k,2,e)*FFX + U(i,j,k,3,e)*FFY
c    &                       +  U(i,j,k,4,e)*FFZ)/ U(i,j,k,1,e)
            enddo
         enddo
      enddo

      return
      end 

