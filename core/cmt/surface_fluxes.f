C> @file surface_fluxes.f Routines for surface terms on RHS.
C> \ingroup isurf
C> @{
C> overwrite beginning of /CMTSURFLX/ with -[[U]] for viscous terms
      subroutine fillujumpu
!-----------------------------------------------------------------------
! JH091319 Yes, I know things like llf_euler already did this, but
!          1.) I decided not to dedicate memory to [[U]] in addition
!              to primitive variables and
!          2.) It isn't very modular to store [[U]] from split forms
!              that need LLF. We would still need to do this again for
!              entropy-stable fluxes!
! JH092319 Dirichlet BC for viscous fluxes go here
!-----------------------------------------------------------------------
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
      include 'INPUT'

! JH070219 "heresize" and "hdsize" come from a failed attempt at managing
!          memory in CMTSURFLX by redeclaration that was abandoned before
!          the two-point split form. They need to be taken care of in
!          CMTSIZE and consistent with the desired subroutine
      common /CMTSURFLX/ fatface(heresize),uplus(hdsize)
      real fatface,uplus
      integer eq
      character*32 cname

      nfq=lx1*lz1*2*ldim*nelt
      nstate = toteq
! where different things live
      iwm =1
      iwp =iwm+nstate*nfq  ! duplicate one conserved variable at a time for jumps in LLF
                               ! tuck it before jph=nqq
!     i_cvars=(ju1-1)*nfq+1
      i_cvars=1
      do eq=1,toteq
         call faceu(eq,fatface(i_cvars))
! JH091319 still need to do something about phi; does jph=nqq still?
!        call invcol2(fatface(i_cvars),fatface(iwm+nfq*(jph-1)),nfq)
         i_cvars=i_cvars+nfq
      enddo
      call face_state_commo(fatface(iwm),uplus,nfq,nstate,dg_hndl)
      call sub3(fatface(iwm),uplus,fatface(iwm),nstate*nfq)

C> @}
      return
      end

C> \ingroup isurf
C> @{
C> Restrict and copy face data and compute inviscid numerical flux 
C> \f$\oint \mathbf{H}^{c\ast}\cdot\mathbf{n}dA\f$ on face points
!     subroutine fluxes_full_field(fstab,parameter_vector,fsharp)
! still not sure how to abstract these three functions out
!     subroutine fluxes_full_field_kepec
!     iwp =iwm+(nstate-1)*nfq  ! you're going to try to be clever later
!     iflx=iwm+nstate*nfq

      subroutine fluxes_full_field_kg
!-----------------------------------------------------------------------
! JH041419 Kennedy and Gruber (2008) split form adapted to DGSEM and
!          fully vectorized for state -> parameter vector -> gs_op ->
!          flux. Local Lax-Friedrichs for stabilization/penalizing jumps.
!-----------------------------------------------------------------------
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
      include 'INPUT'

! JH070219 "heresize" and "hdsize" come from a failed attempt at managing
!          memory in CMTSURFLX by redeclaration that was abandoned before
!          the two-point split form. They need to be taken care of in
!          CMTSIZE and consistent with the desired subroutine
      common /CMTSURFLX/ fatface(heresize),notyet(hdsize)
      real fatface,notyet
      integer eq
      character*32 cname

      nfq=lx1*lz1*2*ldim*nelt
      nstate = nqq
! where different things live
      iwm =1
      iwp =iwm+nstate*nfq  ! duplicate one conserved variable at a time for jumps in LLF
                               ! tuck it before jph=nqq
      iflx=iwp+nstate*nfq

      call rzero(fatface(iflx),nfq*toteq)

! just multipy by {{phi}}
! until we can verify correct multiphase two-point fluxes
!     call faceu(1,fatface(iwm+nfq*(jrhof-1)))

! stabilization first
! LLF for equations 1 through 4 (mass & momentum)
      call fillq(jux, vx,    fatface(iwm),fatface(iwp))
      call fillq(juy, vy,    fatface(iwm),fatface(iwp))
      call fillq(juz, vz,    fatface(iwm),fatface(iwp))
      call fillq(jph, phig,  fatface(iwm),fatface(iwp))
      call fillq(jsnd,csound,fatface(iwm),fatface(iwp))

      call llf_euler_vec(fatface(iwm),fatface(iwp),fatface(iflx),nstate)

      call fillq(jpr, pr,    fatface(iwm),fatface(iwp))
! ONLY needed by Kennedy-Gruber (2008) as written. This is done in fstab
! for Chandrashekar (2013), but overwrites jsnd for KG and friends.
      call fillq(jrhof,vtrans,fatface(iwm),fatface(iwp))
! q- -> z-. Kennedy-Gruber, Pirozzoli, and most energy-
!           conserving fluxes have z=q, so I just divide total energy by
!           U1 here since Kennedy-Gruber needs E

      call rhoe_to_e(fatface(iwm),nfq,nstate)

! z- -> z^, which is {{z}} for Kennedy-Gruber, Pirozzoli, and some parts
!           of other energy-conserving fluxes.
      call dg_face_avg(fatface(iwm),nfq,nstate,dg_hndl)

! z^ -> F#. Some parameter-vector stuff can go here too as long as it's all
!           local to a given element.
      call kennedygruber_vec(fatface(iwm),fatface(iflx),nstate,toteq)

!     i_cvars=iwm!(ju1-1)*nfq+1
!     do eq=1,toteq
!        call faceu(eq,fatface(i_cvars))
!        i_cvars=i_cvars+nfq
!     enddo

C> @}
! Now do all fluxes for all boundaries, both F# and stabilized
      call InviscidBC(fatface(iflx))

      return
      end

!-----------------------------------------------------------------------


      subroutine faceu(ivar,yourface)
! get faces of conserved variables stored contiguously
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'
      integer e
      real yourface(lx1,lz1,2*ldim,nelt)

      do e=1,nelt
         call full2face_cmt(1,lx1,ly1,lz1,iface_flux(1,e),
     >                      yourface(1,1,1,e),u(1,1,1,ivar,e))
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine fillq(jvar,field,wminus,yourface)
      include 'SIZE'
      include 'DG'

      integer jvar! intent(in)
      real field(lx1,ly1,lz1,nelt)! intent(in)
!     real, intent(out)wminus(7,lx1*lz1*2*ldim*nelt) ! gs_op no worky
      real wminus(lx1*lz1*2*ldim*nelt,*)! intent(out)
      real yourface(lx1*lz1*2*ldim*nelt)
      integer e,f

      call full2face_cmt(nelt,lx1,ly1,lz1,iface_flux,yourface,field)

      do i=1,ndg_face ! need to make sure ndg_face=lx1*lz1*2*ldim*nelt
         wminus(i,jvar)=yourface(i)
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine dg_face_avg(mine,nf,nstate,handle)

! JH110818 A lot of entropy-stable fluxes have a product of averages.
!          mine: starting parameter vector of nstate quantities,
!                quantity outermost, on piles of faces with nf points

      integer handle,nf,nstate ! intent(in)
      real mine(*)

      ntot=nf*nstate
      call cmult(mine,0.5,ntot)
!-----------------------------------------------------------------------
! operation flag is second-to-last arg, an integer
!                                                1 ==> +
      call fgslib_gs_op_fields(handle,mine,nf,nstate,1,1,0)
      return
      end

!-----------------------------------------------------------------------

      subroutine face_state_commo(mine,yours,nf,nstate,handle)

! JH060414 if we ever want to be more intelligent about who gets what,
!          who gives what and who does what, this is the place where all
!          that is done. At the very least, gs_op may need the transpose
!          flag set to 1. Who knows. Everybody duplicates everything for
!          now.
! JH070714 figure out gs_op_fields, many, vec, whatever (and the
!          corresponding setup) to get this done for the transposed
!          ordering of state variables. I want variable innermost, not
!          grid point.

      integer handle,nf,nstate ! intent(in)
      real yours(*),mine(*)

      ntot=nf*nstate
      call copy(yours,mine,ntot)
!-----------------------------------------------------------------------
! operation flag is second-to-last arg, an integer
!                                                1 ==> +
      call fgslib_gs_op_fields(handle,yours,nf,nstate,1,1,0)
      call sub2 (yours,mine,ntot)
      return
      end

!-----------------------------------------------------------------------

      subroutine avg_and_jump(avg,jump,scratch,nf,nstate,handle)

! JH011419 Get the most out of every gs_op. pile of faces of U in avg
!          gets overwritten by {{U}}. jump gets [[U]].

      integer handle,nf,nstate ! intent(in)
      real avg(*) !, intent(inout) ! U- on input, {{U}} on output
      real jump(*),scratch(*) !, intent(out) ! jump becomes [[U]]

      ntot=nf*nstate
      call copy(scratch,avg,ntot)
!-----------------------------------------------------------------------
! operation flag is second-to-last arg, an integer
!                                                1 ==> +
      call fgslib_gs_op_fields(handle,scratch,nf,nstate,1,1,0)
      call copy(jump,scratch,ntot) ! jump = U- + U+
      const=-2.0
      call add2s2(jump,avg,const,ntot) ! jump = (U- + U+) - 2*U- = U+-U-
      call cmult(scratch,0.5,ntot)
      call copy(avg,scratch,ntot) ! avg={{scratch}}
      return
      end

!-----------------------------------------------------------------------

      subroutine face_flux_commo(flux1,flux2,nf,neq,handle)
! JH060514 asymmetric transposed gs_op, gs_unique magic may be needed if
!          we ever decide to avoid redundancy. For now, this routine
!          doesn't need to do anything.
      integer ntot,handle
      real flux1(*),flux2(*)
! JH061814 It doesn't need to do anything, but a sanity check would be
!          wise.
      return
      end

!-------------------------------------------------------------------------------

      subroutine sequential_flux(flux,wminus,wplus,uminus,uplus,
     >                           jaminus,japlus,
     >                           fluxfunction,nstate,npt)
! Calls two-point flux functions one point at a time for npt points for which both
! points are given. Mostly intended to allow quantity-innermost volume flux
! functions to be used where needed for boundary points too, after *bc routines
! provide Dirichlet ``rind'' states in wplus and uplus.
      include 'SIZE'
      include 'CMTSIZE'
      real flux(toteq,npt),wminus(nstate,npt),wplus(nstate,npt),
     >   jaminus(3,npt),japlus(3,npt),uminus(toteq,npt),uplus(toteq,npt)
      external fluxfunction

      do i=1,npt
         call fluxfunction(flux(1,i),uminus(1,i),uplus(1,i),wminus(1,i),
     >                     wplus(1,i),jaminus(1,i),japlus(1,i))
      enddo

      return
      end

!-------------------------------------------------------------------------------
! OBSOLETE ROUTINES
!-------------------------------------------------------------------------------

      subroutine InviscidFlux(wminus,wplus,flux,nstate,nflux)
!-------------------------------------------------------------------------------
! JH091514 A fading copy of RFLU_ModAUSM.F90 from RocFlu
!-------------------------------------------------------------------------------

!#ifdef SPEC
!      USE ModSpecies, ONLY: t_spec_type
!#endif
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'GEOM' ! for unx
      include 'CMTDATA' ! do we need this without outflsub? It has toteq now
      include 'TSTEP' ! for ifield?
      include 'DG'

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real wminus(lx1*lz1,2*ldim,nelt,nstate),
     >     wplus(lx1*lz1,2*ldim,nelt,nstate),
     >     flux(lx1*lz1,2*ldim,nelt,nflux)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface
      parameter (lfd=lxd*lzd)
! JH111815 legacy rocflu names.
!
! nx,ny,nz : outward facing unit normal components
! fs       : face speed. zero until we have moving grid
! jaco_c   : fdim-D GLL grid Jacobian
! nm       : jaco_c, fine grid
!
! State on the interior (-, "left") side of the face
! rl       : density
! ul,vl,wl : velocity
! tl       : temperature
! al       : sound speed
! pl       : pressure, then phi
! el      : rho*cp
! State on the exterior (+, "right") side of the face
! rr       : density
! ur,vr,wr : velocity
! tr       : temperature
! ar       : sound speed
! pr       : pressure
! er      : rho*cp

      COMMON /SCRNS/ nx(lfd), ny(lfd), nz(lfd), rl(lfd), ul(lfd),
     >               vl(lfd), wl(lfd), pl(lfd), tl(lfd), al(lfd),
     >               el(lfd),rr(lfd), ur(lfd), vr(lfd), wr(lfd),
     >               pr(lfd),tr(lfd), ar(lfd),er(lfd),phl(lfd),fs(lfd),
     >               jaco_f(lfd),flx(lfd,toteq),jaco_c(lx1*lz1)
      real nx, ny, nz, rl, ul, vl, wl, pl, tl, al, el, rr, ur, vr, wr,
     >                pr,tr, ar,er,phl,fs,jaco_f,flx,jaco_c

!     REAL vf(3)
      real nTol
      character*132 deathmessage
      common /nekcb/ cb
      character*3 cb

      nTol = 1.0E-14

      fdim=ldim-1
      nface = 2*ldim
      nxz   = lx1*lz1
      nxzd  = lxd*lzd
      ifield= 1 ! You need to figure out the best way of dealing with
                ! this variable

!     if (outflsub)then
!        call maxMachnumber
!     endif
      do e=1,nelt
      do f=1,nface
! diagnostic
           call facind (kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,f)

! JH021717 Finally corrected BC wrongheadedness. Riemann solver handles
!          all fluxes with even the slightest Dirichlet interpretation,
!          and BC fill wplus sanely for the Riemann solver to provide
!          a good flux for weak BC.
! JH111715 now with dealiased surface integrals. I am too lazy to write
!          something better

         if (lxd.gt.lx1) then
            call map_faced(nx,unx(1,1,f,e),lx1,lxd,fdim,0)
            call map_faced(ny,uny(1,1,f,e),lx1,lxd,fdim,0)
            call map_faced(nz,unz(1,1,f,e),lx1,lxd,fdim,0)

            call map_faced(rl,wminus(1,f,e,jrhof),lx1,lxd,fdim,0)
            call map_faced(ul,wminus(1,f,e,jux),lx1,lxd,fdim,0)
            call map_faced(vl,wminus(1,f,e,juy),lx1,lxd,fdim,0)
            call map_faced(wl,wminus(1,f,e,juz),lx1,lxd,fdim,0)
            call map_faced(pl,wminus(1,f,e,jpr),lx1,lxd,fdim,0)
            call map_faced(tl,wminus(1,f,e,jthm),lx1,lxd,fdim,0)
            call map_faced(al,wminus(1,f,e,jsnd),lx1,lxd,fdim,0)
            call map_faced(el,wminus(1,f,e,jen),lx1,lxd,fdim,0)

            call map_faced(rr,wplus(1,f,e,jrhof),lx1,lxd,fdim,0)
            call map_faced(ur,wplus(1,f,e,jux),lx1,lxd,fdim,0)
            call map_faced(vr,wplus(1,f,e,juy),lx1,lxd,fdim,0)
            call map_faced(wr,wplus(1,f,e,juz),lx1,lxd,fdim,0)
            call map_faced(pr,wplus(1,f,e,jpr),lx1,lxd,fdim,0)
            call map_faced(tr,wplus(1,f,e,jthm),lx1,lxd,fdim,0)
            call map_faced(ar,wplus(1,f,e,jsnd),lx1,lxd,fdim,0)
            call map_faced(er,wplus(1,f,e,jen),lx1,lxd,fdim,0)

            call map_faced(phl,wminus(1,f,e,jph),lx1,lxd,fdim,0)

            call invcol3(jaco_c,area(1,1,f,e),w2m1,nxz)
            call map_faced(jaco_f,jaco_c,lx1,lxd,fdim,0) 
            call col2(jaco_f,wghtf,nxzd)
         else

            call copy(nx,unx(1,1,f,e),nxz)
            call copy(ny,uny(1,1,f,e),nxz)
            call copy(nz,unz(1,1,f,e),nxz)

            call copy(rl,wminus(1,f,e,jrhof),nxz)
            call copy(ul,wminus(1,f,e,jux),nxz)
            call copy(vl,wminus(1,f,e,juy),nxz)
            call copy(wl,wminus(1,f,e,juz),nxz)
            call copy(pl,wminus(1,f,e,jpr),nxz)
            call copy(tl,wminus(1,f,e,jthm),nxz)
            call copy(al,wminus(1,f,e,jsnd),nxz)
            call copy(el,wminus(1,f,e,jcpf),nxz)

            call copy(rr,wplus(1,f,e,jrhof),nxz)
            call copy(ur,wplus(1,f,e,jux),nxz)
            call copy(vr,wplus(1,f,e,juy),nxz)
            call copy(wr,wplus(1,f,e,juz),nxz)
            call copy(pr,wplus(1,f,e,jpr),nxz)
            call copy(tr,wplus(1,f,e,jthm),nxz)
            call copy(ar,wplus(1,f,e,jsnd),nxz)
            call copy(er,wplus(1,f,e,jen),nxz)

            call copy(phl,wminus(1,f,e,jph),nxz)

            call copy(jaco_f,jface(1,1,f,e),nxz) 
         endif
         call rzero(fs,nxzd) ! moving grid stuff later

         call AUSM_FluxFunction(nxzd,nx,ny,nz,jaco_f,fs,rl,ul,vl,wl,pl,
     >                          al,tl,rr,ur,vr,wr,pr,ar,tr,flx,el,er)

         do j=1,toteq
            call col2(flx(1,j),phl,nxzd)
         enddo

         if (lxd.gt.lx1) then
            do j=1,toteq
               call map_faced(flux(1,f,e,j),flx(1,j),lx1,lxd,fdim,1)
            enddo
         else
            do j=1,toteq
               call copy(flux(1,f,e,j),flx(1,j),nxz)
            enddo
         endif

      enddo
      enddo

      end

!-----------------------------------------------------------------------

      subroutine surface_integral_full(vol,flux)
! Integrate surface fluxes for an entire field. Add contribution of flux
! to volume according to add_face2full_cmt
      include 'SIZE'
      include 'GEOM'
      include 'DG'
      include 'CMTDATA'
      real vol(lx1*ly1*lz1*nelt),flux(*)
      integer e,f

! weak form until we get the time loop rewritten
!     onem=-1.0
!     ntot=lx1*lz1*2*ldim*nelt
!     call cmult(flux,onem,ntot)
! weak form until we get the time loop rewritten
      call add_face2full_cmt(nelt,lx1,ly1,lz1,iface_flux,vol,flux)

      return
      end

!-------------------------------------------------------------------------------

      subroutine diffh2graduf(e,eq,graduf)
! peels off diffusiveH into contiguous face storage via restriction operator R
! for now, stores gradU.n for QQT in igu. ONLY 5 FLUXES, for {{AgradU.n}}
! 
      include  'SIZE'
      include  'DG' ! iface
      include  'CMTDATA'
      include  'GEOM'
      integer e,eq
      real graduf(lx1*lz1*2*ldim,nelt,toteq)
      common /scrns/ hface(lx1*lz1,2*ldim)
     >              ,normal(lx1*ly1,2*ldim)
      real hface, normal

      integer f

      nf    = lx1*lz1*2*ldim*nelt
      nfaces=2*ldim
      nxz   =lx1*lz1
      nxzf  =nxz*nfaces
      nxyz  = lx1*ly1*lz1

      call rzero(graduf(1,e,eq),nxzf) !   . dot nhat -> overwrites beginning of flxscr
      do j =1,ldim
         if (j .eq. 1) call copy(normal,unx(1,1,1,e),nxzf)
         if (j .eq. 2) call copy(normal,uny(1,1,1,e),nxzf)
         if (j .eq. 3) call copy(normal,unz(1,1,1,e),nxzf)
         call full2face_cmt(1,lx1,ly1,lz1,iface_flux,hface,diffh(1,j)) 
         call addcol3(graduf(1,e,eq),hface,normal,nxzf)
      enddo
      call col2(graduf(1,e,eq),area(1,1,1,e),nxzf)

      return
      end

!-------------------------------------------------------------------------------

      subroutine diffh2face(e,eq,diffhf)
! peels off diffusiveH into contiguous face storage via restriction operator R
! for final central flux in BR1. ALL 15 FLUXES, for {{AgradU}}.n
      include  'SIZE'
      include  'DG' ! iface
      include  'CMTDATA'
      include  'GEOM'
      integer e,eq
! must not exceed hdsize
! faster dot product {{Hd}}.n
      real diffhf(lx1*lz1*2*ldim,nelt,ldim,toteq)
! vs allowing diffhf(:,:,:,1) to overwrite [[U]]
!     real diffhf(lx1*lz1*2*ldim,nelt,toteq,ldim)
      common /scrns/ hface(lx1*lz1,2*ldim)
     >              ,normal(lx1*ly1,2*ldim)
      real hface, normal

      integer f

      nf    = lx1*lz1*2*ldim*nelt
      nfaces=2*ldim
      nxz   =lx1*lz1
      nxzf  =nxz*nfaces
      nxyz  = lx1*ly1*lz1

      do j =1,ldim
         call full2face_cmt(1,lx1,ly1,lz1,iface_flux,hface,diffh(1,j)) 
         call copy(diffhf(1,e,j,eq),hface,nxzf)
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine igu_cmt(flxscr,gdudxk,wminus)
! Hij^{d*}
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'

      real flxscr(lx1*lz1*2*ldim*nelt,toteq)
      real gdudxk(lx1*lz1*2*ldim,nelt,toteq)
      real wminus(lx1*lz1,2*ldim,nelt,nqq)
      real const
      integer e,eq,f

      nxz = lx1*lz1
      nfaces=2*ldim
      nxzf=nxz*nfaces
      nfq =lx1*lz1*nfaces*nelt
      ntot=nfq*toteq

      call copy (flxscr,gdudxk,ntot) ! save AgradU.n
      const = 0.5
      call cmult(gdudxk,const,ntot)
!-----------------------------------------------------------------------
! supa huge gs_op to get {{AgradU}}
! operation flag is second-to-last arg, an integer
!                                                   1 ==> +
      call fgslib_gs_op_fields(dg_hndl,gdudxk,nfq,toteq,1,1,0)
!-----------------------------------------------------------------------
      call sub2  (flxscr,gdudxk,ntot) ! overwrite flxscr with
                                      !           -
                                      ! (AgradU.n)  - {{AgradU.n}}
! [v] changes character on boundaries, so we actually need
! 1. (AgradU.n)- on Dirichlet boundaries
      call igu_dirichlet(flxscr,gdudxk)
! 2. (Fbc.n)- on Neumann boundaries
      call bcflux(flxscr,gdudxk,wminus)
      call chsign(flxscr,ntot) ! needs to change with sign changes

      return
      end

!-----------------------------------------------------------------------

      subroutine igu_dirichlet(flux,agradu)
! Acts on ALL boundary faces because I'm lazy. SET NEUMANN BC AFTER THIS
! CALL. BCFLUX IS PICKIER ABOUT THE BOUNDARY FACES IT ACTS ON.
      include 'SIZE'
      include 'CMTSIZE'
      include 'TOTAL'
      integer e,eq,f
      real flux(lx1*lz1,2*ldim,nelt,toteq)
      real agradu(lx1*lz1,2*ldim,nelt,toteq)
      character*3 cb2

      nxz=lx1*lz1
      nfaces=2*ldim

      ifield=1
      do e=1,nelt
         do f=1,nfaces
            cb2=cbc(f, e, ifield)
            if (cb2.ne.'E  '.and.cb2.ne.'P  ') then ! cbc bndy.
! all Dirichlet conditions result in IGU being
! strictly one-sided, so we undo 0.5*QQT
! UNDER THE ASSUMPTIONS THAT
! 1. agradu's actual argument is really gdudxk AND
! 2. IT HAS ALREADY BEEN MULTIPLIED BY 0.5
! 3. gs_op has not changed it at all.
! overwriting flux with it and and multiplying it 2.0 should do the trick
               do eq=1,toteq
!                  call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
!! in fact, that copy should not be necessary at all. TEST WITHOUT IT
!                  call cmult(flux(1,f,e,eq),2.0,nxz)
! JH112216 This may be better because agradu (without the factor of 1/2) is
!          needed for some Neumann conditions (like adiabatic walls)
                   call cmult(agradu(1,f,e,eq),2.0,nxz)
                   call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
               enddo
            endif
         enddo
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine br1primary(flux,gdudxk)
! Hij^{d*}=JA{{Hij_d}}.n
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'

      real flux(lx1*lz1*2*ldim*nelt,toteq)
!     real gdudxk(lx1*lz1*2*ldim,nelt,toteq)
      real gdudxk(lx1*lz1*2*ldim*nelt,ldim,toteq)
      real const
      integer e,eq,f

      nxz = lx1*lz1
      nfaces=2*ldim
      nxzf=nxz*nfaces
      nfq =lx1*lz1*nfaces*nelt
      ntot=nfq*toteq

      const = 0.5
      call cmult(gdudxk,const,ntot*ldim)
!-----------------------------------------------------------------------
! supa huge gs_op to get {{AgradU}}
! operation flag is second-to-last arg, an integer
!                                                   1 ==> +
      call fgslib_gs_op_fields(dg_hndl,gdudxk,nfq,toteq*ldim,1,1,0)
      do eq=1,toteq ! {{AgradU}}.n
         call agradu_normal_flux(flux(1,eq),gdudxk(1,1,eq))
      enddo
      call br1bc(flux)
      call chsign(flux,ntot) ! needs to change with sign changes

      return
      end

!-----------------------------------------------------------------------

      subroutine agradu_normal_flux(flux,graduf)
! JA{{AgradU}}.n
      include  'SIZE'
      include 'INPUT'
      include  'CMTDATA'
      include  'GEOM'
      integer e
      integer f
      real graduf(lx1*lz1*2*ldim,nelt,ldim)
      real flux(lx1*lz1*2*ldim,nelt)
      common /SCRNS/ jscr(lfq)
      real jscr

      nfaces=2*ldim
      nxz=lx1*lz1
      nxzf  =nxz*nfaces
      nf=nxzf*nelt

! I don't know what to do with volume fraction phi, and this is my first guess
!     call col3(jscr,jface,z(1,jph),nf) ! Jscr=JA*{{\phi_g}}
! GET phi in HERE SOMEDAY SOON
      call copy(jscr,jface,nf)
! zero out jscr at boundary faces; gs_op is degenerate there.
!     call bcmask_cmt(jscr) ! RECONCILE THIS WITH br1bc!!!

      nxyz  = lx1*ly1*lz1

      do e =1,nelt
         call col3(flux(1,e),graduf(1,e,1),unx(1,1,1,e),nxzf)
      enddo
      do e =1,nelt
         call addcol3(flux(1,e),graduf(1,e,2),uny(1,1,1,e),nxzf)
      enddo
      if (if3d) then
         do e =1,nelt
            call addcol3(flux(1,e),graduf(1,e,3),unz(1,1,1,e),nxzf)
         enddo
      endif
      call col2(flux,jscr,nf)

      return
      end

!-----------------------------------------------------------------------

      subroutine br1bc(flux)
!-----------------------------------------------------------------------
! Two boundary condition operations remain even for artificial viscosity
! 1. (AgradU)- on Dirichlet boundaries has been multiplied by 1/2, but
! needs to be doubled because gs_op has added zeros to 1/2(AgradU)
! 2. (Fbc.n)- on Neumann boundaries via bcflux_br1
      include 'SIZE'
      include 'CMTSIZE'
      include 'TOTAL'
      integer e,eq,f
      real flux(lx1*lz1,2*ldim,nelt,toteq)
      character*3 cb2

      nxz=lx1*lz1
      nfaces=2*ldim

      ifield=1
      do e=1,nelt
         do f=1,nfaces
            cb2=cbc(f, e, ifield)
            if (cb2.ne.'E  '.and.cb2.ne.'P  ') then ! cbc bndy.
! all Dirichlet conditions result in {{}} being
! UNDER THE ASSUMPTIONS THAT
! 1. agradu's actual argument is really gdudxk AND
! 2. IT HAS ALREADY BEEN MULTIPLIED BY 0.5
! 3. gs_op has not changed it at all.
! overwriting flux with agradu and and multiplying it 2.0 should do the trick
               do eq=1,toteq
! JH092719 already dotted with normal. so just multiply by 2
                   call cmult(flux(1,f,e,eq),2.0,nxz)
               enddo
               call bcflux_br1(flux,f,e)
            endif
         enddo
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine bcflux_br1(flux,f,e)
! JH092719. Placeholder for Neumann conditions in an actual Navier-Stokes
!           solver. For artificial viscosity, we want zero viscous fluxes
!           at all boundaries.
      include 'SIZE'
      include 'CMTSIZE'
      integer f,e
      integer eq
      real flux(lx1*lz1,2*ldim,nelt,toteq)
      logical ifadiabatic
      data ifadiabatic /.false./ 
      return
      end

!-----------------------------------------------------------------------

      subroutine strong_sfc_flux(flux,vflx,e,eq)
! JH052818 dedicated routine for stripping faces from an array vflx of
!          physical fluxes on GLL faces and storing them in the flux
!          pile of faces for a call to surface_integral_full
      include 'SIZE'
      include 'CMTSIZE'
      include 'INPUT'
      include 'GEOM' ! for unx
      include 'DG'
      
      parameter (lf=lx1*lz1*2*ldim)
      COMMON /SCRNS/ yourface(lf),normal(lf)
      real yourface,normal
      real flux(lx1*lz1,2*ldim,nelt,toteq) ! out
      real vflx(lx1*ly1*lz1,ldim)          ! in
      integer e,eq
      integer f

      nxz =lx1*lz1
      nxzf=nxz*2*ldim
      nfaces=2*ldim

      call rzero(flux(1,1,e,eq),nxzf)

      do j=1,ldim
         if (j .eq. 1) call copy(normal,unx(1,1,1,e),nxzf)
         if (j .eq. 2) call copy(normal,uny(1,1,1,e),nxzf)
         if (j .eq. 3) call copy(normal,unz(1,1,1,e),nxzf)
         call full2face_cmt(1,lx1,ly1,lz1,iface_flux(1,e),yourface,
     >                      vflx(1,j))
         call col2(yourface,normal,nxzf)
!        call add2(flux(1,1,e,eq),yourface,nxzf) ! needed for +sign in RK loop
         call sub2(flux(1,1,e,eq),yourface,nxzf)
      enddo

!     call col2(flux(1,1,e,eq),area(1,1,1,e),nxzf)
      do f=1,nfaces
         call col2(flux(1,f,e,eq),w2m1,nxz)
      enddo

      return
      end

!-----------------------------------------------------------------------
! JUNKYARD
!-----------------------------------------------------------------------

      subroutine fluxes_full_field_chold
!-----------------------------------------------------------------------
! Chandrashekar KEPEC flux done in the most redundant and expensive way
! possible. not recommended
!-----------------------------------------------------------------------
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
      include 'INPUT'

      common /CMTSURFLX/ fatface(heresize),notyet(hdsize)
      real fatface,notyet
      integer eq
      character*32 cname

      nstate = nqq
      nfq=lx1*lz1*2*ldim*nelt
! where different things live
      iwm =1
      iwp =iwm+nstate*nfq
      iflx=iwp+nstate*nfq

      call rzero(fatface(iflx),nfq*toteq)
 
      call fillq(jrhof,vtrans,fatface(iwm),fatface(iwp))
      call fillq(jux, vx,    fatface(iwm),fatface(iwp))
      call fillq(juy, vy,    fatface(iwm),fatface(iwp))
      call fillq(juz, vz,    fatface(iwm),fatface(iwp))
      call fillq(jpr, pr,    fatface(iwm),fatface(iwp))
      call fillq(jsnd,csound,fatface(iwm),fatface(iwp))
      call fillq(jph, phig,  fatface(iwm),fatface(iwp))

      call llf_euler_vec(fatface(iwm),fatface(iwp),fatface(iflx),nstate)

! Duplicate everything
      call face_state_commo(fatface(iwm),fatface(iwp),nfq,nstate
     >                     ,dg_hndl)

! Now do all fluxes for all boundaries, both F# and stabilized
!     call InviscidBC(fatface(iflx))
! just for periodic test cases
      call KEPEC_duplicated(fatface(iwm),fatface(iwp),fatface(iflx))
!     call InviscidBC(fatface(iflx)) ! roll back to remember where you did this for KG

      return
      end

      subroutine fluxes_full_field_old
!-----------------------------------------------------------------------
! JH060314 First, compute face fluxes now that we have the primitive variables
! JH091514 renamed from "surface_fluxes_inviscid" since it handles all fluxes
!          that we compute from variables stored for the whole field (as
!          opposed to one element at a time).
! JH070918 redone for two-point fluxes 
!-----------------------------------------------------------------------
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
      include 'INPUT'

      common /CMTSURFLX/ fatface(heresize),notyet(hdsize)
      real fatface,notyet
      integer eq
      character*32 cname

      nstate = nqq
      nfq=lx1*lz1*2*ldim*nelt
! where different things live
      iwm =1
      iwp =iwm+nstate*nfq
      iflx=iwp+nstate*nfq
 
      call fillq(jrhof,vtrans,fatface(iwm),fatface(iwp))
      call fillq(jux, vx,    fatface(iwm),fatface(iwp))
      call fillq(juy, vy,    fatface(iwm),fatface(iwp))
      call fillq(juz, vz,    fatface(iwm),fatface(iwp))
      call fillq(jpr, pr,    fatface(iwm),fatface(iwp))
      call fillq(jthm,t,     fatface(iwm),fatface(iwp))
      call fillq(jsnd,csound,fatface(iwm),fatface(iwp))
      call fillq(jph, phig,  fatface(iwm),fatface(iwp))
      call fillq(jcvf,vtrans(1,1,1,1,jcv),fatface(iwm),fatface(iwp))
      call fillq(jcpf,vtrans(1,1,1,1,jcp),fatface(iwm),fatface(iwp))
!     call fillq(jen,vtrans(1,1,1,1,jen),fatface(iwm),fatface(iwp))
      call fillq(jenf,vtrans(1,1,1,1,jen),fatface(iwm),fatface(iwp))
      call fillq(jmuf, vdiff(1,1,1,1,jmu), fatface(iwm),fatface(iwp))
      call fillq(jkndf,vdiff(1,1,1,1,jknd),fatface(iwm),fatface(iwp))
      call fillq(jlamf,vdiff(1,1,1,1,jlam),fatface(iwm),fatface(iwp))

      i_cvars=(ju1-1)*nfq+1
      do eq=1,toteq
         call faceu(eq,fatface(i_cvars))
! JH080317 at least get the product rule right until we figure out how
!          we want the governing equations to look
         call invcol2(fatface(i_cvars),fatface(iwm+nfq*(jph-1)),nfq)
         i_cvars=i_cvars+nfq
      enddo
      call face_state_commo(fatface(iwm),fatface(iwp),nfq,nstate
     >                     ,dg_hndl)
      call InviscidBC(fatface(iwm),fatface(iwp),nstate)
!     call InviscidFlux(fatface(iwm),fatface(iwp),fatface(iflx)
      call InviscidFluxRot(fatface(iwm),fatface(iwp),fatface(iflx)
     >                 ,nstate,toteq)
      return
      end

!-----------------------------------------------------------------------

      subroutine InviscidFluxRot(wminus,wplus,flux,nstate,nflux)
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'GEOM' ! for unx
      include 'CMTDATA' ! do we need this without outflsub?
      include 'TSTEP' ! for ifield?
      include 'DG'

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real wminus(lx1*lz1,2*ldim,nelt,nstate),
     >     wplus(lx1*lz1,2*ldim,nelt,nstate),
     >     flux(lx1*lz1,2*ldim,nelt,nflux)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface
      parameter (lfd=lxd*lzd)
! JH111815 legacy rocflu names.
!
! nx,ny,nz : outward facing unit normal components
! fs       : face speed. zero until we have moving grid
! jaco_c   : fdim-D GLL grid Jacobian
! nm       : jaco_c, fine grid
!
! State on the interior (-, "left") side of the face
! rl       : density
! ul,vl,wl : velocity
! tl       : temperature
! al       : sound speed
! pl       : pressure, then phi
! el      : rho*cp
! State on the exterior (+, "right") side of the face
! rr       : density
! ur,vr,wr : velocity
! tr       : temperature
! ar       : sound speed
! pr       : pressure
! er      : rho*cp

      COMMON /SCRNS/ nx(lfd), ny(lfd), nz(lfd), rl(lfd), ul(lfd),
     >               vl(lfd), wl(lfd), pl(lfd), tl(lfd), al(lfd),
     >               el(lfd),rr(lfd), ur(lfd), vr(lfd), wr(lfd),
     >               pr(lfd),tr(lfd), ar(lfd),er(lfd),phl(lfd),fs(lfd),
     >               jaco_f(lfd),flx(lfd,toteq),jaco_c(lx1*lz1)
      real nx, ny, nz, rl, ul, vl, wl, pl, tl, al, el, rr, ur, vr, wr,
     >                pr,tr, ar,er,phl,fs,jaco_f,flx,jaco_c

!     REAL vf(3)
      real nTol
      character*132 deathmessage
      common /nekcb/ cb
      character*3 cb

      nTol = 1.0E-14

      fdim=ldim-1
      nface = 2*ldim
      nxz   = lx1*lz1
      nxzd  = lxd*lzd
      ifield= 1 ! You need to figure out the best way of dealing with
                ! this variable

      do e=1,nelt
      do f=1,nface

! JH111218 Test 2-point fluxes in rotated form.

         call copy(nx,unx(1,1,f,e),nxz)
         call copy(ny,uny(1,1,f,e),nxz)
         call copy(nz,unz(1,1,f,e),nxz)

         call copy(rl,wminus(1,f,e,jrhof),nxz)

         if(if3d) then
            call vdot3(ul,wminus(1,f,e,jux),wminus(1,f,e,juy),
     >                    wminus(1,f,e,juz),nx,ny,nz,nxz)
            call vdot3(vl,wminus(1,f,e,jux),wminus(1,f,e,juy),
     >wminus(1,f,e,juz),t1x(1,1,f,e),t1y(1,1,f,e),t1z(1,1,f,e),nxz)
            call vdot3(wl,wminus(1,f,e,jux),wminus(1,f,e,juy),
     >wminus(1,f,e,juz),t2x(1,1,f,e),t2y(1,1,f,e),t2z(1,1,f,e),nxz)
            call vdot3(ur,wplus(1,f,e,jux),wplus(1,f,e,juy),
     >                    wplus(1,f,e,juz),nx,ny,nz,nxz)
            call vdot3(vr,wplus(1,f,e,jux),wplus(1,f,e,juy),
     >wplus(1,f,e,juz),t1x(1,1,f,e),t1y(1,1,f,e),t1z(1,1,f,e),nxz)
            call vdot3(wr,wplus(1,f,e,jux),wplus(1,f,e,juy),
     >wplus(1,f,e,juz),t2x(1,1,f,e),t2y(1,1,f,e),t2z(1,1,f,e),nxz)
         else
            call vdot2(ul,wminus(1,f,e,jux),wminus(1,f,e,juy),
     >                    nx,ny,nxz)
            call vdot2(vl,wminus(1,f,e,jux),wminus(1,f,e,juy),
     >                 t1x(1,1,f,e),t1y(1,1,f,e),nxz)
            call rzero(wl,nxz)
            call vdot2(ur,wplus(1,f,e,jux),wplus(1,f,e,juy),
     >                    nx,ny,nxz)
            call vdot2(vr,wplus(1,f,e,jux),wplus(1,f,e,juy),
     >                 t1x(1,1,f,e),t1y(1,1,f,e),nxz)
            call rzero(wr,nxz)
         endif

         call copy(pl,wminus(1,f,e,jpr),nxz)
         call copy(tl,wminus(1,f,e,jthm),nxz)
         call copy(al,wminus(1,f,e,jsnd),nxz)
         call copy(el,wminus(1,f,e,jcpf),nxz)

         call copy(rr,wplus(1,f,e,jrhof),nxz)
         call copy(pr,wplus(1,f,e,jpr),nxz)
         call copy(tr,wplus(1,f,e,jthm),nxz)
         call copy(ar,wplus(1,f,e,jsnd),nxz)
         call copy(er,wplus(1,f,e,jcpf),nxz)

         call copy(phl,wminus(1,f,e,jph),nxz)

         call copy(jaco_f,jface(1,1,f,e),nxz) 
         do i=1,nxz
            el(i)=el(i)+0.5*(ul(i)**2+vl(i)**2+wl(i)**2)
            er(i)=er(i)+0.5*(ur(i)**2+vr(i)**2+wr(i)**2)
         enddo

         call KGrotFluxFunction(nxzd,jaco_f,rl,ul,vl,wl,pl,
     >                          al,tl,rr,ur,vr,wr,pr,ar,tr,flx,el,er)

         do j=1,toteq
            call col2(flx(1,j),phl,nxzd)
         enddo
! rotate back
         call copy(flux(1,f,e,1),flx(1,1),nxz)
         if (if3d) then
            call vdot3(flux(1,f,e,2),flx(1,2),flx(1,3),flx(1,4),
     >                                 nx,t1x(1,1,f,e),t2x(1,1,f,e),nxz)
            call vdot3(flux(1,f,e,3),flx(1,2),flx(1,3),flx(1,4),
     >                                 ny,t1y(1,1,f,e),t2y(1,1,f,e),nxz)
            call vdot3(flux(1,f,e,4),flx(1,2),flx(1,3),flx(1,4),
     >                                 nz,t1z(1,1,f,e),t2z(1,1,f,e),nxz)
         else
            call vdot2(flux(1,f,e,2),flx(1,2),flx(1,3),
     >                                 nx,t1x(1,1,f,e),nxz)
            call vdot2(flux(1,f,e,3),flx(1,2),flx(1,3),
     >                                 ny,t1y(1,1,f,e),nxz)
         endif
         call copy(flux(1,f,e,5),flx(1,5),nxz)

      enddo
      enddo

      end

!-----------------------------------------------------------------------

      subroutine gtu_wrapper(fatface)
      include 'SIZE'
      real fatface(*)
      integer eq
               !                   -
      iuj=iflx ! overwritten with U -{{U}}
!-----------------------------------------------------------------------
!                          /     1  T \
! JH082316 imqqtu computes | I - -QQ  | U for all 5 conserved variables
!                          \     2    /
! which I now make the following be'neon-billboarded assumption:
!***********************************************************************
! ASSUME CONSERVED VARS U1 THROUGH U5 ARE CONTIGUOUSLY STORED
! SEQUENTIALLY IN /CMTSURFLX/ i.e. that ju2=ju1+1, etc.
! CMTDATA BETTA REFLECT THIS!!!
!***********************************************************************
C> res1+=\f$\int_{\Gamma} \{\{\mathbf{A}^{\intercal}\nabla v\}\} \cdot \left[\mathbf{U}\right] dA\f$
! JH070918 conserved variables done here.
      i_cvars=(ju1-1)*nfq+1
      do eq=1,toteq
         call faceu(eq,fatface(i_cvars))
! JH080317 at least get the product rule right until we figure out how
!          we want the governing equations to look
         call invcol2(fatface(i_cvars),fatface(iwm+nfq*(jph-1)),nfq)
         i_cvars=i_cvars+nfq
      enddo
      ium=(ju1-1)*nfq+iwm
      iup=(ju1-1)*nfq+iwp
      call   imqqtu(fatface(iuj),fatface(ium),fatface(iup))
      call   imqqtu_dirichlet(fatface(iuj),fatface(iwm),fatface(iwp))

      if (1 .eq. 2) then
      call igtu_cmt(fatface(iwm),fatface(iuj),graduf) ! [[u]].{{gradv}}
C> res1+=\f$\int \left(\nabla v\right) \cdot \left(\mathbf{H}^c+\mathbf{H}^d\right)dV\f$ 
C> for each equation (inner), one element at a time (outer)
      endif
      return
      end
