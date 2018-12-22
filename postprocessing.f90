!====================================================================
      SUBROUTINE DUDTPOSTPROC(kkx,kky,kkz,ux,uy,uz,uold,dir,id)

! Calculations of fluid acceleration Du/Dt
      USE params
!      IMPLICIT NONE
      INCLUDE 'mpif.h'

      REAL,DIMENSION      (llx,lly,lz)  ::  ux,uy,uz,uold
      REAL,DIMENSION      (llx,lly,lz)  ::  wx1,wy1,wz1,acce
      REAL,DIMENSION      (llx,lly,lz)  ::  kkx,kky,kkz
      REAL,DIMENSION      (llx,lly,lz)  ::  tmph
      REAL				::  pmin1,pmin,pmax1,pmax, &
					    dudt1,dudt2,dudt1_t,   &
					    dudt2_t
      REAL                              ::  lpmax(3),co1(3),co1_t(3)
      INTEGER				::  id,dir,ierr

      if(dir.eq.1) tmph=ux
      if(dir.eq.2) tmph=uy
      if(dir.eq.3) tmph=uz
      wx1 = kkx * tmph
      wy1 = kky * tmph
      wz1 = kkz * tmph
      acce = wx1
      wx1(1:lx , 1:lly-1:2 , :) = -acce(1:lx , 2:lly:2   , :)
      wx1(1:lx , 2:lly:2   , :) =  acce(1:lx , 1:lly-1:2 , :)
      acce = wy1
      wy1(1:lx , 1:lly-1:2 , :) = -acce(1:lx , 2:lly:2   , :)
      wy1(1:lx , 2:lly:2   , :) =  acce(1:lx , 1:lly-1:2 , :)
      acce = wz1
      wz1(1:lx , 1:lly-1:2 , :) = -acce(1:lx , 2:lly:2   , :)
      wz1(1:lx , 2:lly:2   , :) =  acce(1:lx , 1:lly-1:2 , :)
      acce = (tmph-uold)/dt + ux*wx1 + uy*wy1 + uz*wz1
      CALL mpifft3DCR (acce,id)

      dudt1 = sum(acce(1:lx,1:ly,:))
      dudt2 = sum((acce(1:lx,1:ly,:))**2)
      CALL MPI_REDUCE (dudt1,dudt1_t,1,MPI_REAL8,MPI_SUM,0,        &
                                   mpi_comm_world,ierr)
      CALL MPI_REDUCE (dudt2,dudt2_t,1,MPI_REAL8,MPI_SUM,0,        &
                                   mpi_comm_world,ierr)
      dudt1_t = dudt1_t / real(nx*ny*nz)
      dudt2_t = sqrt ( (dudt2_t-dudt1_t**2)/real(nx*ny*nz) )

      pmax1=MAXVAL(abs(acce(1:lx,1:ly,:)))
      pmin1=MINVAL(abs(acce(1:lx,1:ly,:)))
      CALL MPI_REDUCE(pmax1,pmax,1,MPI_REAL8,MPI_MAX,0, &
                                   mpi_comm_world,ierr)
      CALL MPI_REDUCE(pmin1,pmin,1,MPI_REAL8,MPI_MIN,0, &
                                   mpi_comm_world,ierr)

      if(id.eq.0) then
        lpmax(1) = 0.9 * pmax
        lpmax(2) = 0.8 * pmax
        lpmax(3) = 0.7 * pmax
      endif
      CALL MPI_BCAST (lpmax,3,MPI_REAL8,0,mpi_comm_world,ierr)
      co1(1) = float( count(abs(acce(1:lx,1:ly,:)).gt.lpmax(1)) )
      co1(2) = float( count(abs(acce(1:lx,1:ly,:)).gt.lpmax(2)) )
      co1(3) = float( count(abs(acce(1:lx,1:ly,:)).gt.lpmax(3)) )
      CALL MPI_ALLREDUCE (co1,co1_t,3,MPI_REAL8,MPI_SUM,&
                                   mpi_comm_world,ierr)

      if(id.eq.0) write(42+dir,601) time,pmin,pmax,dudt1_t,           &
                                    dudt2_t,co1_t(1)/float(nx*ny*nz), &
                                    co1_t(2)/float(nx*ny*nz),         &
                                    co1_t(3)/float(nx*ny*nz)

601   format ( 2x,8(1pe15.5))

      RETURN
      END SUBROUTINE DUDTPOSTPROC

!==================================================================
      SUBROUTINE POSTPROC(ux,uy,uz,uxr,uyr,uzr,wxc,wyc,wzc, &
                          kkx,kky,kkz,k2l,ik2l,rnul,id)

      USE params
!      IMPLICIT NONE
      INCLUDE 'mpif.h'

      REAL,DIMENSION     (llx,lly,lz)  ::  ux,uy,uz
      REAL,DIMENSION     (llx,lly,lz)  ::  uxr,uyr,uzr
      REAL,DIMENSION     (llx,lly,lz)  ::  wxc,wyc,wzc
      REAL,DIMENSION     (llx,lly,lz)  ::  wx1,wy1,wz1,vdiv
      REAL,DIMENSION     (llx,lly,lz)  ::  tmph,tmp1h,k2l
      INTEGER,DIMENSION  (llx,lly,lz)  ::  ik2l
      REAL,DIMENSION     (llx,lly,lz)  ::  kkx,kky,kkz
      REAL			       ::  a3,a3t,a2,a2t,vskew,vflat
      REAL			       ::  rnul,XINTLS,dissp,qq,CFL, &
					   RESOL,REE,TMSE,XL,tk,eta, &
					   UPRIME,vk
      REAL                             ::  e_t,e_ts,eks,ek,et,pmax,  &
					   pmax1,dx
      REAL			       ::  vdiv2,vdiv1,vort2_t,vort1_t, &
					   divvt_t,divvt,vdiv2_t,       &
					   vdiv1_t,vort2,vort1
      INTEGER			       ::  i,kkk,id,ierr

! Calculations of skewness and flatness of velocity gradient
      wx1 = kkx * ux
      wy1 = kky * uy
      wz1 = kkz * uz
      tmph = wx1
      wx1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wy1
      wy1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wy1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wz1
      wz1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wz1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)

!     Parallel FFT
      CALL mpifft3DCR (wx1,id)
      CALL mpifft3DCR (wy1,id)
      CALL mpifft3DCR (wz1,id)

      tmph = wx1**2 + wy1**2 + wz1**2
      tmph = tmph / 3.
      a2   = sum ( tmph(1:lx,1:ly,:) )
      CALL MPI_REDUCE (a2,a2t,1,MPI_REAL8,MPI_SUM,0,     &
                                mpi_comm_world,ierr)

      tmph  = wx1**3 + wy1**3 + wz1**3
      tmph  = tmph / 3.
      a3   = sum ( tmph(1:lx,1:ly,:) )
      CALL MPI_REDUCE (a3,a3t,1,MPI_REAL8,MPI_SUM,0,     &
                                mpi_comm_world,ierr)

      vskew = a3t / (nx*ny*nz) / (a2t/(nx*ny*nz))**1.5

      tmph = wx1**4 + wy1**4 + wz1**4
      tmph = tmph / 3.
      a3   = sum ( tmph(1:lx,1:ly,:) )
      CALL MPI_REDUCE (a3,a3t,1,MPI_REAL8,MPI_SUM,0,     &
                                mpi_comm_world,ierr)

      vflat = a3t / (nx*ny*nz) / (a2t/(nx*ny*nz))**2

! Calculations of energy and dissipation spectra. And other
! important statistics. 

      if ( id.eq.0 ) then
        XINTLS = 0.0
        dissp  = 0.0
        qq     = 0.0
      endif

      tmph = ux*ux + uy*uy + uz*uz
      if(idlocalY.eq.0) then
        tmph(:,1,:) = 0.5 * tmph(:,1,:)
        tmph(:,2,:) = 0.5 * tmph(:,2,:)
      endif
      tmp1h = 2. * rnul * tmph * k2l

!Write out the total energy in K space
      ek = sum(tmph(1:lx,:,:))
      CALL MPI_REDUCE (ek,e_t,1,MPI_REAL8,MPI_SUM,0,      &
                                 mpi_comm_world,ierr)
!      if ( id.eq.0 ) write(*,*) jstep, time, e_t

      if ( id.eq.0 ) then
        if ( mod(jstep,ispec).eq.0 ) then
          kkk = jstep/ispec - 1
          write (20,*) kkk - 1
        endif
      end if

      do i = 1, nek
        ek  = 0.0
        e_t = 0.0
        ek  = sum ( tmph(1:lx,:,:), mask=(ik2l(1:lx,:,:).eq.i) )
        CALL MPI_REDUCE (ek,e_t,1,MPI_REAL8,MPI_SUM,0,       &
                                     mpi_comm_world,ierr)
        eks  = 0.0
        e_ts = 0.0
        eks  = sum ( tmp1h(1:lx,:,:), mask=(ik2l(1:lx,:,:).eq.i) )
        CALL MPI_REDUCE (eks,e_ts,1,MPI_REAL8,MPI_SUM,0,     &
                                     mpi_comm_world,ierr)
        if ( id.eq.0 ) then
          if ( mod(jstep,ispec).eq.0 ) write(20,*) i, e_t, e_ts
          qq     = qq + e_t
          dissp  = dissp + e_ts
          XINTLS = XINTLS + e_t/real(i)
        end if
      end do

      PI2   = 8.*atan(1.)
      if ( id.eq.0 ) then
        eta    = ( rnul**3 / dissp )**.25
        tk     = ( rnul / dissp )**.5
        vk     = ( rnul * dissp )**.25
        UPRIME = sqrt ( 2.0 * qq / 3. )
        TMSE   = sqrt ( 15. * rnul * UPRIME**2 / dissp )
        REE    = UPRIME * TMSE / rnul
        XL     = UPRIME**3 / dissp
        ET     = XL / UPRIME
        RESOL  = real(akmax) * ETA
        XINTLS = XINTLS * (PI2/2.) / 2. / UPRIME**2
      endif

      DX    = PI2 / REAL(ny)
      tmph  = abs(uxr) + abs(uyr) + abs(uzr)
      pmax1 = MAXVAL (TMPH)
      CALL MPI_REDUCE(pmax1,pmax,1,MPI_REAL8,MPI_MAX,0, &
                                   mpi_comm_world,ierr)
      if ( id.eq.0 ) then
        CFL = pmax * dt / dx
        if ( jstep.eq.0 ) CFL = pmax * dt_h / dx
        write(40,601) time,eta,TMSE,XINTLS,tk,ET,XL,UPRIME
        write(41,601) time,qq,dissp,REE,RESOL,CFL,vskew,vflat
        write(*,*) 'REE,CFL,RESOL,vskew,vflat',REE,CFL,RESOL,   &
                     vskew,vflat
      end if
601   format ( 2x,8(1pe15.5) )

       if ( divstat.eqv..TRUE.) then
         vdiv  = kkx*real(ux) + kky*real(uy) + kkz*real(uz)
         tmph = vdiv
         vdiv(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
         vdiv(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
         CALL mpifft3DCR (vdiv,id)

!-------------------------------------------------------------------
!       compute the mean and rms of vdiv
        vdiv1 = 0.0
        vdiv2 = 0.0
        vdiv1 = sum (vdiv(1:lx,1:ly,:)) 
        vdiv2 = sum ((vdiv(1:lx,1:ly,:))**2)
        CALL MPI_REDUCE (vdiv1,vdiv1_t,1,MPI_REAL8,MPI_SUM,0,        &
                                      mpi_comm_world,ierr)
        CALL MPI_REDUCE (vdiv2,vdiv2_t,1,MPI_REAL8,MPI_SUM,0,        &
                                      mpi_comm_world,ierr)
        vdiv1_t = vdiv1_t / real(nx*ny*nz)
        vdiv2_t = sqrt ( (vdiv2_t-vdiv1_t**2)/real(nx*ny*nz) )

!-------------------------------------------------------------------
!       compute the mean and rms of vorticity
        vort1 = 0.0
        vort2 = 0.0
        tmph   = sqrt( (wxc)**2 + (wyc)**2 +  (wzc)**2 )
        vort1 = sum ( tmph(1:lx,1:ly,:) ) 
        vort2 = sum ( tmph(1:lx,1:ly,:)**2 )
        CALL MPI_REDUCE (vort1,vort1_t,1,MPI_REAL8,MPI_SUM,0,        &
                                       mpi_comm_world,ierr)
        CALL MPI_REDUCE (vort2,vort2_t,1,MPI_REAL8,MPI_SUM,0,        &
                                       mpi_comm_world,ierr)
        vort1_t = vort1_t / real(nx*ny*nz)
        vort2_t = sqrt ( (vort2_t-vort1_t**2) / real(nx*ny*nz) )

!-------------------------------------------------------------------
!       cross correlation
        divvt = 0.0
        tmph   = tmph*vdiv
        divvt = sum ( tmph(1:lx,1:ly,:))
        CALL MPI_REDUCE (divvt,divvt_t,1,MPI_REAL8,MPI_SUM,0,        &
                                       mpi_comm_world,ierr)
        divvt_t = (divvt_t - vdiv1_t*vort1_t) / real(2*nx*ny*nz)

        if( id.eq.0 ) then
          write (42,602) time, vdiv1_t, vdiv2_t, vort1_t, vort2_t, divvt_t
        end if

       endif

602    format ( 2x,7(1pe15.5) )

      RETURN
      END SUBROUTINE POSTPROC

!====================================================================
      SUBROUTINE TRANSPOSTPROC(ux,uy,uz,uwx,uwy,uwz,kk2,id)

! Energy transfer calculations
      USE params
!      IMPLICIT NONE
      INCLUDE 'mpif.h'

      REAL,DIMENSION     (llx,lly,lz)  ::  ux,uy,uz
      REAL,DIMENSION     (llx,lly,lz)  ::  uwx,uwy,uwz
!      REAL,DIMENSION     (llx,lly,lz)  ::  kk2

      INTEGER,DIMENSION     (llx,lly,lz)  ::  kk2

      REAL,DIMENSION     (llx,lly,lz)  ::  tmph
      INTEGER			       ::  i,id,jjj,ierr
      REAL			       ::  ek,e_t

      tmph = ux*uwx + uy*uwy + uz*uwz
      if(idlocalY.eq.0) then
        tmph(:,1,:) = 0.5 * tmph(:,1,:)
        tmph(:,2,:) = 0.5 * tmph(:,2,:)
      endif

      if ( id.eq.0 ) then
        jjj = jstep/ispec - 1
        write(22,*) jjj
      endif

      do i = 1, nek
        ek  = 0.
        e_t = 0.
!        ek = sum ( tmph(1:lx,:,:),mask=(abs(sqrt(kk2(1:lx,:,:)) &
!                                               -i-0.49999).lt.0.5))
        ek  = sum ( tmph(1:lx,:,:), mask=(kk2(1:lx,:,:).eq.i) )


        CALL MPI_REDUCE (ek,e_t,1,MPI_REAL8,MPI_SUM,0,mpi_comm_world,ierr)
        if ( id.eq.0 ) write(22,*) i, e_t
      enddo

      RETURN
      END SUBROUTINE TRANSPOSTPROC

!==================================================================

