!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE wavenumber (kx,ky,kz,k2,id)

      USE params
!      IMPLICIT NONE

      REAL, DIMENSION (llx,lly,lz)  :: kx,ky,kz,k2
      INTEGER			    :: i,j,k,ig,jg
      INTEGER                       :: id

      do i = 2,(ly+lyext),2
        ig = (ly*idlocalY) + i
        kx(:,i-1,:) = (ig/2) - 1
        kx(:,i,:) = (ig/2) - 1
      enddo

      do j = 1, lz
       jg = (lz*idlocalZ)+j
       if (jg.lt.ny/2+2) then
        ky(:,:,j) = jg - 1
       else
        ky(:,:,j) = -(ny+1-jg)
       endif
      enddo

      do k = 1, lx
       if (k.lt.nz/2+2) then
        kz(k,:,:) = k - 1
       else
        kz(k,:,:) = -(nz+1-k)
       endif
      enddo

      k2   = kx*kx + ky*ky + kz*kz
      if ( id.eq.0 ) then
        k2(1,1:2,1)   = 1.e-5
      endif

! Putting zeros (or ones for k2) at the additional padding on the x_matrix 
! direction. This padding is needed for the physical space calculation but
! not on the spectral space case
      kx(lx+1:llx,:,:) = 0.0
      ky(lx+1:llx,:,:) = 0.0
      kz(lx+1:llx,:,:) = 0.0
      k2(lx+1:llx,:,:) = 1.0

      RETURN
      END SUBROUTINE wavenumber

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE initial_flow (ux,uy,uz,iseed_m0,ik2l,enuml,  &
                               kkx,kky,kkz,kk2,force,id)

!      this line corrects irand functionality for intel compiler.    
!      USE IFPORT   
      USE params
!      IMPLICIT NONE
      INCLUDE 'mpif.h' 

      REAL, DIMENSION (llx,lly,lz)  :: ux,uy,uz,tmp
      REAL, DIMENSION (llx,lly,lz)  :: kkx,kky,kkz,kk2
      REAL, DIMENSION (nek)         :: enuml
      REAL                          :: force(2),ek,e_t,vol,e1,akp,u0,ccc

      INTEGER, DIMENSION (llx,lly,lz)  :: ik2l
      INTEGER                          :: iseed_m0,iseed(nproc)
      INTEGER                          :: i,id,ierr,irr,irand

        pi  = 4.*atan(1.)
! -----------------------------------------------------------------
!       Generate random number seed for each process
 
        if ( id.eq.0 ) then
!          CALL srand (iseed_m0)
          do i = 1, nproc
!            iseed(i) = irand()
            iseed(i) = 24 
          end do
        end if

        CALL MPI_BCAST (iseed,nproc,MPI_INTEGER,0,mpi_comm_world,ierr)
        irr = iseed(id+1)
        CALL srand (irr)

! ----------------------------------------------------------------
!       Pre-initial Gaussian flow field

        CALL gaussian (ux)
        CALL gaussian (uy)
        CALL gaussian (uz)

        CALL symmetrize (ux,id)
        CALL symmetrize (uy,id)
        CALL symmetrize (uz,id)

        tmp = ( kkx*ux + kky*uy + kkz*uz ) / kk2
        ux  = ux - kkx*tmp
        uy  = uy - kky*tmp
        uz  = uz - kkz*tmp

! -----------------------------------------------------------------
!       scale the spectrum to a specific shape
!       (1) K^-5/3 with a fixed E(1)

        tmp = ux*ux + uy*uy + uz*uz
        if(idlocalY.eq.0) then
           tmp(:,1,:) = 0.5 * tmp(:,1,:)
           tmp(:,2,:) = 0.5 * tmp(:,2,:)
        endif

        do i = 1, nek
          ek  = 0.0
          e_t = 0.0
          ek  = sum ( tmp, mask=(ik2l.eq.i) )
          CALL MPI_ALLREDUCE (ek,e_t,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

          vol = 4./3. * pi * ( (float(i)+0.5)**3 - (float(i)-0.5)**3 )

!         (1) K^-5/3 with a fixed E(1)
          ccc = force(1) * 4./3. * pi * ( 1.5**3 - 0.5**3 ) / enuml(1)
          e1  = sqrt ( ccc/e_t*enuml(i)/vol ) / ( real(i) )**0.833333

!          (2) Gaussian spectrum
!          akp   = 4.75683
!          u0    = 1.0
!          ccc = 16.*sqrt(2./pi)/akp
!          ccc = ccc*4./3.*pi*((real(i)+0.5)**3-(real(i)-0.5)**3)/enuml(i)
!          e1=u0*sqrt(ccc/e_t)*(real(i)/akp)**2*exp(-(real(i)/akp)**2)

          where ( ik2l.eq.i )
            ux = ux * e1
            uy = uy * e1
            uz = uz * e1
          end where
        enddo

        CALL dealiasing (ux,ik2l)
        CALL dealiasing (uy,ik2l)
        CALL dealiasing (uz,ik2l)

      RETURN
      END SUBROUTINE initial_flow

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE gaussian (u)

      USE IFPORT 
      USE params
!      IMPLICIT NONE

      REAL, DIMENSION (llx,lly,lz) :: u
      REAL                         :: t1,t2
      INTEGER                      :: i,j,k

      u = 0.0

      do i = 1, lx
       do j = 1, lly,2
        do k = 1, lz
          t1 = rand()
          t2 = rand()
          if ( t1.le.1.e-10 ) t1 = 1.e-10
          if ( t2.le.1.e-10 ) t2 = 1.e-10
          t2 = 2*pi*t2
          u(i,j,k)   = sqrt ( -2.0 * log(t1)) * cos(t2)
          u(i,j+1,k) = sqrt ( -2.0 * log(t1)) * sin(t2)
        end do
       end do
      end do

      RETURN
      END SUBROUTINE gaussian

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
!     2/3 rule

      SUBROUTINE dealiasing (vx,ik2)

      USE params
!      IMPLICIT NONE

      REAL,    DIMENSION (llx,lly,lz)   :: vx
      INTEGER, DIMENSION (llx,lly,lz)   :: ik2

      where ( ik2.gt.nek ) vx = 0.0

      RETURN
      END SUBROUTINE dealiasing

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE supforDET (vx,vy,vz,k2,ik2,force,id)

      USE params
!      IMPLICIT NONE
      INCLUDE 'mpif.h'

      REAL,    DIMENSION (llx,lly,lz)   :: vx, vy, vz
      REAL,    DIMENSION (llx,lly,lz)   :: k2, tmp
      INTEGER, DIMENSION (llx,lly,lz)   :: ik2
      REAL                              :: force(2), ff, ek
      INTEGER			        :: is, id, ierr

      tmp = vx*vx + vy*vy + vz*vz
      if(idlocalY.eq.0) then
         tmp(:,1,:) = 0.5 * tmp(:,1,:)
         tmp(:,2,:) = 0.5 * tmp(:,2,:)
      endif

      DO is = 1, 2
        ek    = 0.0
        ff    = 0.0
        ek = sum ( tmp, mask=(ik2.eq.is) )
        CALL MPI_ALLREDUCE (ek,ff,1,MPI_REAL8,MPI_SUM,       &
                                    mpi_comm_world,ierr)

        ff = sqrt ( force(is) / ff )
        where ( ik2.eq.is )
           vx = vx*ff
           vy = vy*ff
           vz = vz*ff
        end where
      END DO

      RETURN
      END SUBROUTINE supforDET

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE symmetrize(c,id)

! For this subroutine, I'm assuming 2^n processors along Z_matrix direction

      USE params
!      IMPLICIT NONE
      INCLUDE 'mpif.h'

      REAL,    DIMENSION (llx,lly,lz)   :: c
      REAL,    DIMENSION (llx,2,2:lz)   :: ccS, ccR
      REAL,    DIMENSION (llx,2)        :: ccvS, ccvR
      INTEGER                           :: id, ierr, npSR, isize
      INTEGER                           :: i, ii, k, kk
      INTEGER                           :: status(MPI_STATUS_SIZE,2),req(2)

! Force a zero mean flow
      if (id.eq.0) c(1,1:2,1) = 0.0


! Truncate the n/2 wave number as this is noise usually
! For kx = nx/2
      if(idlocalY.eq.(nprocY-1)) c(:,lly-1:lly,:) = 0.0
! For ky = ny/2
      if(idlocalZ.eq.(nprocZ/2)) c(:,:,1) = 0.0
! For kz = nz/2
      c(nz/2+1,:,:) = 0.0


! Round and other errors are corrected to ensure the physical data is real
! valued. For that, symmetry must be checked and enforced at the plane kx=0. 
! For the other planes, kx>0, the symmetry is warrantied by the FFT.
      IF(idlocalY.eq.0) THEN
        IF(nproc.ne.1) THEN

          isize = 2*llx*(lz-1)
          ccS = c(:,1:2,2:lz)
          npSR = nprocZ - idlocalZ - 1
          npSR = npSR * nprocY
          call mpi_isend(ccS,isize,MPI_REAL8,npSR,1,MPI_COMM_WORLD, &
                         req(1),ierr)
          call mpi_irecv(ccR,isize,MPI_REAL8,npSR,1,MPI_COMM_WORLD, &
                         req(2),ierr)
          call mpi_waitall(2,req,status,ierr)
          do i=1,lx
            ii=nz+2-i
            if (i.eq.1) ii=1 
            do k=2,lz 
              kk = lz+2-k
              c(i,1,k) = 0.5*(c(i,1,k) + ccR(ii,1,kk)) 
              c(i,2,k) = 0.5*(c(i,2,k) - ccR(ii,2,kk))
            enddo
          enddo
    
          if ((idlocalZ.ne.0).and.(idlocalZ.ne.nprocZ/2)) then 
            isize = llx*2
            ccvS = c(:,1:2,1)
            npSR = nprocZ - idlocalZ
            npSR = npSR * nprocY
            call mpi_isend(ccvS,isize,MPI_REAL8,npSR,2,MPI_COMM_WORLD, &
                         req(1),ierr)
            call mpi_irecv(ccvR,isize,MPI_REAL8,npSR,2,MPI_COMM_WORLD, &
                         req(2),ierr)
            call mpi_waitall(2,req,status,ierr)

            do i=1,lx
              ii=nz+2-i
              if (i.eq.1) ii=1
              c(i,1,1) = 0.5*(c(i,1,1) + ccvR(ii,1))
              c(i,2,1) = 0.5*(c(i,2,1) - ccvR(ii,2))
            enddo
          endif

        ELSE

          do i=1,lx
            ii=nz+2-i
            if (i.eq.1) ii=1
            do k=2,lz/2
              kk = ny+2-k
              c(i,1,k) = 0.5*(c(i,1,k) + c(ii,1,kk))
              c(i,2,k) = 0.5*(c(i,2,k) - c(ii,2,kk))
              c(ii,1,kk) =  c(i,1,k)
              c(ii,2,kk) = -c(i,2,k)
            enddo
          enddo

        ENDIF

          if (idlocalZ.eq.0) then
            do i=2,lx/2
              ii=nz+2-i
              c(i,1,1) = 0.5*(c(i,1,1) + c(ii,1,1))
              c(i,2,1) = 0.5*(c(i,2,1) - c(ii,2,1))
              c(ii,1,1) =  c(i,1,1)
              c(ii,2,1) = -c(i,2,1)
            enddo
          endif

      ENDIF

      RETURN
      END SUBROUTINE symmetrize

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE SETWAVERAND (kxr,kyr,kzr)

! This subroutine is to set up the wave number values to be used on the 
! sthocastic forcing scheme. This is done as the this forcing scheme is perform
! on ALL processors. They will be calculating same thing but that way no 
! processor is idle. Thus, ALL processors need all wave number values as they
! can not use their local wave numbers. Forcing is a quick calculation, hard to 
! parallelize.

      USE params
!      IMPLICIT NONE

      INTEGER 			      ::  i,j,k
      REAL			      ::  kxr(nx+2),kyr(ny),kzr(nz)

      do i = 1, nx+2, 2
       kxr(i)   = int(i/2)
       kxr(i+1) = kxr(i)
      end do

      do j = 1, ny
       if ( j.lt.ny/2+2 ) then
        kyr(j) = j - 1
       else
        kyr(j) = -(ny+1-j)
       endif
      end do

      do k = 1, nz
       if ( k.lt.nz/2+2 ) then
        kzr(k) = k - 1
       else
        kzr(k) = -(nz+1-k)
       endif
      end do

      END SUBROUTINE SETWAVERAND

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE supforSTH (ux,uy,uz,iseedf,ivf,iyf,kxr,kyr,kzr,  &
                            a1r,a2r,a3r,b1r,b2r,b3r,id)

      USE params
!      IMPLICIT NONE
      INCLUDE 'mpif.h'

      REAL, DIMENSION (llx,lly,lz)   :: ux, uy, uz
      REAL      		     :: fr(6,5,5,3),kxr(nx+2),kyr(ny),kzr(nz)
      REAL                           :: a1r(6,5,5),a2r(6,5,5),a3r(6,5,5), &
					b1r(6,5,5),b2r(6,5,5),b3r(6,5,5)
      INTEGER                        :: k,kk,j,jj,i,ii,mc,ip,id
      INTEGER   		     :: iseedf, iyf, ivf(32)
      REAL      		     :: cont1, rmag, rph, divb, radm, smallx,  &
					xkf,var,tf,t1,t2,t3,ranff

!     Use this subroutine to time step the random forcing and...
!     generate the random, incompressible, real-valued body...
!     force (a1r,a2r,a3r) for the next time level.
!     Calculations are in spectral space.

      smallx = 1.0e-18
      xkf = sqrt(8.0)
      var = 447.31
      tf  = 0.038

!     Generate white noise
      do mc=1,3
        do  k=1,5
          do  j=1,5
            do  i=1,6,2
              ip=i+1
              cont1 = ranff(iseedf,iyf,ivf)
              if ( cont1.lt.1.e-6 ) cont1 = 1.e-6
              rmag = sqrt ( -4.0 * var * dt * alog(cont1) / tf )
              rph =  ranff (iseedf,iyf,ivf)
              fr(i ,j,k,mc) = rmag * cos(pi2*rph)
              fr(ip,j,k,mc) = rmag * sin(pi2*rph)
            end do
          end do
        end do
      end do

!     Overwrite components for k1=0 to ensure real forcing
!     in physical space with no spatial mean.

      do mc = 1, 3
        do k = 1, 5
          kk = 7 - k
          if ( k.eq.1 ) kk = 1
          do j = 2, 3
            jj = 7 - j
            fr(1,jj,kk,mc) =  fr(1,j,k,mc)
            fr(2,jj,kk,mc) = -fr(2,j,k,mc)
          end do
        end do
        do k = 2, 3
          kk = 7 - k
          fr(1,1,kk,mc) =  fr(1,1,k,mc)
          fr(2,1,kk,mc) = -fr(2,1,k,mc)
        end do
        fr(1:2,1,1,mc)  = 0.0
      end do

!     Euler step for Langevin type process with the...
!     noise forcing. Then obtain incompressible part...
!     and save as (a1r,a2r,a3r). Possible wave nos. are...
!     k1=0,1,2 and k2,k3 = 0,1,2,-2,-1 .

      DO k = 1, 5
        kk = k
        if ( k.gt.3 ) kk = k + nz - 5
        DO j = 1, 5
          jj = j
          if ( j.gt.3 ) jj = j + ny - 5
          DO i = 1, 6
            radm = sqrt ( kxr(i)**2 + kyr(jj)**2 + kzr(kk)**2 )
            IF ( radm.lt.xkf ) THEN
              b1r(i,j,k) = b1r(i,j,k) * (1.0-dt/tf) + fr(i,j,k,1)
              b2r(i,j,k) = b2r(i,j,k) * (1.0-dt/tf) + fr(i,j,k,2)
              b3r(i,j,k) = b3r(i,j,k) * (1.0-dt/tf) + fr(i,j,k,3)

              t1 = kxr(i)  * b1r(i,j,k)
              t2 = kyr(jj) * b2r(i,j,k)
              t3 = kzr(kk) * b3r(i,j,k)

              divb = (t1+t2+t3) / ( (radm+smallx)**2 )

              if ( radm.eq.0. ) divb = 0.0

              a1r(i,j,k) = b1r(i,j,k) - kxr(i)  * divb
              a2r(i,j,k) = b2r(i,j,k) - kyr(jj) * divb
              a3r(i,j,k) = b3r(i,j,k) - kzr(kk) * divb
            ELSE
              a1r(i,j,k) = 0.0
              a2r(i,j,k) = 0.0
              a3r(i,j,k) = 0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!     Eliminate forcing of a mean flow at zero wave no.
      a1r(1:2,1,1) = 0.0
      a2r(1:2,1,1) = 0.0
      a3r(1:2,1,1) = 0.0

!     Ensure that the forcing is real-valued. This step...
!     should only clean up rounding errors if fr is ok.

      DO k = 1, 5
        kk = 7 - k
        if ( k.eq.1 ) kk = 1
        do j = 2, 3
          jj = 7 - j
          a1r(1,jj,kk) =  a1r(1,j,k)
          a1r(2,jj,kk) = -a1r(2,j,k)
          a2r(1,jj,kk) =  a2r(1,j,k)
          a2r(2,jj,kk) = -a2r(2,j,k)
          a3r(1,jj,kk) =  a3r(1,j,k)
          a3r(2,jj,kk) = -a3r(2,j,k)
        enddo
      ENDDO

      do k = 2, 3
        kk = 7 - k
        a1r(1,1,kk) =  a1r(1,1,k)
        a1r(2,1,kk) = -a1r(2,1,k)
        a2r(1,1,kk) =  a2r(1,1,k)
        a2r(2,1,kk) = -a2r(2,1,k)
        a3r(1,1,kk) =  a3r(1,1,k)
        a3r(2,1,kk) = -a3r(2,1,k)
      enddo

! NOW, let's add the forcing term to advance the velocity
! at the substep using Euler scheme which is effective.
! The forcing is added to the spectral domain velocity only
! in a few mesh nodes where the wave number is lower than
! SQRT(8). Those nodes are 3 at the begining and 2 at the end
! of the kz nd ky directions, and 6 at the begining of kx
! direction (considering that the complex pair are saved
! one next to the other.
! THIS SUBROUTINE WILL ASSUME THAT THE DOMAIN DECOMPOSITION
! WAS DONE SUCH AS THE MENTIONED NODES ARE INSIDE ONLY ONE 
! PROCESSOR. This means that the number of processors along
! the Y-direction should be: nprocY > nx/6; and  the number
! of processors along Z-direction should be: nprocZ > nx/3.

! NOTE: a1r,a2r & a3r are based on a1r(kx,ky,kz) while velocity
! is based on u(kz,kx,ky)

      IF (idlocalY.eq.0) THEN
        IF (idlocalZ.eq.0) THEN
          do k=1,5
            kk = k
            if(k.gt.3) kk = k+nz-5
            do j=1,3
              jj = j
              do i=1,6
                ux(kk,i,jj) = ux(kk,i,jj) + dt*a1r(i,j,k)
                uy(kk,i,jj) = uy(kk,i,jj) + dt*a2r(i,j,k)
                uz(kk,i,jj) = uz(kk,i,jj) + dt*a3r(i,j,k)
              enddo
            enddo
          enddo
        ENDIF

        IF (idlocalZ.eq.(nprocZ-1)) THEN
          do k=1,5
            kk = k
            if(k.gt.3) kk = k+nz-5
            do j=4,5
              jj = j+lz-5
              do i=1,6
                ux(kk,i,jj) = ux(kk,i,jj) + dt*a1r(i,j,k)
                uy(kk,i,jj) = uy(kk,i,jj) + dt*a2r(i,j,k)
                uz(kk,i,jj) = uz(kk,i,jj) + dt*a3r(i,j,k)
              enddo
            enddo
          enddo
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE supforSTH

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

