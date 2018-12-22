!========================================================================
!========================================================================
!   Subroutine for MPI parallel FFT scheme of Ayala (2010)
!   based on a 2D decomposition of the domain 
!   We used:
!
!           MPICH (portable MPI) version 1.2.0 
!           www-unix.mcs.anl.gov/mpi/mpich
!
!           Personal version of Fortran 90 (V3.4N5), from 
!           Pacific Sierra Research www.psrv.com
!
!           (MPI and F90 are actually libraries based on GNU g77 native 
!           Linux Fortran 77 compiler)
!
!           FFTW Fastest Fourier Transform in the West version 2.1.3
!           www.fftw.org
!
!   This version uses a coupling communication scheme (see mpitranspose
!   subroutine). The processors communicate in a cyclic permutation scheme
!
!========================================================================
!========================================================================
!   Set FFTW plans
!   Plans are passed to calling programs through COMMON in
!   included file 'plans.h'
!=======================================================================

      SUBROUTINE set_plan1D_RC (nx)

      INCLUDE 'fftw_include.h'
      INCLUDE 'plans.h'

      INTEGER, Intent (in)     :: nx

      call rfftwnd_f77_create_plan_(plan_RC,1,nx, &
                FFTW_REAL_TO_COMPLEX,FFTW_MEASURE+FFTW_IN_PLACE)

      call rfftwnd_f77_create_plan_(plan_CR,1,nx, &
                FFTW_COMPLEX_TO_REAL,FFTW_MEASURE+FFTW_IN_PLACE)

      END SUBROUTINE set_plan1D_RC
!=======================================================================

      SUBROUTINE set_plan1D_CC (ny)

      INCLUDE 'fftw_include.h'
      INCLUDE 'plans.h'

      INTEGER, Intent (in)     :: ny

      call fftw_f77_create_plan_(plan_F,ny,       &
                 FFTW_FORWARD,FFTW_MEASURE+FFTW_OUT_OF_PLACE)

      call fftw_f77_create_plan_(plan_B,ny,       &
                 FFTW_BACKWARD,FFTW_MEASURE+FFTW_OUT_OF_PLACE)

      END SUBROUTINE set_plan1D_CC
!=======================================================================

      SUBROUTINE destroy_plan1D_RC

       INCLUDE 'plans.h'

       call rfftwnd_f77_destroy_plan_(plan_RC)
       call rfftwnd_f77_destroy_plan_(plan_CR)

       END SUBROUTINE destroy_plan1D_RC
!=======================================================================

      SUBROUTINE destroy_plan1D_CC

      INCLUDE 'plans.h'

      call fftw_f77_destroy_plan_(plan_F)
      call fftw_f77_destroy_plan_(plan_B)

      END SUBROUTINE destroy_plan1D_CC
!=======================================================================

      SUBROUTINE  mpifft3DRC (v,id)

      USE params

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'plans.h'

      INTEGER, Intent (in)             :: id
      REAL,    Intent (inout)          :: v(lx+2,ly+lyext,lz)
      COMPLEX                          :: vaux1(lx),vaux2(lx)
      INTEGER                          :: k, i, ik

!    1D real to complex in x

      do k=1,lz
        do i=1,ly
          call rfftwnd_f77_one_real_to_complex_(plan_RC,v(:,i,k),0)
        enddo
      enddo

      v=v/float(lx)

      call mpitransposeXY(v,id)

!    1D fft complex-to-complex in y
      do k=1,lz
        do i=1,ly+lyext,2
         do ik=1,lx
         vaux1(ik)=cmplx(v(ik,i,k),v(ik,i+1,k))
         enddo
         call fftw_f77_one_(plan_F,vaux1,vaux2)
         do ik=1,lx
         v(ik,i,k)=real(vaux2(ik))
         v(ik,i+1,k)=aimag(vaux2(ik))
         enddo
        enddo
      enddo

      v=v/float(lx)

      call mpitransposeYZ(v,id)

!    1D fft complex-to-complex in z
      do k=1,lz
        do i=1,ly+lyext,2
         do ik=1,lx
         vaux1(ik)=cmplx(v(ik,i,k),v(ik,i+1,k))
         enddo
         call fftw_f77_one_(plan_F,vaux1,vaux2)
         do ik=1,lx
         v(ik,i,k)=real(vaux2(ik))
         v(ik,i+1,k)=aimag(vaux2(ik))
         enddo
        enddo
      enddo

      v=v/float(lx)

      RETURN
      END SUBROUTINE  mpifft3DRC
!========================================================================

! 3D Complex to Real FFT

      SUBROUTINE  mpifft3DCR (v,id)

      USE params

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'plans.h'

      INTEGER, Intent (in)             :: id
      REAL,    Intent (inout)          :: v(lx+2,ly+lyext,lz)
      COMPLEX                          :: vaux1(lx),vaux2(lx)
      INTEGER                          :: k, i, ik

      do k=1,lz
        do i=1,ly+lyext,2
         do ik=1,lx
         vaux1(ik)=cmplx(v(ik,i,k),v(ik,i+1,k))
         enddo
         call fftw_f77_one_(plan_B,vaux1,vaux2)
         do ik=1,lx
         v(ik,i,k)=real(vaux2(ik))
         v(ik,i+1,k)=aimag(vaux2(ik))
         enddo
        enddo
      enddo

      call mpitransposeYZ(v,id)

      do k=1,lz
        do i=1,ly+lyext,2
         do ik=1,lx
         vaux1(ik)=cmplx(v(ik,i,k),v(ik,i+1,k))
         enddo
         call fftw_f77_one_(plan_B,vaux1,vaux2)
         do ik=1,lx
         v(ik,i,k)=real(vaux2(ik))
         v(ik,i+1,k)=aimag(vaux2(ik))
         enddo
        enddo
      enddo

      call mpitransposeXY(v,id)

      do k=1,lz
        do i=1,ly
          call rfftwnd_f77_one_complex_to_real_(plan_CR,v(:,i,k),0)
        enddo
      enddo

       RETURN
       END SUBROUTINE  mpifft3DCR
!======================================================================

      SUBROUTINE mpitransposeXY(v,id)

      USE params 
      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, Intent (in)          :: id
      REAL, Intent (inout)          :: v(lx+2,ly+lyext,lz)
      REAL                        :: tmp(ly+lyext,lx,lz)
      REAL                        :: tmp1(lx/nprocY+2,ly,lz)
      REAL                        :: tmp2(lx/nprocY+2,ly,lz)
      INTEGER            ::  isize,nzpsend,nzprecv
      INTEGER            ::  status(MPI_STATUS_SIZE,2),req(2)
      INTEGER            ::  lllx,i,js,j,j1,ks,k,k1,ierr,llxh,ii

        lllx=lx/nprocY
        isize=(lllx+2)*ly*lz

        do i=1,nprocY-1
          nzpsend=mod(idlocalY+i,nprocY)
          nzprecv=mod(idlocalY+nprocY-i,nprocY)
          nzpsend=nzpsend+int(id/nprocY)*nprocY
          nzprecv=nzprecv+int(id/nprocY)*nprocY
          js=(nzpsend-int(id/nprocY)*nprocY)*lllx

          llxh=lllx
          if ((mod(nzpsend,nprocY)+1).eq.nprocY) llxh=lllx+2
          do ii=1,lz
            do k=1,ly
              do j=1,llxh
                j1=js+j
                tmp1(j,k,ii)=v(j1,k,ii)
              enddo
            enddo
          enddo

          call mpi_isend(tmp1,isize,MPI_REAL8,nzpsend,i,MPI_COMM_WORLD, &
                         req(1), ierr)
          call mpi_irecv(tmp2,isize,MPI_REAL8,nzprecv,i,MPI_COMM_WORLD, &
                         req(2),ierr)
          call mpi_waitall(2,req,status,ierr)

          ks=(nzprecv-int(id/nprocY)*nprocY)*ly

          llxh=lllx
          if ((mod(id,nprocY)+1).eq.nprocY) llxh=lllx+2
          do ii=1,lz
            do k=1,ly
              k1=ks+k
              do j=1,llxh
                tmp(j,k1,ii)=tmp2(j,k,ii)
              end do
            enddo
          enddo
        end do

        ks=mod(id,nprocY)*ly
        js=mod(id,nprocY)*lllx
        llxh=lllx
        if ((mod(id,nprocY)+1).eq.nprocY) llxh=lllx+2
        do i=1,lz
          do k=1,ly
            k1=ks+k
            do j=1,llxh
              j1=js+j
              tmp(j,k1,i)=v(j1,k,i)
            end do
          end do
        end do

        do i=1,lz
          do j=1,lx
            do k=1,ly+lyext
               v(j,k,i)=tmp(k,j,i)
             end do
          end do
        end do

        RETURN
        END SUBROUTINE mpitransposeXY
!===================================================================

      SUBROUTINE mpitransposeYZ(v,id)

      USE params
      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, Intent (in)          :: id
      REAL, Intent (inout)          :: v(lx+2,ly+lyext,lz)  
      REAL                        :: tmp(lx/nprocZ,ly+lyext,lz*nprocZ)
      REAL                        :: tmp1(lx/nprocZ,ly+lyext,lz)
      REAL                        :: tmp2(lx/nprocZ,ly+lyext,lz)
      INTEGER            ::  isize,nzpsend,nzprecv
      INTEGER            ::  status(MPI_STATUS_SIZE,2),req(2)
      INTEGER            ::  lllx,i,js,j,j1,ks,k,k1,ierr,ii

        lllx=lx/nprocZ
        isize=lllx*(ly+lyext)*lz

        do i=1,nprocZ-1
          nzpsend=mod(idlocalZ+i,nprocZ)
          nzprecv=mod(idlocalZ+nprocZ-i,nprocZ)
          nzpsend=id+(nzpsend-int(id/nprocY))*nprocY
          nzprecv=id+(nzprecv-int(id/nprocY))*nprocY

          js=((nzpsend-id)/nprocY+INT(id/nprocY))*lllx
          do ii=1,lz
            do k=1,ly+lyext
              do j=1,lllx
                j1=js+j
                tmp1(j,k,ii)=v(j1,k,ii)
              enddo
            enddo
          enddo

          call mpi_isend(tmp1,isize,MPI_REAL8,nzpsend,i,MPI_COMM_WORLD, &
                         req(1), ierr)
          call mpi_irecv(tmp2,isize,MPI_REAL8,nzprecv,i,MPI_COMM_WORLD, & 
                         req(2),ierr)
          call mpi_waitall(2,req,status,ierr)

          ks=((nzprecv-id)/nprocY+INT(id/nprocY))*lz
          do k=1,lz
            k1=ks+k
            do j=1,ly+lyext
              do ii=1,lllx
                tmp(ii,j,k1)=tmp2(ii,j,k)
              end do
            enddo
          enddo
        end do

        ks=INT(id/nprocY)*lz
        js=INT(id/nprocY)*lllx
        do k=1,lz
          k1=ks+k
          do i=1,ly+lyext
            do j=1,lllx
              j1=js+j
              tmp(j,i,k1)=v(j1,i,k)
            end do
          end do
        end do

        do j=1,lx
          do i=1,ly+lyext
            do k=1,lz
               v(j,i,k)=tmp(k,i,j)
            end do
          end do
        end do

        RETURN
        END SUBROUTINE mpitransposeYZ
!===================================================================

