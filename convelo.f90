!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE SAVEOMP (wx,wy,wz,ox,oy,oz,destination,id)
        
!     This subroutine saves velocity in one file for OpenMP runs.

      USE params 
      INCLUDE 'mpif.h'

      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   :: wxyzr, wxyzrt
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   :: vorticity
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   :: v_tmp1

      REAL,    DIMENSION (nx+2,ly,lz)           :: wx, wy, wz
      REAL,    DIMENSION (nx+2,ly,lz)           :: ox, oy, oz

      INTEGER, DIMENSION (MPI_STATUS_SIZE)      :: status
      INTEGER                                   :: idy, idz, ls, ierr, error 

      CHARACTER*50  destination
      CHARACTER*70  fvort

      ALLOCATE ( vorticity (nx+2,ny,nz) )
      ALLOCATE ( wxyzr (nx,ly,lz)       )

      vorticity (nx+1,:,:) = 0.
      vorticity (nx+2,:,:) = 0.

      ls = nx * ly * lz

      if (id.eq.0) write(*,*) '==== converting data to OpenMP form ===='

!     for vx

      wxyzr (1:nx,1:ly,1:lz) = wx (1:nx,1:ly,1:lz) 

      IF (id.eq.0) THEN

       vorticity (1:nx,1:ly,1:lz) = wxyzr

       do n=1, nproc-1

        idy = mod (n, nprocY)
        idz = int ( n / nprocY )

        ALLOCATE ( wxyzrt(nx,ly,lz) )

        CALL MPI_RECV(wxyzrt,ls,MPI_REAL8,n,111,MPI_COMM_WORLD,status,ierr)

        vorticity (1:nx,(1+idy*ly):(ly+idy*ly),(1+idz*lz):(lz+idz*lz)) = wxyzrt(1:nx,1:ly,1:lz)

        DEALLOCATE ( wxyzrt )

       enddo

      ELSE
       CALL MPI_SEND (wxyzr,ls,MPI_REAL8,0,111,MPI_COMM_WORLD,status,ierr)
      ENDIF

      CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

      IF (id.eq.0) THEN

        ALLOCATE ( v_tmp1 (nx,ny,nz) )
                                                                                                                                                              
        do i = 1, nx
         do k = 1, nz
          v_tmp1 (k,:,i) = vorticity (i,:,k)
         enddo
        enddo
                                                                                                                                                              
        do i = 1, nx
         do j = 1, ny
          vorticity (j,i,:) = v_tmp1 (i,j,:)
         enddo
        enddo
                                                                                                                                                              
        DEALLOCATE ( v_tmp1 )

       write(*,*) id, MAXVAL(real(vorticity(:,:,:)))
       write(*,*) 'vx'
       write(*,*) '========================================='

!       do ijk=1,64
!        write(*,*) vorticity (64,21,ijk)
!       enddo

       fvort = TRIM(destination)//'vx'

       open(10,file=fvort,status='unknown',form='unformatted')
       write(10) vorticity
       close(10)

      ENDIF

      CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

!      for vy

       wxyzr (1:nx,1:ly,1:lz) = wy (1:nx,1:ly,1:lz) 

       IF (id.eq.0) THEN

        vorticity (1:nx,1:ly,1:lz) = wxyzr

        do n=1, nproc-1

         idy = mod (n, nprocY)
         idz = int ( n / nprocY )

         ALLOCATE ( wxyzrt(nx,ly,lz) )

         CALL MPI_RECV(wxyzrt,ls,MPI_REAL8,n,111,MPI_COMM_WORLD,status,ierr)

         vorticity (1:nx,(1+idy*ly):(ly+idy*ly),(1+idz*lz):(lz+idz*lz)) = wxyzrt(1:nx,1:ly,1:lz)

         DEALLOCATE ( wxyzrt )

        enddo

       ELSE
        CALL MPI_SEND (wxyzr,ls,MPI_REAL8,0,111,MPI_COMM_WORLD,status,ierr)
       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

       IF (id.eq.0) THEN
                                                                                                                                                              
        ALLOCATE ( v_tmp1 (nx,ny,nz) )
                                                                                                                                                              
        do i = 1, nx
         do k = 1, nz
          v_tmp1 (k,:,i) = vorticity (i,:,k)
         enddo
        enddo
                                                                                                                                                              
        do i = 1, nx
         do j = 1, ny
          vorticity (j,i,:) = v_tmp1 (i,j,:)
         enddo
        enddo
                                                                                                                                                              
        DEALLOCATE ( v_tmp1 )

        write(*,*) id, MAXVAL(real(vorticity(:,:,:)))
        write(*,*) 'vy'
        write(*,*) '========================================='

!        do ijk=1,64
!        write(*,*) vorticity (64,21,ijk)
!        enddo

        fvort = TRIM(destination)//'vy'

        open(10,file=fvort,status='unknown',form='unformatted')
        write(10) vorticity
        close(10)

       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

!      for vz

       wxyzr (1:nx,1:ly,1:lz) = wz (1:nx,1:ly,1:lz) 

       IF (id.eq.0) THEN

        vorticity (1:nx,1:ly,1:lz) = wxyzr

        do n=1, nproc-1
                                                                                                                                                              
         idy = mod (n, nprocY)
         idz = int ( n / nprocY )

         ALLOCATE ( wxyzrt(nx,ly,lz) )

         CALL MPI_RECV(wxyzrt,ls,MPI_REAL8,n,111,MPI_COMM_WORLD,status,ierr)

         vorticity (1:nx,(1+idy*ly):(ly+idy*ly),(1+idz*lz):(lz+idz*lz)) = wxyzrt(1:nx,1:ly,1:lz)

         DEALLOCATE ( wxyzrt )

        enddo

       ELSE
        CALL MPI_SEND (wxyzr,ls,MPI_REAL8,0,111,MPI_COMM_WORLD,status,ierr)
       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

       IF (id.eq.0) THEN
                                                                                                                                                              
        ALLOCATE ( v_tmp1 (nx,ny,nz) )
                                                                                                                                                              
        do i = 1, nx
         do k = 1, nz
          v_tmp1 (k,:,i) = vorticity (i,:,k)
         enddo
        enddo
                                                                                                                                                              
        do i = 1, nx
         do j = 1, ny
          vorticity (j,i,:) = v_tmp1 (i,j,:)
         enddo
        enddo
                                                                                                                                                              
        DEALLOCATE ( v_tmp1 )

        write(*,*) id, MAXVAL(real(vorticity(:,:,:)))
        write(*,*) 'vz'
        write(*,*) '========================================='

!        do ijk=1,64
!        write(*,*) vorticity (64,21,ijk)
!        enddo

        fvort = TRIM(destination)//'vz'

        open(10,file=fvort,status='unknown',form='unformatted')
        write(10) vorticity
        close(10)

       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

!      for ox

       wxyzr (1:nx,1:ly,1:lz) = ox (1:nx,1:ly,1:lz) 

       IF (id.eq.0) THEN

        vorticity (1:nx,1:ly,1:lz) = wxyzr

        do n=1, nproc-1
                                                                                                                                                              
         idy = mod (n, nprocY)
         idz = int ( n / nprocY )
                                                                                                                                                              
         ALLOCATE ( wxyzrt(nx,ly,lz) )

         CALL MPI_RECV(wxyzrt,ls,MPI_REAL8,n,111,MPI_COMM_WORLD,status,ierr)

         vorticity (1:nx,(1+idy*ly):(ly+idy*ly),(1+idz*lz):(lz+idz*lz)) = wxyzrt(1:nx,1:ly,1:lz)

         DEALLOCATE ( wxyzrt )

        enddo

       ELSE
        CALL MPI_SEND (wxyzr,ls,MPI_REAL8,0,111,MPI_COMM_WORLD,status,ierr)
       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

       IF (id.eq.0) THEN
                                                                                                                                                              
        ALLOCATE ( v_tmp1 (nx,ny,nz) )
                                                                                                                                                              
        do i = 1, nx
         do k = 1, nz
          v_tmp1 (k,:,i) = vorticity (i,:,k)
         enddo
        enddo
                                                                                                                                                              
        do i = 1, nx
         do j = 1, ny
          vorticity (j,i,:) = v_tmp1 (i,j,:)
         enddo
        enddo
                                                                                                                                                              
        DEALLOCATE ( v_tmp1 )

        write(*,*) id, MAXVAL(real(vorticity(:,:,:)))
        write(*,*) 'ox'
        write(*,*) '========================================='

!        do ijk=1,64
!        write(*,*) vorticity (64,21,ijk)
!        enddo

        fvort = TRIM(destination)//'ox'

        open(10,file=fvort,status='unknown',form='unformatted')
        write(10) vorticity
        close(10)

       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

!      for oy

       wxyzr (1:nx,1:ly,1:lz) = oy (1:nx,1:ly,1:lz) 

       IF (id.eq.0) THEN

        vorticity (1:nx,1:ly,1:lz) = wxyzr

        do n=1, nproc-1
                                                                                                                                                              
         idy = mod (n, nprocY)
         idz = int ( n / nprocY )

         ALLOCATE ( wxyzrt(nx,ly,lz) )

         CALL MPI_RECV(wxyzrt,ls,MPI_REAL8,n,111,MPI_COMM_WORLD,status,ierr)

         vorticity (1:nx,(1+idy*ly):(ly+idy*ly),(1+idz*lz):(lz+idz*lz)) = wxyzrt(1:nx,1:ly,1:lz)

         DEALLOCATE ( wxyzrt )

        enddo

       ELSE
        CALL MPI_SEND (wxyzr,ls,MPI_REAL8,0,111,MPI_COMM_WORLD,status,ierr)
       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

       IF (id.eq.0) THEN

        ALLOCATE ( v_tmp1 (nx,ny,nz) )
                                                                                                                                                              
        do i = 1, nx
         do k = 1, nz
          v_tmp1 (k,:,i) = vorticity (i,:,k)
         enddo
        enddo
                                                                                                                                                              
        do i = 1, nx
         do j = 1, ny
          vorticity (j,i,:) = v_tmp1 (i,j,:)
         enddo
        enddo
                                                                                                                                                              
        DEALLOCATE ( v_tmp1 )
                                                                                                                                                              
        write(*,*) id, MAXVAL(real(vorticity(:,:,:)))
        write(*,*) 'oy'
        write(*,*) '========================================='

!        do ijk=1,64
!        write(*,*) vorticity (64,21,ijk)
!        enddo

        fvort = TRIM(destination)//'oy'
        open(10,file=fvort,status='unknown',form='unformatted')
        write(10) vorticity
        close(10)

       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

!      for oz

       wxyzr (1:nx,1:ly,1:lz) = oz (1:nx,1:ly,1:lz) 

       IF (id.eq.0) THEN

        vorticity (1:nx,1:ly,1:lz) = wxyzr

        do n=1, nproc-1
                                                                                                                                                              
         idy = mod (n, nprocY)
         idz = int ( n / nprocY )

         ALLOCATE ( wxyzrt(nx,ly,lz) )

         CALL MPI_RECV (wxyzrt,ls,MPI_REAL8,n,111,MPI_COMM_WORLD,status,ierr)

         vorticity (1:nx,(1+idy*ly):(ly+idy*ly),(1+idz*lz):(lz+idz*lz)) = wxyzrt(1:nx,1:ly,1:lz)

         DEALLOCATE ( wxyzrt )

        enddo

       ELSE
        CALL MPI_SEND (wxyzr,ls,MPI_REAL8,0,111,MPI_COMM_WORLD,status,ierr)
       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

       IF (id.eq.0) THEN

        ALLOCATE ( v_tmp1 (nx,ny,nz) )

        do i = 1, nx
         do k = 1, nz
          v_tmp1 (k,:,i) = vorticity (i,:,k)
         enddo
        enddo

        do i = 1, nx
         do j = 1, ny
          vorticity (j,i,:) = v_tmp1 (i,j,:)
         enddo
        enddo

        DEALLOCATE ( v_tmp1 )

        write(*,*) id, MAXVAL(real(vorticity(:,:,:)))
        write(*,*) 'oz'
        write(*,*) '========================================='
                                                                                                                                                              
!        do ijk=1,64
!        write(*,*) vorticity (64,21,ijk)
!        enddo

        fvort = TRIM(destination)//'oz'
        open(10,file=fvort,status='unknown',form='unformatted')
        write(10) vorticity
        close(10)

       ENDIF

       CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

       DEALLOCATE ( vorticity )
       DEALLOCATE ( wxyzr     )

      END SUBROUTINE SAVEOMP

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
