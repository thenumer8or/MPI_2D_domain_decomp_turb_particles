!==================================================================
      SUBROUTINE input_output(directvelo,ux,uy,uz,ox,oy,oz,idump,id,idirec)

      USE params
      IMPLICIT NONE

      REAL,DIMENSION   (llx,lly,lz)  ::  ux,uy,uz
      REAL,DIMENSION   (llx,lly,lz)  ::  ox,oy,oz
      INTEGER                        ::  i1d,i2d,i3d,i4d,i4dump,i5dump,id,idump,idirec

      CHARACTER*62 fnm1
      CHARACTER*50 directvelo

      i1d=int(id/1000)
      i2d=int((id-i1d*1000)/100)
      i3d=int((id-i1d*1000-i2d*100)/10)
      i4d=id-i1d*1000-i2d*100-i3d*10

!      i1d=int(id/100)
!      i2d=int((id-i1d*100)/10)
!      i3d=id-i1d*100-i2d*10

      if(idump.lt.10) then
        fnm1=TRIM(directvelo)//'vel2D'//char(idump+48)//'.'                &
              //char(i1d+48)//char(i2d+48)//char(i3d+48)//char(i4d+48)
      else
        i4dump=int(idump/10)
        i5dump=mod(idump,10)
        fnm1=TRIM(directvelo)//'vel2D'//char(i4dump+48)//char(i5dump+48)//'.'     &
              //char(i1d+48)//char(i2d+48)//char(i3d+48)//char(i4d+48)
      endif

      if(idirec.eq.1) then 
        open(10,file=fnm1,status='unknown',form='unformatted')
        read(10)ux
        read(10)uy
        read(10)uz
        read(10)ox
        read(10)oy
        read(10)oz
        close(10)
      elseif(idirec.eq.2) then
        open(10,file=fnm1,status='unknown',form='unformatted')
        write(10)ux
        write(10)uy
        write(10)uz
        write(10)ox
        write(10)oy
        write(10)oz
        close(10)
      endif
 

      RETURN
      END SUBROUTINE input_output

!====================================================================
      SUBROUTINE input_outputf (directfor2,b1r,b2r,b3r,iseedff,iv,iy,idump,id,idirec)

      USE params
      IMPLICIT NONE

      REAL                         ::  b1r(6,5,5),b2r(6,5,5),b3r(6,5,5)
      INTEGER                      ::  iseedff,iy,iv(32)
      INTEGER                      ::  idump,id,idirec,i4dump,i5dump
      CHARACTER*62 fnm1
      CHARACTER*50 directfor2

      if (id.ne.0) RETURN

      if(idump.lt.10) then
        fnm1=TRIM(directfor2)//'force'//'.'//char(idump+48)
      else
        i4dump=int(idump/10)
        i5dump=mod(idump,10)
        fnm1=TRIM(directfor2)//'force'//'.'//char(i4dump+48)//char(i5dump+48)
      endif

      if(idirec.eq.1) then
        open(10,file=fnm1,status='unknown',form='unformatted')
        read(10)iseedff
        read(10)iv
        read(10)iy
        read(10)b1r,b2r,b3r
        close(10)
      elseif(idirec.eq.2) then
        open(10,file=fnm1,status='unknown',form='unformatted')
        write(10)iseedff
        write(10)iv
        write(10)iy
        write(10)b1r,b2r,b3r
        close(10)
      endif

      RETURN
      END SUBROUTINE input_outputf

!====================================================================

      SUBROUTINE LOAD_PARTICLE (id,destination)
                                                                                                                                                                  
      USE params
!      IMPLICIT NONE
                                                                                                                                                                  
      CHARACTER*50 destination
      CHARACTER*70 fnm1

      i1d=int(id/1000)
      i2d=int((id-i1d*1000)/100)
      i3d=int((id-i1d*1000-i2d*100)/10)
      i4d=id-i1d*1000-i2d*100-i3d*10

!      i1d = int(id/100)
!      i2d = int((id-i1d*100)/10)
!      i3d = id-i1d*100-i2d*10

      fnm1 = TRIM(destination)//'part.'//char(i1d+48)//char(i2d+48)//char(i3d+48)//char(i4d+48)
                                                                                                                                                                  
      open(12,file=fnm1,status='unknown',form='unformatted')
      read(12) nps
      read(12) yp(:,1:nps), fvp(:,0:2,1:nps), fdp(:,0:3,1:nps)
      close(12)

      RETURN
      END SUBROUTINE LOAD_PARTICLE
                                                                                                                                                                  
!====================================================================

      SUBROUTINE SAVE_PARTICLE (id,destination)
                                                                                                                                                                  
      USE params
!      IMPLICIT NONE
                                                                                                                                                                  
      CHARACTER*50 destination
      CHARACTER*70 fnm1

      i1d=int(id/1000)
      i2d=int((id-i1d*1000)/100)
      i3d=int((id-i1d*1000-i2d*100)/10)
      i4d=id-i1d*1000-i2d*100-i3d*10

!      i1d = int(id/100)
!      i2d = int((id-i1d*100)/10)
!      i3d = id-i1d*100-i2d*10

      fnm1 = TRIM(destination)//'part.'//char(i1d+48)//char(i2d+48)//char(i3d+48)//char(i4d+48)
                                                                                                                                                                  
      open(12,file=fnm1,status='unknown',form='unformatted')
      write(12) nps
      write(12) yp(:,1:nps), fvp(:,0:2,1:nps), fdp(:,0:3,1:nps)
      close(12)

      RETURN
      END SUBROUTINE SAVE_PARTICLE

!====================================================================

      SUBROUTINE SAVE_PART_GLOBAL (id,destination)
                                                                                                                                                             
      USE params
      INCLUDE 'mpif.h'
                                                                                                                                                             
      CHARACTER*50 destination
      CHARACTER*70 fnm1, fnm2
                                                                                                                                                             
      REAL, ALLOCATABLE, DIMENSION    (:,:)     :: data, wxyzrt 
      INTEGER, DIMENSION (MPI_STATUS_SIZE)      :: status
      INTEGER                                   :: ierr, error, npsgather(0:nproc-1) 

      ALLOCATE ( data (1:8, npart) )

      data (:,:) = 0.

      CALL MPI_GATHER (nps,1,MPI_INTEGER,npsgather,1,MPI_INTEGER,0,mpi_comm_world,ierr)

      if (id.eq.0) write(*,*) '==== gathering particle data to a single file ===='
!      if (id.eq.0) write(*,*) 'npsgather = ', npsgather

      ipart1 = 0 
      ipart2 = 0 

      IF (id.eq.0) THEN

       do i = 1, nps

        if ( yp(7,i) .eq. radius(1) ) then
         ipart1 = ipart1 + 1
         data (1:8,ipart1) = yp(1:8,i)
        elseif ( yp(7,i) .eq. radius(2) ) then
         ipart2 = ipart2 + 1
         data (1:8,npart/2+ipart2) = yp(1:8,i)
        endif

       enddo

       do n=1, nproc-1

        ALLOCATE ( wxyzrt(8,npsgather(n)) )

        CALL MPI_RECV(wxyzrt,8*npsgather(n),MPI_REAL8,n,111,MPI_COMM_WORLD,status,ierr)

        do i = 1, npsgather(n)
                                                                                                                                                             
         if ( wxyzrt(7,i) .eq. radius(1) ) then
          ipart1 = ipart1 + 1
          data (1:8,ipart1) = wxyzrt(1:8,i)
         elseif ( wxyzrt(7,i) .eq. radius(2) ) then
          ipart2 = ipart2 + 1
          data (1:8,npart/2+ipart2) = wxyzrt(1:8,i)
         endif
                                                                                                                                                             
        enddo

         DEALLOCATE ( wxyzrt )

       enddo

      ELSE
       CALL MPI_SEND (yp(1:8,1:nps),8*nps,MPI_REAL8,0,111,MPI_COMM_WORLD,status,ierr)
      ENDIF

      CALL MPI_BARRIER ( MPI_COMM_WORLD, error )

      IF (id.eq.0) THEN

       fnm1 = TRIM(destination)//'part_loc'
       open(12,file=fnm1,status='unknown',form='unformatted')
       write(12) data (1:3,1:npart)
       close(12)

       write(999,6158) data (1:8,1:npart)
6158   format ( 8(1pe25.15) )

       fnm2 = TRIM(destination)//'part_vel'
       open(13,file=fnm2,status='unknown',form='unformatted')
       write(13) data (4:6,1:npart)
       close(13)
                                                                                                                                                             
      ENDIF

      DEALLOCATE ( data )

      RETURN
      END SUBROUTINE SAVE_PART_GLOBAL
                                                                                                                                                             
!====================================================================
