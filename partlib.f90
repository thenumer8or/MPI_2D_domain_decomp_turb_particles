!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

! This file contains the subroutines required to track particles in turbulence in Lagrangian way in a 2D-decomposed domain.
! Coded by: hparish.
! Implementation is reviewed and verified by: oayala and brosa.
! Questions: Please email h.parish@uci.edu
! Final version committed to SVN: Aug. 2015.
! The 3D domain-decomposition implementation is pending.

! Main flow solver is turb.f90
! Other accessory files that are needed to run this code are in the current folder. 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE PART_INIT (id, part_in)

!     This subroutine initializes and sets up particle location

        USE params
        INCLUDE 'mpif.h'

        REAL, ALLOCATABLE, DIMENSION  (:,:)    ::  location, location2

        INTEGER                                ::  ierr, error
        INTEGER                                ::  iseedp, iyp, ivp(NTAB)
        INTEGER, DIMENSION (MPI_STATUS_SIZE)   ::  status

        CHARACTER*50                           ::   part_in

        nuniform = int ( real(npart) / real(nproc) )
        msize    = int ( 20 * nuniform )

        ncdx = nx 
        ncdy = nx / nprocY 
        ncdz = nx / nprocZ 
        nbin = 180 

        ncollmax = 100 

        ALLOCATE ( location (3,npart)  )
        ALLOCATE ( location2 (3,npart) )

        ALLOCATE ( lhnode (msize,3)    )
        ALLOCATE ( pindex (msize,8)    )
        ALLOCATE ( yp (8,msize)        )
        ALLOCATE ( vp (3,msize)        )
        ALLOCATE ( bg (msize,6,3)      )
        ALLOCATE ( enst_p(msize)       )
        ALLOCATE ( disp_p(msize)       )

        ALLOCATE ( dmove(3,msize)      )
        ALLOCATE ( fvp(3,0:2,msize)    )
        ALLOCATE ( fdp(3,0:3,msize)    )

        ALLOCATE ( removepart(msize,0:4)        )
        ALLOCATE ( ncollpart(0:4)               )
        ALLOCATE ( nirbin(nset,nset,0:nbin-1)   )
        ALLOCATE ( wrbin(nset,nset,0:nbin-1)    )
        ALLOCATE ( wxbin(nset,nset,0:nbin-1)    )
        ALLOCATE ( wybin(nset,nset,0:nbin-1)    )
        ALLOCATE ( wzbin(nset,nset,0:nbin-1)    )
        ALLOCATE ( urbin(nset,nset,0:nbin-1)    )
        ALLOCATE ( thetabin(nset,nset,0:nbin-1) )
        ALLOCATE ( wrtanbin(nset,nset,0:nbin-1) )
        ALLOCATE ( collist (20, ncollmax)       )

        IF ( HDI_included ) THEN

         ALLOCATE ( pertvel(3,msize)      )
         ALLOCATE ( Axreshaped(3*msize,1) )
         ALLOCATE ( breshaped(3*msize,1)  )
!         ALLOCATE ( nei_num(msize)       )
!         ALLOCATE ( nei_data(7,40,msize) )

         pertvel = 0.0

        ENDIF

        nirbin   = 0
        wrbin    = 0.0
        wxbin    = 0.0
        wybin    = 0.0
        wzbin    = 0.0
        urbin    = 0.0
        thetabin = 0.0
        wrtanbin = 0.0

!       parameters for random generator for particle relocation
        iseedr = -id * 10
        idum2  = 0
        ivr    = 0
        iyr    = 0

!       if (id.eq.0) write(*,*) 'New particles - draw location'
!       iseedp = 290811 changed this value to initiate part at different locations.
        iseedp = 290710
        iseedp = -iseedp
        ivp(:) = 0
        iyp    = 0

!-------------------------------------------------------------------
!      scaling variables, input data are in cgs. 

       etk_CLOUD  = sqrt ( sqrt ( viscosity_CLOUD ** 3 / ediss_CLOUD ) )
       tauk_CLOUD = sqrt ( viscosity_CLOUD / ediss_CLOUD )
       vk_CLOUD   = etk_CLOUD / tauk_CLOUD
                                                                                                                                                                     
       vk_DNS     = etk_DNS    / tauk_DNS

       vk_ratio   = vk_CLOUD   / vk_DNS
       etk_ratio  = etk_CLOUD  / etk_DNS
       tauk_ratio = tauk_CLOUD / tauk_DNS

       radius(1:nset) = radius(1:nset) * ( etk_DNS / etk_CLOUD ) / 10000.0
       radmax         = MAXVAL ( radius(1:nset) )
                                                                                                                                                                     
       tau_CLOUD_const  = 2. * rho_water * etk_ratio**2 / ( 9. * rho_air * viscosity_CLOUD )
       Sv_const         = tau_CLOUD_const * gravity / vk_ratio
       gravity_dns      = gravity * tauk_ratio / vk_ratio
       tau_CLOUD_const2 = tau_CLOUD_const / tauk_ratio
                                                                                                                                                                     
       if ( nset .eq. 1 ) then
        ucharac = Sv_const * radius(1) ** 2
       elseif ( nset .gt. 1 ) then
        ucharac = Sv_const * abs ( radius(1) ** 2 - radius(nset) ** 2 )
       endif

!-------------------------------------------------------------------

        d_coll  = 20.0 * radmax 
        nx_coll = 1 + int ( d_coll * real(ncdx) /  hx ) 
        ny_coll = 1 + int ( d_coll * real(ncdy) /  hy ) 
        nz_coll = 1 + int ( d_coll * real(ncdz) /  hz ) 

        nx_max = nx_coll 
        ny_max = ny_coll 
        nz_max = nz_coll

        if ( HDI_included ) then

         d_hdi  = HDI_trunc * radmax 
         nx_hdi = 1 + int ( d_hdi * real(ncdx) /  hx ) 
         ny_hdi = 1 + int ( d_hdi * real(ncdy) /  hy ) 
         nz_hdi = 1 + int ( d_hdi * real(ncdz) /  hz ) 

         nx_max = MAX ( nx_coll, nx_hdi )
         ny_max = MAX ( ny_coll, ny_hdi )
         nz_max = MAX ( nz_coll, nz_hdi )

        endif

        if ( id .eq. 0 ) write (*,*)'ncdx, ncdy, ncdz, d_coll, d_hdi, nx_max, ny_max, nz_max = ', ncdx, ncdy, ncdz, d_coll, d_hdi, nx_max, ny_max, nz_max 

        area1   = hy * hz
        area2   = ( hy + d_hdi ) * ( hz + d_hdi )
        aratio  = area2 / area1
        newsize = int ( 30.0 * aratio * nuniform )
 
        ALLOCATE ( head ( 1:ncdx, (1-ny_max):(ncdy+ny_max), (1-nz_max):(ncdz+nz_max) ) )
        ALLOCATE ( list ( newsize ) )
        ALLOCATE ( vpn (3,newsize)  )

        if ( HDI_included ) then
         ALLOCATE ( ypn (21,newsize) )
        else
         ALLOCATE ( ypn (15,newsize) )
        endif

!-----  COLLISION: check if slab thickness > maximum radius of gathering statistics.

        if ( id .EQ. 0 ) then
          if ( ( hz .GT. d_coll ) .AND. ( hy .GT. d_coll ) ) then
            write(*,*) 'd_coll, hy, hz = ', d_coll, hy, hz
            write(*,*) 'code is able to handle collision statistics.'
          else
            write(*,*) 'd_coll, hy, hz = ', d_coll, hy, hz
            write(*,*) 'code is not able to handle collision statistics. increase cell width.'
          endif
        endif
 
        if ( ( hz .LE. d_coll ) .OR. ( hy .LE. d_coll ) ) then
         CALL mpi_finalize(ierr)
         STOP
        endif

!-----  COLLISION: check if cell size > maximum collision radius.
                                                                                                                                                                      
        if ( id .EQ. 0 ) then
          if ( ( hz .GT. 2.0*radmax ) .AND. ( hy .GT. 2.0*radmax ) ) then
            write(*,*) 'Rcollmax, hy, hz = ', 2.0*radmax, hy, hz
            write(*,*) 'code is able to handle collision detection'
          else
            write(*,*) 'Rcollmax, hy, hz = ', 2.0*radmax, hy, hz
            write(*,*) 'code is not able to handle collision detection. increase cell width.'
          endif
        endif
                                                                                                                                                                      
        if ( ( hz .LE. 2.0*radmax ) .OR. ( hy .LE. 2.0*radmax ) ) then
         CALL mpi_finalize(ierr)
         STOP
        endif

!-----  HDI: check if slab thickness > maximum truncation radius.
                                                                                                                                                                      
        IF ( HDI_included ) THEN
                                                                                                                                                                      
        if ( id .EQ. 0 ) then
          if ( ( hz .GT. d_hdi ) .AND. ( hy .GT. d_hdi ) ) then
            write(*,*) 'd_hdi, hy, hz = ', d_hdi, hy, hz
            write(*,*) 'code is able to handle HDI'
          else
            write(*,*) 'd_hdi, hy, hz = ', d_hdi, hy, hz
            write(*,*) 'code is not able to handle HDI. increase subdomain size.'
          endif
        endif
                                                                                                                                                                      
        if ( ( hz .LE. d_hdi ) .OR. ( hy .LE. d_hdi ) ) then
         CALL mpi_finalize(ierr)
         STOP
        endif
                                                                                                                                                                      
       ENDIF

!-------------------------------------------------------------------
!      Set up particle location
                                                                                                                                                                 
       IF ( new_particle ) THEN

        if ( id.eq.0 ) then

         do ih = 1, npart
          location(1,ih) = ranpp (iseedp,iyp,ivp)
          location(2,ih) = ranpp (iseedp,iyp,ivp)
          location(3,ih) = ranpp (iseedp,iyp,ivp)
         enddo

        endif

        if ( id.eq.0 ) then
!============================================================================
!       location2(1:3,:) = pi2 * location(1:3,:) 

!       open(65,file='/ptmp/hparish/2d/part.dat',status='unknown',form='unformatted')
!       write(65) location2
               
!       open(66,file='/ptmp/hparish/2d/part_location.dat',status='unknown')
!       write(66,6059) location2
!6059   format ( 2x,3(1pe25.15) )
               
!       close(66)
!       close(65)
!============================================================================
        endif

!       CALL mpi_finalize(ierr)
!       STOP

        CALL MPI_BCAST (location, 3*npart, MPI_REAL8, 0, mpi_comm_world, ierr)

        yp(:,:) = 0.
 
        do ih = 1, npart

         idy    = FLOOR ( location(2,ih) * real(nprocY) ) + 1
         idz    = FLOOR ( location(3,ih) * real(nprocZ) ) + 1
         idpart = ( idz - 1 ) * nprocY + idy - 1

         if ( idpart .eq. id ) then
 
          nps       = nps + 1
          yp(1,nps) = pi2 * location(1,ih)
          yp(2,nps) = pi2 * location(2,ih)
          yp(3,nps) = pi2 * location(3,ih)
 
          i1        = FLOOR ( (ih-1) * real(nset) / real(npart) ) + 1
          yp(7,nps) = REAL ( radius(i1) )
          yp(8,nps) = REAL ( ih )
 
         endif

        enddo

       ELSE
                                                                                                                                                                  
        CALL LOAD_PARTICLE (id,part_in)
        if (id .eq. 0) write(*,*) 'particle location is read from previous MPI_2d run' 

!        index_yp = 463 + id
!        write (index_yp,6058) yp
!6058    format ( 2x,7(1pe25.15) )
                                                                                                                                                                  
       ENDIF

       CALL NEIGHBORS (id,id_e,id_ne,id_n,id_nw,id_w,id_sw,id_s,id_se)
       nei_id = (/ id_e, id_ne, id_n, id_nw, id_w, id_sw, id_s, id_se /)

       DEALLOCATE ( location  )
       DEALLOCATE ( location2 )
 
       RETURN
       END SUBROUTINE PART_INIT 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE PREPARE_INTERP (id)
 
!     This subroutine sets up interpolation factors bg for the particles and
!     locate grid points. lhnode is an index for the nearest
!     grid point coordinate to the left of the particle.
!     value of lhnode is between 0 and n-1.

      USE params
      INCLUDE 'mpif.h'

      stime = mpi_wtime()
 
      dxyz(1) = pi2 / real(nx)
      dxyz(2) = pi2 / real(ny)
      dxyz(3) = pi2 / real(nz)
                                                                                                                                                                       
       do ij = 1, 3
        do ih = 1, nps
                                                                                                                                                                       
          node          = FLOOR ( yp(ij,ih) / dxyz(ij) )
          lhnode(ih,ij) = node
                                                                                                                                                                       
          z  = ( yp(ij,ih) - REAL(node) * dxyz(ij) ) / dxyz(ij)
          z2 = z * z
          z3 = z2 * z
          z4 = z3 * z
          z5 = z4 * z
                                                                                                                                                                       
          bg(ih,1,ij) = ( 6.0 * z - 5.0 * z2 - 5.0 * z3 + 5.0 * z4 - z5 ) / 120.0
          bg(ih,2,ij) = ( -12.0 * z + 16.0 * z2 - z3 - 4.0 * z4 + z5 ) / 24.0
          bg(ih,3,ij) = ( 12.0 - 4.0 * z - 15.0 * z2 + 5.0 * z3 + 3.0 * z4 - z5 ) / 12.0
          bg(ih,4,ij) = ( 12.0 * z + 8.0 * z2 - 7.0 * z3 - 2.0 * z4 + z5 ) / 12.0
          bg(ih,5,ij) = ( -6.0 * z - z2 + 7.0 * z3 + z4 - z5 ) / 24.0
          bg(ih,6,ij) = ( 4.0 * z - 5.0 * z3 + z5 ) / 120.0
                                                                                                                                                                       
        end do
       end do

      etime          = mpi_wtime()
      time_interp_cp = time_interp_cp + etime - stime

      RETURN
      END SUBROUTINE PREPARE_INTERP

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE INTERPOLATE (v, icpt, id)

!     This subroutine does the 6-pt Lagrange interpolation of...
!     the data field u using interpolation factors bg.
                                                                                                                                                                       
      USE params
      INCLUDE 'mpif.h'
                                                                                                                                                                       
      INTEGER                              :: icpt, id
      INTEGER, DIMENSION (6)               :: ix, iy, iz
      INTEGER                              :: status1(MPI_STATUS_SIZE), error, ierr

      REAL, DIMENSION  (llx,lly,lz)        :: v
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: vext
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: v1s, v1r, v2s, v2r, v3s, v3r, v4s, v4r
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: v5s, v5r, v6s, v6r, v7s, v7r, v8s, v8r 

      stime = mpi_wtime()

      ALLOCATE ( v1s  (lx,3,lz)  )
      ALLOCATE ( v1r  (lx,3,lz)  )
      ALLOCATE ( v2s  (lx,3,3)   )
      ALLOCATE ( v2r  (lx,3,3)   )
      ALLOCATE ( v3s  (lx,ly,3)  )
      ALLOCATE ( v3r  (lx,ly,3)  )
      ALLOCATE ( v4s  (lx,2,3)   )
      ALLOCATE ( v4r  (lx,2,3)   )
      ALLOCATE ( v5s  (lx,2,lz)  )
      ALLOCATE ( v5r  (lx,2,lz)  )
      ALLOCATE ( v6s  (lx,2,2)   )
      ALLOCATE ( v6r  (lx,2,2)   )
      ALLOCATE ( v7s  (lx,ly,2)  )
      ALLOCATE ( v7r  (lx,ly,2)  )
      ALLOCATE ( v8s  (lx,3,2)   )
      ALLOCATE ( v8r  (lx,3,2)   )
      ALLOCATE ( vext (lx,ly+5,lz+5) )

!-------------------------------------------------------------------

      etime          = mpi_wtime()
      time_interp_cp = time_interp_cp + etime - stime
      stime          = mpi_wtime()

      ls1 = 3 * lz * lx
      ls2 = 3 * 3  * lx
      ls3 = 3 * ly * lx
      ls4 = 2 * 3  * lx
      ls5 = 2 * lz * lx
      ls6 = 2 * 2  * lx
      ls7 = 2 * ly * lx
      ls8 = 2 * 3  * lx

      v1s ( 1:lx, 1:3, 1:lz ) = v (1:lx, 1:3, 1:lz)
      v2s ( 1:lx, 1:3, 1:3 )  = v (1:lx, 1:3, 1:3) 
      v3s ( 1:lx, 1:ly, 1:3 ) = v (1:lx, 1:ly, 1:3)
      v4s ( 1:lx, 1:2, 1:3 )  = v (1:lx, (ly-1):ly, 1:3)
      v5s ( 1:lx, 1:2, 1:lz ) = v (1:lx, (ly-1):ly, 1:lz)
      v6s ( 1:lx, 1:2, 1:2 )  = v (1:lx, (ly-1):ly, (lz-1):lz) 
      v7s ( 1:lx, 1:ly, 1:2 ) = v (1:lx, 1:ly, (lz-1):lz)
      v8s ( 1:lx, 1:3, 1:2 )  = v (1:lx, 1:3, (lz-1):lz)

! right
      CALL MPI_SENDRECV (v1s,ls1,MPI_REAL8,id_w,311,v1r,ls1,&
             MPI_REAL8,id_e,311,MPI_COMM_WORLD,status1,ierr1)

! left 
      CALL MPI_SENDRECV (v5s,ls5,MPI_REAL8,id_e,315,v5r,ls5,&
             MPI_REAL8,id_w,315,MPI_COMM_WORLD,status1,ierr5)

! top 
      CALL MPI_SENDRECV (v3s,ls3,MPI_REAL8,id_s,313,v3r,ls3,&
             MPI_REAL8,id_n,313,MPI_COMM_WORLD,status1,ierr3)

! bottom 
      CALL MPI_SENDRECV (v7s,ls7,MPI_REAL8,id_n,317,v7r,ls7,&
             MPI_REAL8,id_s,317,MPI_COMM_WORLD,status1,ierr7)

! se corner
      CALL MPI_SENDRECV (v8s,ls8,MPI_REAL8,id_nw,318,v8r,ls8,&
             MPI_REAL8,id_se,318,MPI_COMM_WORLD,status1,ierr8)

! ne corner
      CALL MPI_SENDRECV (v2s,ls2,MPI_REAL8,id_sw,312,v2r,ls2,&
             MPI_REAL8,id_ne,312,MPI_COMM_WORLD,status1,ierr2)

! nw corner
      CALL MPI_SENDRECV (v4s,ls4,MPI_REAL8,id_se,314,v4r,ls4,&
             MPI_REAL8,id_nw,314,MPI_COMM_WORLD,status1,ierr4)

! sw corner
      CALL MPI_SENDRECV (v6s,ls6,MPI_REAL8,id_ne,316,v6r,ls6,&
             MPI_REAL8,id_sw,316,MPI_COMM_WORLD,status1,ierr6)

      vext ( 1:lx, 3:(ly+2), 3:(lz+2) )           = v   (1:lx, 1:ly, 1:lz)
      vext ( 1:lx, (ly+3):(ly+5), 3:(lz+2) )      = v1r (1:lx, 1:3, 1:lz)
      vext ( 1:lx, (ly+3):(ly+5), (lz+3):(lz+5) ) = v2r (1:lx, 1:3, 1:3)
      vext ( 1:lx, 3:(ly+2), (lz+3):(lz+5) )      = v3r (1:lx, 1:ly, 1:3)
      vext ( 1:lx, 1:2, (lz+3):(lz+5) )           = v4r (1:lx, 1:2, 1:3)
      vext ( 1:lx, 1:2, 3:(lz+2) )                = v5r (1:lx, 1:2, 1:lz)
      vext ( 1:lx, 1:2, 1:2 )                     = v6r (1:lx, 1:2, 1:2)
      vext ( 1:lx, 3:(ly+2), 1:2 )                = v7r (1:lx, 1:ly, 1:2)
      vext ( 1:lx, (ly+3):(ly+5), 1:2 )           = v8r (1:lx, 1:3, 1:2)

!      CALL MPI_BARRIER (MPI_COMM_WORLD, error)

      etime          = mpi_wtime()
      time_interp_cm = time_interp_cm + etime - stime
      stime          = mpi_wtime()

!-------------------------------------------------------------------

      mn1 = lx

      vp(icpt, 1:nps) = 0.0

      DO ih0 = 1, nps

        jx    = lhnode(ih0,1) - 2
        ix(1) = ( mod (jx, mn1) + mn1*(1-isign(1,jx))/2 ) + 1

        if ( ix(1) .le. (mn1-5) ) then
          ix(2) = ix(1) + 1
          ix(3) = ix(2) + 1
          ix(4) = ix(3) + 1
          ix(5) = ix(4) + 1
          ix(6) = ix(5) + 1
        else 
          ix(2) = mod (ix(1), mn1) + 1
          ix(3) = mod (ix(2), mn1) + 1
          ix(4) = mod (ix(3), mn1) + 1
          ix(5) = mod (ix(4), mn1) + 1
          ix(6) = mod (ix(5), mn1) + 1
        endif 

        iy(1) = lhnode(ih0,2) - indy * ly + 1
        iy(2) = iy(1) + 1
        iy(3) = iy(2) + 1
        iy(4) = iy(3) + 1
        iy(5) = iy(4) + 1
        iy(6) = iy(5) + 1

        iz(1) = lhnode(ih0,3) - indz * lz + 1
        iz(2) = iz(1) + 1
        iz(3) = iz(2) + 1
        iz(4) = iz(3) + 1
        iz(5) = iz(4) + 1
        iz(6) = iz(5) + 1
       
        do  ndx = 1, 6
         ix1 = ix(ndx)
         do  ndy = 1, 6
          iy1 = iy(ndy)
          do  ndz = 1, 6
           iz1 = iz(ndz)
           vp(icpt,ih0) = vp(icpt,ih0) + vext(ix1,iy1,iz1) * bg(ih0,ndx,1) * bg(ih0,ndy,2) * bg(ih0,ndz,3)
          enddo
         enddo
        enddo

      END DO


      DEALLOCATE ( v1s, v1r, v2s, v2r, v3s, v3r, v4s, v4r, vext )
      DEALLOCATE ( v5s, v5r, v6s, v6r, v7s, v7r, v8s, v8r )

      etime          = mpi_wtime()
      time_interp_cp = time_inter_cp + etime - stime

      RETURN
      END SUBROUTINE INTERPOLATE 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE NEIGHBORS (id,id1,id2,id3,id4,id5,id6,id7,id8)

!     This subroutine finds the processor id's of 8 neighboring domains.
                                                                                                                                                                       
      USE params
      INCLUDE 'mpif.h'

      INTEGER       :: id, id1, id2, id3, id4, id5, id6, id7, id8

      indy = mod (id, nprocY)
      indz = int ( id / nprocY )

      id1 = indz * nprocY + mod (indy+1, nprocY) 
      id5 = indz * nprocY + mod (indy+nprocY-1, nprocY)

      id3 = mod (indz+1, nprocZ) * nprocY + indy
      id7 = mod (indz+nprocZ-1, nprocZ) * nprocY + indy

      id2 = mod (indz+1, nprocZ) * nprocY + mod (indy+1, nprocY)
      id4 = mod (indz+1, nprocZ) * nprocY + mod (indy+nprocY-1, nprocY)

      id8 = mod (indz+nprocZ-1, nprocZ) * nprocY + mod (indy+1, nprocY)
      id6 = mod (indz+nprocZ-1, nprocZ) * nprocY + mod (indy+nprocY-1, nprocY)

      RETURN
      END SUBROUTINE NEIGHBORS

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

       SUBROUTINE PART_HISTORY_FVP (id)

!      This subroutine updates particle history using fvp. 

       USE params
       INCLUDE 'mpif.h'

!      Note: new_particle will be set to .FALSE. in PART_HISTORY_FDP after all initialization is complete.

       IF ( new_particle ) THEN

        yp(4,1:nps) = Sv_const * yp(7,1:nps)**2
        yp(5,1:nps) = 0.0 
        yp(6,1:nps) = 0.0

!        yp(4,1:nps) = vp(1,1:nps) + Sv_const * yp(7,1:nps)**2
!        yp(5,1:nps) = vp(2,1:nps)
!        yp(6,1:nps) = vp(3,1:nps)

        fvp(:,0,1:nps) = yp(4:6,1:nps)
!        fdp(:,0,1:nps) = vp(:,1:nps)

!        if ( HDI_included ) fdp(:,0,1:nps) = fdp(:,0,1:nps) + pertvel(:,1:nps)

        fvp(:,1,1:nps) = fvp(:,0,1:nps)
        fvp(:,2,1:nps) = fvp(:,0,1:nps)
!        fdp(:,1,1:nps) = fdp(:,0,1:nps)
!        fdp(:,2,1:nps) = fdp(:,0,1:nps)
!        fdp(:,3,1:nps) = fdp(:,0,1:nps)

!        new_particle   = .FALSE.

       ELSE

        fvp(:,0,1:nps) = yp(4:6,1:nps)
!        fdp(:,0,1:nps) = vp(1:3,1:nps)

!        if ( HDI_included ) fdp(:,0,1:nps) = fdp(:,0,1:nps) + pertvel(:,1:nps)
                                      
       ENDIF

      RETURN
      END SUBROUTINE PART_HISTORY_FVP 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                                                                                                                                                             
       SUBROUTINE PART_HISTORY_FDP (id)
                                                                                                                                                             
!      This subroutine updates particle history using fdp.
                                                                                                                                                             
       USE params
       INCLUDE 'mpif.h'
                                                                                                                                                             
       IF ( new_particle ) THEN
                                                                                                                                                             
        fdp(:,0,1:nps) = vp(:,1:nps)
        if ( HDI_included ) fdp(:,0,1:nps) = vp(:,1:nps) + pertvel(:,1:nps)
                                                                                                                                                             
        fdp(:,1,1:nps) = fdp(:,0,1:nps)
        fdp(:,2,1:nps) = fdp(:,0,1:nps)
        fdp(:,3,1:nps) = fdp(:,0,1:nps)
                        
        new_particle   = .FALSE.
                        
       ELSE
           
        fdp(:,0,1:nps) = vp(1:3,1:nps)
        if ( HDI_included ) fdp(:,0,1:nps) = vp(:,1:nps) + pertvel(:,1:nps)
                                  
       ENDIF

      RETURN
      END SUBROUTINE PART_HISTORY_FDP

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE PART_ADVANCE (id)
                                                                                                                                                             
!     This subroutine advances the particle position by solving particle equation of motion.
                                                                                                                                                             
      USE params
      INCLUDE 'mpif.h'
                                                                                                                                                             
      REAL    :: FSN, Vrel

       stime = mpi_wtime()

       DO ih = 1, nps

!       double check implementation of FSN for non-linear drag.                     
        FSN = 1.
                
        IF ( ndrag .eq. 1 ) THEN

         Vrel = sqrt ( (yp(4,ih)-vp(1,ih))**2 + (yp(5,ih)-vp(2,ih))**2 + (yp(6,ih)-vp(3,ih))**2 )
         if ( HDI_included ) then
          Vrel = sqrt ( (yp(4,ih)-vp(1,ih)-pertvel(1,ih))**2 + (yp(5,ih)-vp(2,ih)-pertvel(2,ih))**2 + (yp(6,ih)-vp(3,ih)-pertvel(3,ih))**2 )
         endif

         FSN  = 1. + 0.15 * ( Vrel * 2. * yp(7,ih) / rnu ) ** 0.687

        ENDIF
             
        do ic = 1, 3
                    
         unitv = 0.0
                    
         if ( ic .eq. 1 ) unitv = 1.0
                                   
         taup = tau_CLOUD_const2 * yp(7,ih)**2
         dtau = dt / (24.0*taup) * FSN
                                           
         tempv  = 19.0*fvp(ic,0,ih) - 5.0*fvp(ic,1,ih) + fvp(ic,2,ih)
         tempd  = 55.0*fdp(ic,0,ih) - 59.0*fdp(ic,1,ih) + 37.0*fdp(ic,2,ih) - 9.0*fdp(ic,3,ih)
                                           
         wp_vel = Sv_const * yp(7,ih)**2
                                        
         yp(3+ic,ih) = ( yp(3+ic,ih) + dt * unitv * wp_vel / taup * FSN       &
                     - dtau * tempv + dtau * tempd ) / ( 1.0 + 9.0 * dtau )
                                           
        enddo


       ENDDO

!-------------------------------------------------------------------
!      Advance the particle position
                                           
       dtc = dt / 24.0
                      
       dmove(:,1:nps) = dtc * ( 9.0 * yp(4:6,1:nps) + 19.0 * fvp(:,0,1:nps)  &
                      - 5.0 * fvp(:,1,1:nps) + fvp(:,2,1:nps))

!-------------------------------------------------------------------
!      Update particle history
                                           
!       if ( nps.gt.0 ) then
!         fvp(:,2,1:nps) = fvp(:,1,1:nps)
!         fvp(:,1,1:nps) = fvp(:,0,1:nps)
!         fdp(:,3,1:nps) = fdp(:,2,1:nps)
!         fdp(:,2,1:nps) = fdp(:,1,1:nps)
!         fdp(:,1,1:nps) = fdp(:,0,1:nps)
!       endif

!       open (9834,file='fort.9834',status='unknown',position='append')

!       IF ( istep .eq. 1 ) THEN
!        do ijkl = 1, nps
!         write(9834,1134) id, istep, nps, yp(:,ijkl)
!        enddo
!       ENDIF
!       close (9834)
!1134   format ( 3I8, 8(1pe18.10) )

!       CALL COLLISION_DET (id)

!       if ( nps.gt.0 ) then
!         fvp(:,2,1:nps) = fvp(:,1,1:nps)
!         fvp(:,1,1:nps) = fvp(:,0,1:nps)
!         fdp(:,3,1:nps) = fdp(:,2,1:nps)
!         fdp(:,2,1:nps) = fdp(:,1,1:nps)
!         fdp(:,1,1:nps) = fdp(:,0,1:nps)
!       endif

!-------------------------------------------------------------------
!      Update particle position 

!       yp(1:3,1:nps) = yp(1:3,1:nps) + dmove(1:3,1:nps)

       etime           = mpi_wtime()
       time_partadv_cp = time_partadv_cp + etime - stime
                                        
      RETURN
      END SUBROUTINE PART_ADVANCE 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE UPDATE_PART_POSITION (id)

!     This subroutine updates particle position.

      USE params
      INCLUDE 'mpif.h'

      stime = mpi_wtime()

!     Update particle history

      if ( nps.gt.0 ) then
        fvp(:,2,1:nps) = fvp(:,1,1:nps)
        fvp(:,1,1:nps) = fvp(:,0,1:nps)
        fdp(:,3,1:nps) = fdp(:,2,1:nps)
        fdp(:,2,1:nps) = fdp(:,1,1:nps)
        fdp(:,1,1:nps) = fdp(:,0,1:nps)
      endif

!     Update particle position
 
      yp(1:3,1:nps) = yp(1:3,1:nps) + dmove(1:3,1:nps)
                           
      etime           = mpi_wtime()
      time_partadv_cp = time_partadv_cp + etime - stime
                           
      RETURN
      END SUBROUTINE UPDATE_PART_POSITION
                           
!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE COMPLETE_YPN (id)

!     This subroutine adds dmove to ypn.

      USE params 
      INCLUDE 'mpif.h'

      INTEGER   :: req(2), nsize, error, ierr, partind
      INTEGER   :: status1(MPI_STATUS_SIZE), status2(MPI_STATUS_SIZE,2)

      REAL,    ALLOCATABLE, DIMENSION (:,:,:)    :: temp_p, temp_b 

      stime = mpi_wtime()

      nsize = int ( 20.0 * nuniform )

      ALLOCATE ( temp_p  (6,nsize,8) )
      ALLOCATE ( temp_b  (6,nsize,8) )

      ypn (9:11, 1:nps) = dmove(1:3, 1:nps) 
      ypn (4:6, 1:nps)  = yp(4:6, 1:nps) 

      DO jid = 1, 8 
       DO jik = 1, ih_p(jid) 

        partind              = pindex (jik, jid)
        temp_p (1:3,jik,jid) = dmove(1:3, partind)  
        temp_p (4:6,jik,jid) = yp(4:6, partind)  

       ENDDO
      ENDDO

      etime           = mpi_wtime()
      time_partadv_cp = time_partadv_cp + etime - stime
      stime           = mpi_wtime()

      CALL MPI_BARRIER (MPI_COMM_WORLD, error)

!-----Communicate the particles close to subdomain boundaries

      DO jid = 1, 8
                                                                                                                                                                    
        ids = nei_id ( jid )
        idr = nei_id ( mod (jid+3, 8) + 1 )
                                                                                                                                                                    
        if ( ih_p(jid) .gt. 0 ) then
         CALL MPI_SEND (temp_p(:,1:ih_p(jid),jid),6*ih_p(jid),MPI_REAL8,ids,213,MPI_COMM_WORLD,ierr)
        endif
                                                                                                                                                                    
        if ( ih_pt(jid) .gt. 0 ) then
         CALL MPI_RECV (temp_b(:,1:ih_pt(jid),jid),6*ih_pt(jid),MPI_REAL8,idr,213,MPI_COMM_WORLD,status1,ierr)
        endif
                                                                                                                                                                    
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAITALL (2,req,status2,ierr)
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).eq.0 ) CALL MPI_WAIT (req(1),status1,ierr)
!        if ( ih_p(jid).eq.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAIT (req(2),status1,ierr)
                                                                                                                                                                    
        if ( ih_pt(jid) .gt. 0 ) then
         
          nps1 = nps + SUM ( ih_pt ( 1 : (jid-1) ) )
          nps2 = nps + SUM ( ih_pt ( 1 : jid ) )
         
          ypn (9:11, (nps1+1):nps2) = temp_b(1:3,1:ih_pt(jid),jid)
          ypn (4:6, (nps1+1):nps2)  = temp_b(4:6,1:ih_pt(jid),jid)
         
        endif
         
      ENDDO

      CALL MPI_BARRIER (MPI_COMM_WORLD, error)

      etime           = mpi_wtime()
      time_partadv_cm = time_partadv_cm + etime - stime

      DEALLOCATE ( temp_p, temp_b ) 

      RETURN
      END SUBROUTINE COMPLETE_YPN 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE PERIODIC (id)

!     This subroutine checks the condition on domain boundaries be periodic
!     also takes care of particles' subdomain change.

      USE params 
      INCLUDE 'mpif.h'

      INTEGER                                 :: req(2), ierr, error
      INTEGER                                 :: status2(MPI_STATUS_SIZE,2)
      INTEGER                                 :: status1(MPI_STATUS_SIZE)
      INTEGER                                 :: ids, idr, nid

      REAL, ALLOCATABLE, DIMENSION (:,:,:)    :: temp_p, temp_b 

      stime = mpi_wtime()

      nsize = int ( 10 * nuniform )

      ALLOCATE ( temp_p (26,nsize,8) )
      ALLOCATE ( temp_b (26,nsize,8) )

      ih0   = 1
      ijk   = 0 
      ih_p  = 0
      ih_pt = 0

!      CALL MPI_BARRIER (MPI_COMM_WORLD, error)

      IF ( nps.gt.0 ) THEN

156     CONTINUE

        IF ( ( ih0 + SUM (ih_p(1:8)) ) .gt. nps ) GOTO 158 

        if ( (yp(2,ih0).gt.(real(indy+1)*hy)) .AND. (yp(3,ih0).le.(real(indz+1)*hz)) .AND. (yp(3,ih0).ge.(real(indz)*hz)) ) then

           nid       = 1  
           ih_p(nid) = ih_p(nid) + 1

           temp_p (1:8, ih_p(nid), nid)   = yp(:,ih0)
           temp_p (9:11, ih_p(nid), nid)  = fvp(:,1,ih0)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,2,ih0)
           temp_p (15:17, ih_p(nid), nid) = fdp(:,1,ih0)
           temp_p (18:20, ih_p(nid), nid) = fdp(:,2,ih0)
           temp_p (21:23, ih_p(nid), nid) = fdp(:,3,ih0)
           temp_p (24:26, ih_p(nid), nid) = dmove(:,ih0)

           ijk = nps - SUM ( ih_p(1:8) ) + 1

           yp(:,ih0)      = yp(:,ijk)
           fvp(:,1:2,ih0) = fvp(:,1:2,ijk)
!           fvp(:,0:2,ih0) = fvp(:,0:2,ijk)
           fdp(:,1:3,ih0) = fdp(:,1:3,ijk)
!           fdp(:,0:3,ih0) = fdp(:,0:3,ijk)
           dmove(:,ih0)   = dmove(:,ijk)

           GOTO 156
        endif

        if ( (yp(2,ih0).gt.(real(indy+1)*hy)) .AND. (yp(3,ih0).gt.(real(indz+1)*hz)) ) then
                                     
           nid       = 2 
           ih_p(nid) = ih_p(nid) + 1
                                                                                                                                                                         
           temp_p (1:8, ih_p(nid), nid)   = yp(:,ih0)
           temp_p (9:11, ih_p(nid), nid)  = fvp(:,1,ih0)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,2,ih0)
           temp_p (15:17, ih_p(nid), nid) = fdp(:,1,ih0)
           temp_p (18:20, ih_p(nid), nid) = fdp(:,2,ih0)
           temp_p (21:23, ih_p(nid), nid) = fdp(:,3,ih0)
           temp_p (24:26, ih_p(nid), nid) = dmove(:,ih0)

           ijk = nps - SUM ( ih_p(1:8) ) + 1
                                                                                                                                                                         
           yp(:,ih0)      = yp(:,ijk)
           fvp(:,1:2,ih0) = fvp(:,1:2,ijk)
!           fvp(:,0:2,ih0) = fvp(:,0:2,ijk)
           fdp(:,1:3,ih0) = fdp(:,1:3,ijk)
!          fdp(:,0:3,ih0) = fdp(:,0:3,ijk)
           dmove(:,ih0)   = dmove(:,ijk)
                                                                                                                                                                         
           GOTO 156
        endif

        if ( (yp(2,ih0).ge.(real(indy)*hy)) .AND. (yp(2,ih0).le.(real(indy+1)*hy)) .AND. (yp(3,ih0).gt.(real(indz+1)*hz)) ) then

           nid       = 3 
           ih_p(nid) = ih_p(nid) + 1
                                                                                                                                                                         
           temp_p (1:8, ih_p(nid), nid)   = yp(:,ih0)
           temp_p (9:11, ih_p(nid), nid)  = fvp(:,1,ih0)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,2,ih0)
           temp_p (15:17, ih_p(nid), nid) = fdp(:,1,ih0)
           temp_p (18:20, ih_p(nid), nid) = fdp(:,2,ih0)
           temp_p (21:23, ih_p(nid), nid) = fdp(:,3,ih0)
           temp_p (24:26, ih_p(nid), nid) = dmove(:,ih0)
                                                                                                                                                              
           ijk = nps - SUM ( ih_p(1:8) ) + 1
                                                                                                                                                                         
           yp(:,ih0)      = yp(:,ijk)
           fvp(:,1:2,ih0) = fvp(:,1:2,ijk)
!           fvp(:,0:2,ih0) = fvp(:,0:2,ijk)
           fdp(:,1:3,ih0) = fdp(:,1:3,ijk)
!           fdp(:,0:3,ih0) = fdp(:,0:3,ijk)
           dmove(:,ih0)   = dmove(:,ijk)

           GOTO 156
        endif

        if ( (yp(2,ih0).lt.(real(indy)*hy)) .AND. (yp(3,ih0).gt.(real(indz+1)*hz)) ) then
                                                                                                                                                                         
           nid       = 4 
           ih_p(nid) = ih_p(nid) + 1
                                                                                                                                                                         
           temp_p (1:8, ih_p(nid), nid)   = yp(:,ih0)
           temp_p (9:11, ih_p(nid), nid)  = fvp(:,1,ih0)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,2,ih0)
           temp_p (15:17, ih_p(nid), nid) = fdp(:,1,ih0)
           temp_p (18:20, ih_p(nid), nid) = fdp(:,2,ih0)
           temp_p (21:23, ih_p(nid), nid) = fdp(:,3,ih0)
           temp_p (24:26, ih_p(nid), nid) = dmove(:,ih0)

           ijk = nps - SUM ( ih_p(1:8) ) + 1

           yp(:,ih0)      = yp(:,ijk)
           fvp(:,1:2,ih0) = fvp(:,1:2,ijk)
!           fvp(:,0:2,ih0) = fvp(:,0:2,ijk)
           fdp(:,1:3,ih0) = fdp(:,1:3,ijk)
!           fdp(:,0:3,ih0) = fdp(:,0:3,ijk)
           dmove(:,ih0)   = dmove(:,ijk)
                                                                                                                                                                         
           GOTO 156
        endif

        if ( (yp(2,ih0).lt.(real(indy)*hy)) .AND. (yp(3,ih0).le.(real(indz+1)*hz)) .AND. (yp(3,ih0).ge.(real(indz)*hz)) ) then

           nid       = 5 
           ih_p(nid) = ih_p(nid) + 1
                                                                                                                                                                         
           temp_p (1:8, ih_p(nid), nid)   = yp(:,ih0)
           temp_p (9:11, ih_p(nid), nid)  = fvp(:,1,ih0)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,2,ih0)
           temp_p (15:17, ih_p(nid), nid) = fdp(:,1,ih0)
           temp_p (18:20, ih_p(nid), nid) = fdp(:,2,ih0)
           temp_p (21:23, ih_p(nid), nid) = fdp(:,3,ih0)
           temp_p (24:26, ih_p(nid), nid) = dmove(:,ih0)
                                                                                                                                                                         
           ijk = nps - SUM ( ih_p(1:8) ) + 1
                                                                                                                                                                         
           yp(:,ih0)      = yp(:,ijk)
           fvp(:,1:2,ih0) = fvp(:,1:2,ijk)
!           fvp(:,0:2,ih0) = fvp(:,0:2,ijk)
           fdp(:,1:3,ih0) = fdp(:,1:3,ijk)
!           fdp(:,0:3,ih0) = fdp(:,0:3,ijk)
           dmove(:,ih0)   = dmove(:,ijk)
                                                                                                                                                                         
           GOTO 156
        endif

        if ( (yp(2,ih0).lt.(real(indy)*hy)) .AND. (yp(3,ih0).lt.(real(indz)*hz)) ) then
                                
           nid       = 6 
           ih_p(nid) = ih_p(nid) + 1
                                                                                                                                                                         
           temp_p (1:8, ih_p(nid), nid)   = yp(:,ih0)
           temp_p (9:11, ih_p(nid), nid)  = fvp(:,1,ih0)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,2,ih0)
           temp_p (15:17, ih_p(nid), nid) = fdp(:,1,ih0)
           temp_p (18:20, ih_p(nid), nid) = fdp(:,2,ih0)
           temp_p (21:23, ih_p(nid), nid) = fdp(:,3,ih0)
           temp_p (24:26, ih_p(nid), nid) = dmove(:,ih0)
 
           ijk = nps - SUM ( ih_p(1:8) ) + 1
                                  
           yp(:,ih0)      = yp(:,ijk)
           fvp(:,1:2,ih0) = fvp(:,1:2,ijk)
!           fvp(:,0:2,ih0) = fvp(:,0:2,ijk)
           fdp(:,1:3,ih0) = fdp(:,1:3,ijk)
!           fdp(:,0:3,ih0) = fdp(:,0:3,ijk)
           dmove(:,ih0)   = dmove(:,ijk)
                                                                                                                                                                         
           GOTO 156
        endif

        if ( (yp(3,ih0).lt.(real(indz)*hz)) .AND. (yp(2,ih0).le.(real(indy+1)*hy)) .AND. (yp(2,ih0).ge.(real(indy)*hy)) ) then
                                                                                                                                                                         
           nid       = 7 
           ih_p(nid) = ih_p(nid) + 1
                                                                                                                                                                         
           temp_p (1:8, ih_p(nid), nid)   = yp(:,ih0)
           temp_p (9:11, ih_p(nid), nid)  = fvp(:,1,ih0)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,2,ih0)
           temp_p (15:17, ih_p(nid), nid) = fdp(:,1,ih0)
           temp_p (18:20, ih_p(nid), nid) = fdp(:,2,ih0)
           temp_p (21:23, ih_p(nid), nid) = fdp(:,3,ih0)
           temp_p (24:26, ih_p(nid), nid) = dmove(:,ih0)

           ijk = nps - SUM ( ih_p(1:8) ) + 1
                                                                                                                                                                         
           yp(:,ih0)      = yp(:,ijk)
           fvp(:,1:2,ih0) = fvp(:,1:2,ijk)
!           fvp(:,0:2,ih0) = fvp(:,0:2,ijk)
           fdp(:,1:3,ih0) = fdp(:,1:3,ijk)
!           fdp(:,0:3,ih0) = fdp(:,0:3,ijk)
           dmove(:,ih0)   = dmove(:,ijk)
                                                                                                                                                                         
           GOTO 156
        endif

        if ( (yp(3,ih0).lt.(real(indz)*hz)) .AND. (yp(2,ih0).gt.(real(indy+1)*hy)) ) then
                                                                                                                                                                         
           nid       = 8 
           ih_p(nid) = ih_p(nid) + 1
                                                                                                                                                                         
           temp_p (1:8, ih_p(nid), nid)   = yp(:,ih0)
           temp_p (9:11, ih_p(nid), nid)  = fvp(:,1,ih0)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,2,ih0)
           temp_p (15:17, ih_p(nid), nid) = fdp(:,1,ih0)
           temp_p (18:20, ih_p(nid), nid) = fdp(:,2,ih0)
           temp_p (21:23, ih_p(nid), nid) = fdp(:,3,ih0)
           temp_p (24:26, ih_p(nid), nid) = dmove(:,ih0)

           ijk = nps - SUM ( ih_p(1:8) ) + 1
                                                                                                                                                                         
           yp(:,ih0)      = yp(:,ijk)
           fvp(:,1:2,ih0) = fvp(:,1:2,ijk)
!           fvp(:,0:2,ih0) = fvp(:,0:2,ijk)
           fdp(:,1:3,ih0) = fdp(:,1:3,ijk)
!           fdp(:,0:3,ih0) = fdp(:,0:3,ijk)
           dmove(:,ih0)   = dmove(:,ijk)
                                                                                                                                                                         
           GOTO 156
        endif

        ih0 = ih0 + 1

        IF ( ( ih0 + SUM (ih_p(1:8)) ) .gt. nps ) GOTO 158 
        GOTO 156

      ENDIF

158   CONTINUE

!     Actual number of particles in the pencil. SUM(ih_p(1:8)) particles are...
!     going to be in the neighbouring subdomains next couple of lines.   

      nps = nps - SUM ( ih_p(1:8) )

!-------------------------------------------------------------------

       etime            = mpi_wtime()
       time_periodic_cp = time_periodic_cp + etime - stime
       stime            = mpi_wtime()

       DO jid = 1, 8

        ids = nei_id ( jid )
        idr = nei_id ( mod (jid+3, 8) + 1 )

        CALL MPI_SENDRECV (ih_p(jid),1,MPI_INTEGER,ids,212,ih_pt(jid),1,MPI_INTEGER,idr,212,MPI_COMM_WORLD,status1,ierr)
                                                                                                                                                                         
        if ( ih_p(jid) .gt. 0 ) then
         CALL MPI_SEND (temp_p(:,1:ih_p(jid),jid),26*ih_p(jid),MPI_REAL8,ids,213,MPI_COMM_WORLD,ierr)
        endif
                                                                                                                                                                         
        if ( ih_pt(jid) .gt. 0 ) then
         CALL MPI_RECV (temp_b(:,1:ih_pt(jid),jid),26*ih_pt(jid),MPI_REAL8,idr,213,MPI_COMM_WORLD,status1,ierr)
        endif
                                                                                                                                                                         
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAITALL (2,req,status2,ierr)
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).eq.0 ) CALL MPI_WAIT (req(1),status1,ierr)
!        if ( ih_p(jid).eq.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAIT (req(2),status1,ierr)

       if ( ih_pt(jid) .gt. 0 ) then

         nps1 = nps + 1
         nps2 = nps + ih_pt(jid)
                          
         yp(:,nps1:nps2)    = temp_b(1:8,1:ih_pt(jid),jid)
         fvp(:,1,nps1:nps2) = temp_b(9:11,1:ih_pt(jid),jid)
         fvp(:,2,nps1:nps2) = temp_b(12:14,1:ih_pt(jid),jid)
         fdp(:,1,nps1:nps2) = temp_b(15:17,1:ih_pt(jid),jid)
         fdp(:,2,nps1:nps2) = temp_b(18:20,1:ih_pt(jid),jid)
         fdp(:,3,nps1:nps2) = temp_b(21:23,1:ih_pt(jid),jid)
         dmove(:,nps1:nps2) = temp_b(24:26,1:ih_pt(jid),jid)
                                              
         nps = nps2

        endif

       ENDDO

       etime            = mpi_wtime()
       time_periodic_cm = time_periodic_cm + etime - stime
       stime            = mpi_wtime()

!-------------------------------------------------------------------
!     Relocate particles considering 3D periodicity

       DO ih = 1, nps

          if ( yp(1,ih) .ge. hx  ) yp(1,ih) = yp(1,ih) - hx
          if ( yp(1,ih) .lt. 0.0 ) yp(1,ih) = yp(1,ih) + hx
          if ( yp(2,ih) .ge. pi2 ) yp(2,ih) = yp(2,ih) - pi2
          if ( yp(2,ih) .lt. 0.0 ) yp(2,ih) = yp(2,ih) + pi2
          if ( yp(3,ih) .ge. pi2 ) yp(3,ih) = yp(3,ih) - pi2
          if ( yp(3,ih) .lt. 0.0 ) yp(3,ih) = yp(3,ih) + pi2

       ENDDO

!       CALL MPI_BARRIER (MPI_COMM_WORLD, error)

       DEALLOCATE ( temp_p )
       DEALLOCATE ( temp_b )

       etime            = mpi_wtime()
       time_periodic_cp = time_periodic_cp + etime - stime

       RETURN
       END SUBROUTINE PERIODIC

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE COLLISION_DET (id)

!     This subroutine detects particle collisions and computes RDF & Wr.
!     After detection, both colliding particles are relocated with the same material properties.

      USE params
      INCLUDE 'mpif.h'

      REAL      :: xa(3), xan(3), va(3), van(3), ax1(3), xta(3)
      REAL      :: xb(3), xbn(3), vb(3), vbn(3), bx1(3), xtb(3)
      REAL      :: xd(3)

      INTEGER   :: req(2), error, ierr
      INTEGER   :: status2(MPI_STATUS_SIZE,2)
      INTEGER   :: status1(MPI_STATUS_SIZE)

      INTEGER, ALLOCATABLE, DIMENSION (:,:)    :: colltype, colltype_t
      INTEGER, ALLOCATABLE, DIMENSION (:)      :: ncollpart_t
      INTEGER, ALLOCATABLE, DIMENSION (:,:)    :: removepart_t
      INTEGER, ALLOCATABLE, DIMENSION (:,:)    :: nir, nir_t
      INTEGER, ALLOCATABLE, DIMENSION (:,:,:)  :: nirbin0, nirbin_t

      REAL,    ALLOCATABLE, DIMENSION (:,:,:)  :: wrtan0, wrtan_t
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)  :: thetabin0, thetabin_t
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)  :: wxbin0, wxbin_t
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)  :: wybin0, wybin_t
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)  :: wzbin0, wzbin_t
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)  :: wrbin0, wrbin_t
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)  :: urbin0, urbin_t
      REAL,    ALLOCATABLE, DIMENSION (:,:)    :: buffer
      REAL,    ALLOCATABLE, DIMENSION (:,:)    :: bufferrecv

      mevent = int ( 10 * nuniform )
      nsize  = int ( 10 * nuniform )

      ALLOCATE ( colltype   (nset,nset)       )
      ALLOCATE ( colltype_t (nset,nset)       )
      ALLOCATE ( nir        (nset,nset)       )
      ALLOCATE ( nir_t      (nset,nset)       )
      ALLOCATE ( buffer     (12,mevent)       )

      ALLOCATE ( nirbin0    (nset,nset,0:nbin-1) )
      ALLOCATE ( nirbin_t   (nset,nset,0:nbin-1) )
      ALLOCATE ( wxbin0     (nset,nset,0:nbin-1) )
      ALLOCATE ( wxbin_t    (nset,nset,0:nbin-1) )
      ALLOCATE ( wybin0     (nset,nset,0:nbin-1) )
      ALLOCATE ( wybin_t    (nset,nset,0:nbin-1) )
      ALLOCATE ( wzbin0     (nset,nset,0:nbin-1) )
      ALLOCATE ( wzbin_t    (nset,nset,0:nbin-1) )
      ALLOCATE ( wrbin0     (nset,nset,0:nbin-1) )
      ALLOCATE ( wrbin_t    (nset,nset,0:nbin-1) )
      ALLOCATE ( urbin0     (nset,nset,0:nbin-1) )
      ALLOCATE ( urbin_t    (nset,nset,0:nbin-1) )
      ALLOCATE ( thetabin0  (nset,nset,0:nbin-1) )
      ALLOCATE ( thetabin_t (nset,nset,0:nbin-1) )
      ALLOCATE ( wrtan0     (nset,nset,0:nbin-1) )
      ALLOCATE ( wrtan_t    (nset,nset,0:nbin-1) )

      ALLOCATE ( ncollpart_t(1:4)             )
      ALLOCATE ( removepart_t(nsize,1:4)      )

      stime   = mpi_wtime()

      ij      = 0
      ijl     = 0
      ncoll   = 0
      nir     = 0
      nir_t   = 0

      ntotpair  = 0
      nirbin0   = 0
      wrbin0    = 0.0
      wxbin0    = 0.0
      wybin0    = 0.0
      wzbin0    = 0.0
      urbin0    = 0.0
      thetabin0 = 0.0
      wrtan0    = 0.0

      colltype    = 0
      ncollpart   = 0
      removepart  = 0
      ncollpart_t = 0
!     removepart and ncollpart are global arrays allocated in PART_INIT. 

!     number of colision of each type.
      icc1 = 0
      icc2 = 0
      icc3 = 0
  
      nps4 = nps + SUM ( ih_pt (1:4) )
      nps5 = nps + SUM ( ih_pt (1:5) )
      nps6 = nps + SUM ( ih_pt (1:6) )
      nps7 = nps + SUM ( ih_pt (1:7) )
      nps8 = nps + SUM ( ih_pt (1:8) )

!-----For each particle ih1 inside the pencil:

      d_stat  = 20.0 * radmax 
      imax    = 1 + int ( d_stat * real(ncdx) /  hx )
      jmax    = 1 + int ( d_stat * real(ncdy) /  hy )
      kmax    = 1 + int ( d_stat * real(ncdz) /  hz )

      DO ih1 = 1, nps

!       ncollip1 is the total number of the collision events for one target...
!       particle (one target particle may collide with couple of particles not just one).

        ncollip1 = 0

!        dist_test = sqrt ( (yp(1,ih1)-0.14)**2 + (yp(2,ih1)-1.2)**2 + (yp(3,ih1)-2.42)**2 )
!        if ( dist_test .le. 0.20 ) write (*,*) 'istep, ih1, yp(1:8,ih1) = ', istep, ih1, yp(1:8,ih1) 

!        if ( yp(8,ih1) .eq. 163305. ) write (*,*) 'id,istep,yp(1:8,ih1),fvp(1:3,0,ih1) = ', id,istep,yp(1:8,ih1),fvp(1:3,0,ih1) 

!       old and new location of the first particle
        xa(:)  = yp(1:3,ih1)
        xan(:) = xa(:) + dmove(:,ih1)

!       old velocity of the first particle
        va(:) = fvp(:, 0, ih1)

        ix = 1 + int ( xa(1) / hx * real(ncdx) )
        iy = 1 + int ( ( xa(2) - real(indy) * hy ) / hy * real(ncdy) )
        iz = 1 + int ( ( xa(3) - real(indz) * hz ) / hz * real(ncdz) )
                                                                                                                                                                      
        if ( ix .gt. ncdx ) ix = ncdx
        if ( iy .gt. ncdy ) iy = ncdy
        if ( iz .gt. ncdz ) iz = ncdz

        DO i = -imax, imax
         DO j = -jmax, jmax
          DO k = -kmax, kmax 

!-------   Find the second particle ih2 in neighbourhood of particle ih1.

           iy2 = iy + j
           IF ( ( iy2 .LT. (1-ny_coll) ) .OR. ( iy2 .GT. (ncdy+ny_coll) ) ) GOTO 201

           iz2 = iz + k
           IF ( ( iz2 .LT. (1-nz_coll) ) .OR. ( iz2 .GT. (ncdz+nz_coll) ) ) GOTO 201

!          eliminate double counting of collisions.
           IF ( iz2 .LT. 1 ) GOTO 201
           IF ( ( iz2 .LE. ncdz ) .AND. ( iy2 .LT. 1 ) ) GOTO 201

           ix2 = ix + i
           sx  = 0.0

           if ( ix2.lt.1 ) then
            ix2 = ncdx + ix2
            sx  = -1.
           elseif ( ix2.gt.ncdx ) then
            ix2 = mod (ix2, ncdx) 
            sx  = 1.
           endif

!          2nd particle's index
           ih2 = head (ix2,iy2,iz2)
           IF ( ( ih2 .LE. ih1 ) .OR. ( ih2 .EQ. 0 ) ) GOTO 201

199        CONTINUE

!          old and new location of the second droplet
           xb(1)  = ypn(1,ih2) + sx*hx
           xb(2)  = ypn(2,ih2)
           xb(3)  = ypn(3,ih2)

           xbn(:) = xb(:) + ypn(9:11,ih2)

!          old velocity of the second particle
           vb(:) = ypn(12:14,ih2) 

           xd    = xa - xb
           dnij0 = sqrt ( xd(1)**2 + xd(2)**2 + xd(3)**2 )

!-------   define the collision radius and inner & outer shell radius.

!           rcoll = ypn(7,ih1) + ypn(7,ih2)
           rcoll = yp(7,ih1) + ypn(7,ih2)

           if ( HDI_included ) then
            rcolla = rcoll
           else
            rcolla = rcoll * ( 1. - st/100. )
           endif

           rcollb = rcoll * ( 1. + st/100. )

!-------------------------------------------------------------------
!  RDF(r)

           IF ( (dnij0.GT.rcolla) .AND. (dnij0.LE.10*rcoll) ) THEN

            if ( dnij0.GT.rcoll ) then

             xd  = xd / dnij0
             Vr1 = DOT_PRODUCT (va,xd)
             Vr2 = DOT_PRODUCT (vb,xd)
             Wr  = Vr1 - Vr2

!             ir  = int ( (dnij0-rcoll) / (0.05*rcoll) )
             ir  = int ( (dnij0-rcoll) / ((9./real(nbin))*rcoll) )
   
             if ( ir .gt. (nbin-1) ) ir = nbin-1 
             if ( ir .lt. 0   )      ir = 0

             CALL FINDINDEX (yp(7,ih1),ypn(7,ih2),index1,index2)

             nirbin0 (index1,index2,ir) = nirbin0 (index1,index2,ir) + 1

             if ( index2 .eq. index1 ) then 
              wxbin0 (index1,index2,ir) = wxbin0 (index1,index2,ir) + va(1) - vb(1)
              wybin0 (index1,index2,ir) = wybin0 (index1,index2,ir) + va(2) - vb(2)
              wzbin0 (index1,index2,ir) = wzbin0 (index1,index2,ir) + va(3) - vb(3)
             endif

             if ( yp(7,ih1) .gt. ypn(7,ih2) ) then 
              wxbin0 (index1,index2,ir) = wxbin0 (index1,index2,ir) + va(1) - vb(1)
              wybin0 (index1,index2,ir) = wybin0 (index1,index2,ir) + va(2) - vb(2)
              wzbin0 (index1,index2,ir) = wzbin0 (index1,index2,ir) + va(3) - vb(3)
             elseif ( yp(7,ih1) .lt. ypn(7,ih2) ) then
              wxbin0 (index1,index2,ir) = wxbin0 (index1,index2,ir) + vb(1) - va(1)
              wybin0 (index1,index2,ir) = wybin0 (index1,index2,ir) + vb(2) - va(2)
              wzbin0 (index1,index2,ir) = wzbin0 (index1,index2,ir) + vb(3) - va(3)
             endif

!            wrbin0 (index1,index2,ir)  = wrbin0 (index1,index2,ir) + abs(Wr)
             wrbin0 (index1,index2,ir)  = wrbin0 (index1,index2,ir) + Wr
!            urbin0 (index1,index2,ir)  = urbin0 (index1,index2,ir) + abs ( DOT_PRODUCT (vpn(:,ih1)-vpn(:,ih2),xd) )
             urbin0 (index1,index2,ir)  = urbin0 (index1,index2,ir) + DOT_PRODUCT (vpn(:,ih1)-vpn(:,ih2),xd)

             aaa    = DOT_PRODUCT (va-vb,xd)  
             bbb    = sqrt ( (va(1)-vb(1))**2 + (va(2)-vb(2))**2 + (va(3)-vb(3))**2 )
             vangle = acos ( abs(aaa/bbb) ) 
             if ( bbb .ne. 0. ) thetabin0 (index1,index2,ir) = thetabin0 (index1,index2,ir) + vangle 

             if ( bbb .ne. 0. ) ccc = ( bbb * sin(vangle) )**2 / dnij0 
             if ( bbb .ne. 0. ) wrtan0 (index1,index2,ir) = wrtan0 (index1,index2,ir) + ccc

            endif 

!-------------------------------------------------------------------
!  RDF and relative velocity at contact

            if ( dnij0 .LE. rcollb ) then

             if ( dnij0 .LE. rcoll ) then

              xd  = xd / dnij0
              Vr1 = DOT_PRODUCT (va,xd)
              Vr2 = DOT_PRODUCT (vb,xd)
              Wr  = Vr1 - Vr2

              CALL FINDINDEX (yp(7,ih1),ypn(7,ih2),index1,index2)
             endif 

!             if ( ( abs(Wr) .gt. 1.0 ) .and. ( index1 .eq. 1 ) .and. ( index2 .eq. 1 ) ) then
!              write (*,*) 'id,istep,ih1,yp(8,ih1),ih2,ypn(8,ih2),nps,dnij0,va(:),vb(:) =',id,istep,ih1,yp(8,ih1),ih2,ypn(8,ih2),nps,dnij0,va(:),vb(:)
!             endif

!             write (*,*) 'id,ypn(1:8,ih1),pertvel(1:3,ih1),ypn(1:8,ih2),pertvel(1:3,ih2),rcolla,dnij0,rcollb =',id,ypn(1:8,ih1),pertvel(1:3,ih1),ypn(1:8,ih2),pertvel(1:3,ih2),rcolla,dnij0,rcollb 

             nir (index1,index2) = nir (index1,index2) + 1

             ij = ij + 1
             buffer(1,ij)    = real(index1)
             buffer(2,ij)    = real(index2)
             buffer(3:5,ij)  = va(:)
             buffer(6:8,ij)  = vb(:)
             buffer(9:11,ij) = xd(:)
             buffer(12,ij)   = Wr

            endif

           ENDIF        

!-------------------------------------------------------------------
!          restrict collision counting to neighboring cells.
           IF ( ( abs(i) .GT. 1 ) .OR. ( abs(j) .GT. 1 ) .OR. ( abs(k) .GT. 1 ) ) GOTO 201

!-------------------------------------------------------------------
!  coarse check

!           dmove_a = sqrt ( ypn(9,ih1)**2 + ypn(10,ih1)**2 + ypn(11,ih1)**2 )
           dmove_a = sqrt ( dmove(1,ih1)**2 + dmove(2,ih1)**2 + dmove(3,ih1)**2 )
           dmove_b = sqrt ( ypn(9,ih2)**2 + ypn(10,ih2)**2 + ypn(11,ih2)**2 )

           tmove = rcoll + dmove_a + dmove_b

           IF ( dnij0 .GT. tmove ) GOTO 91

           xd    = xan - xbn 
           dnij1 = sqrt ( xd(1)**2 + xd(2)**2 + xd(3)**2 )

           IF ( (dnij0 .LT. rcoll) .AND. (dnij1 .GT. rcoll) ) GOTO 91

! searching for type I collisions
           if ( (dnij0.GT.rcoll) .AND. (dnij1.LE.rcoll) ) then
            icc1     = icc1 + 1
            ncollip1 = ncollip1 + 1

!            if ( (yp(7,ih1).eq.radius(1)) .AND. (ypn(7,ih2).eq.radius(1)) ) then
!             write (*,*) 'type 1: id,yp(1:8,ih1),ypn(1:8,ih2),dnij0,rcoll,dnij1 =',id,yp(1:8,ih1),ypn(1:8,ih2),dnij0,rcoll,dnij1
!             write (*,*) 'type 1: id,fvp(:,0,ih1),fvp(:,1,ih1),fvp(:,2,ih1),ypn(12:14,ih2) =',id,fvp(:,0,ih1),fvp(:,1,ih1),fvp(:,2,ih1),ypn(12:14,ih2)
!             write (*,*) 'type 1: id,dmove(:,ih1),ypn(9:11,ih2) =',id,dmove(:,ih1),ypn(9:11,ih2)
!            endif

            GOTO 92
           endif

           ax1(:) = dmove(:,ih1)
           van(:) = yp(4:6,ih1)
!           ax1(:) = ypn(9:11,ih1)
!           van(:) = ypn(4:6,ih1)

           bx1(:) = ypn(9:11,ih2)
           vbn(:) = ypn(4:6,ih2)

           do icheck = 1, 19 

            tfac  = 0.05 * real(icheck)
            tfac1 = dt * tfac
            tfac2 = tfac**2
            tfac3 = tfac2 * tfac

            xta = xa + va*tfac1 + (3.*ax1-(2.*va+van)*dt)*tfac2 + ((van+va)*dt-2.*ax1)*tfac3
            xtb = xb + vb*tfac1 + (3.*bx1-(2.*vb+vbn)*dt)*tfac2 + ((vbn+vb)*dt-2.*bx1)*tfac3

            xd   = xta - xtb 
            dnij = sqrt ( xd(1)**2 + xd(2)**2 + xd(3)**2 )
 
! searching for type II collisions
            if ( (dnij0.GT.rcoll) .AND. (dnij.LT.rcoll) ) then
             icc2     = icc2 + 1
             ncollip1 = ncollip1 + 1

             GOTO 92
            endif

! Searching for type III collisions
            if ( (dnij1.LE.rcoll) .AND. (dnij.GT.rcoll) ) then
             icc3     = icc3 + 1
             ncollip1 = ncollip1 + 1

             GOTO 92
            endif 

           enddo 

           GOTO 91

92         CONTINUE

!----      Count the collision

           CALL FINDINDEX (yp(7,ih1),ypn(7,ih2),index1,index2)
           colltype(index1,index2) = colltype(index1,index2) + 1

!          list ih2 particles based on their host subdomain and original particle index.
!          there should be 4 list for neighboring subdomains. 

           if ( ih2 .le. nps )  then

            ncollpart(0)                = ncollpart(0) + 1
            removepart (ncollpart(0),0) = ih2

           elseif ( ( ih2 .gt. nps4 ) .and. ( ih2 .le. nps5 ) )  then

            ncollpart(1) = ncollpart(1) + 1
            removepart (ncollpart(1),1) = int ( ypn(15,ih2) )

           elseif ( ( ih2 .gt. nps5 ) .and. ( ih2 .le. nps6 ) )  then

            ncollpart(2) = ncollpart(2) + 1
            removepart (ncollpart(2),2) = int ( ypn(15,ih2) )

           elseif ( ( ih2 .gt. nps6 ) .and. ( ih2 .le. nps7 ) )  then

            ncollpart(3) = ncollpart(3) + 1
            removepart (ncollpart(3),3) = int ( ypn(15,ih2) )

           elseif ( ( ih2 .gt. nps7 ) .and. ( ih2 .le. nps8 ) )  then

            ncollpart(4) = ncollpart(4) + 1
            removepart (ncollpart(4),4) = int ( ypn(15,ih2) )

           endif

           ntotpair = ntotpair + 1

           collist (1:8, ntotpair)   = yp(1:8,ih1)
           collist (9, ntotpair)     = real ( ih1 )
           collist (10, ntotpair)    = real ( id )
           collist (11:18, ntotpair) = ypn(1:8,ih2)
           collist (19, ntotpair)    = real ( ih2 )
           collist (20, ntotpair)    = real ( id )

!----      go for the next particle in the list 

91         ih2 = list (ih2)
           IF ( (ih2.LE.ih1) .OR. (ih2.EQ.0) ) GOTO 201
           GOTO 199

201        CONTINUE

          ENDDO
         ENDDO
        ENDDO

!----   consider the target particle itself (particle index ih1)

        if ( ncollip1 .NE. 0 ) then

         ncollpart(0)                = ncollpart(0) + 1
         removepart (ncollpart(0),0) = ih1

        endif 

       ENDDO

       etime           = mpi_wtime()
       time_colldet_cp = time_colldet_cp + etime - stime
       stime           = mpi_wtime()

!-------------------------------------------------------------------
! SEND THE RESULTS BACK

       DO jid = 1, 4 
                                                                                                                                                                      
        ids = nei_id ( jid )
        idr = nei_id ( mod (jid+3, 8) + 1 )
                                                                                                                                                                      
        CALL MPI_SENDRECV (ncollpart(jid),1,MPI_INTEGER,ids,212,ncollpart_t(jid),1,MPI_INTEGER,idr,212,&
                       MPI_COMM_WORLD,status1,ierr)
                                                                                                                                                                      
        if ( ncollpart(jid) .gt. 0 ) then
         CALL MPI_SEND (removepart(1:ncollpart(jid),jid),ncollpart(jid),MPI_INTEGER,ids,213,MPI_COMM_WORLD,ierr)
        endif
                                                                                                                                                                      
        if ( ncollpart_t(jid) .gt. 0 ) then
         CALL MPI_RECV (removepart_t(1:ncollpart_t(jid),jid),ncollpart_t(jid),MPI_INTEGER,idr,213,MPI_COMM_WORLD,status1,ierr)
        endif
                                                                                                                                                                      
!        if ( ncollpart(jid).gt.0 .and. ncollpart_t(jid).gt.0 ) CALL MPI_WAITALL (2,req,status2,ierr)
!        if ( ncollpart(jid).gt.0 .and. ncollpart_t(jid).eq.0 ) CALL MPI_WAIT (req(1),status1,ierr)
!        if ( ncollpart(jid).eq.0 .and. ncollpart_t(jid).gt.0 ) CALL MPI_WAIT (req(2),status1,ierr)
                                                                                                                                                                      
        if ( ncollpart_t(jid) .gt. 0 ) then
                                                                                                                                                                      
          ncoll1 = ncollpart(0) + SUM ( ncollpart_t ( 1 : (jid-1) ) )
          ncoll2 = ncollpart(0) + SUM ( ncollpart_t ( 1 : jid ) )
                                                                                                                                                                      
          removepart ((ncoll1+1):ncoll2, 0) = removepart_t(1:ncollpart_t(jid),jid)
                                                                                                                                                                      
         endif
                                                                                                                                                                      
        ENDDO

        etime           = mpi_wtime()
        time_colldet_cm = time_colldet_cm + etime - stime
        stime           = mpi_wtime()

        ncoll = ncollpart(0) + SUM ( ncollpart_t (1:4) )

        if ( ncoll .gt. 0 ) CALL hpsort (removepart(:,0), ncoll)

!        if ( ncollpart(0) .gt. 0 ) write(*,*) 'ncollpart(0), ncoll, removepart(1:ncoll,0) = ', ncollpart(0), ncoll, removepart(1:ncoll,0) 

        etime           = mpi_wtime()
        time_colldet_cp = time_colldet_cp + etime - stime

!-------------------------------------------------------------------
! SAVE DATA FOR RDF(r), WR(r) AND THETA(r) 

       CALL MPI_ALLREDUCE (nirbin0,nirbin_t,nset*nset*nbin,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (wxbin0,wxbin_t,nset*nset*nbin,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (wybin0,wybin_t,nset*nset*nbin,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (wzbin0,wzbin_t,nset*nset*nbin,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (wrbin0,wrbin_t,nset*nset*nbin,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (urbin0,urbin_t,nset*nset*nbin,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (thetabin0,thetabin_t,nset*nset*nbin,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (wrtan0,wrtan_t,nset*nset*nbin,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

       nirbin   = nirbin + nirbin_t
       wxbin    = wxbin + wxbin_t
       wybin    = wybin + wybin_t
       wzbin    = wzbin + wzbin_t
       wrbin    = wrbin + wrbin_t
       urbin    = urbin + urbin_t
       thetabin = thetabin + thetabin_t
       wrtanbin = wrtanbin + wrtan_t

       IF ( mod (istep+1, 5).eq.0 ) THEN

          if ( id.eq.0 ) then

!----      write for all possible particle sets: 1-1, 1-2, 2-2 if there are just 2 sets.

           write (1011,125) nirbin (1,1,0:nbin-1)
           write (1012,125) nirbin (1,2,0:nbin-1)
           write (1022,125) nirbin (2,2,0:nbin-1)

           write (1110,126) wrbin (1,1,0:nbin-1)
           write (1120,126) wrbin (1,2,0:nbin-1)
           write (1220,126) wrbin (2,2,0:nbin-1)

           write (5110,126) wxbin (1,1,0:nbin-1)
           write (5120,126) wxbin (1,2,0:nbin-1)
           write (5220,126) wxbin (2,2,0:nbin-1)

           write (6110,126) wybin (1,1,0:nbin-1)
           write (6120,126) wybin (1,2,0:nbin-1)
           write (6220,126) wybin (2,2,0:nbin-1)

           write (7110,126) wzbin (1,1,0:nbin-1)
           write (7120,126) wzbin (1,2,0:nbin-1)
           write (7220,126) wzbin (2,2,0:nbin-1)

           write (3110,126) urbin (1,1,0:nbin-1)
           write (3120,126) urbin (1,2,0:nbin-1)
           write (3220,126) urbin (2,2,0:nbin-1)

           do ikk = 0, nbin-1 
            if ( nirbin (1,1,ikk) .ne. 0 ) thetabin (1,1,ikk) = thetabin (1,1,ikk) / real ( nirbin (1,1,ikk) )
            if ( nirbin (1,2,ikk) .ne. 0 ) thetabin (1,2,ikk) = thetabin (1,2,ikk) / real ( nirbin (1,2,ikk) )
            if ( nirbin (2,2,ikk) .ne. 0 ) thetabin (2,2,ikk) = thetabin (2,2,ikk) / real ( nirbin (2,2,ikk) )

!            if ( nirbin (1,1,ikk) .ne. 0 ) wrtanbin (1,1,ikk) = wrtanbin (1,1,ikk) / real ( nirbin (1,1,ikk) )
!            if ( nirbin (1,2,ikk) .ne. 0 ) wrtanbin (1,2,ikk) = wrtanbin (1,2,ikk) / real ( nirbin (1,2,ikk) )
!            if ( nirbin (2,2,ikk) .ne. 0 ) wrtanbin (2,2,ikk) = wrtanbin (2,2,ikk) / real ( nirbin (2,2,ikk) )
           enddo

           write (2110,126) thetabin (1,1,0:nbin-1)
           write (2120,126) thetabin (1,2,0:nbin-1)
           write (2220,126) thetabin (2,2,0:nbin-1)

           write (4110,126) wrtanbin (1,1,0:nbin-1)
           write (4120,126) wrtanbin (1,2,0:nbin-1)
           write (4220,126) wrtanbin (2,2,0:nbin-1)

          endif

          nirbin   = 0
          wrbin    = 0.
          wxbin    = 0.
          wybin    = 0.
          wzbin    = 0.
          urbin    = 0.
          thetabin = 0.
          wrtanbin = 0.

       ENDIF

125    format ( <nbin>(I6)      )
126    format ( <nbin>(1pe12.4) )

!-------------------------------------------------------------------
!  SAVE DATA FOR DYNAMIC COLLISION KERNEL

       CALL MPI_ALLREDUCE (icc1, icc1_t, 1, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
       CALL MPI_ALLREDUCE (icc2, icc2_t, 1, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
       CALL MPI_ALLREDUCE (icc3, icc3_t, 1, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)

       CALL MPI_ALLREDUCE (colltype,colltype_t,nset*nset,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)

       itotal = icc1_t + icc2_t + icc3_t

       if ( (id.eq.0) .and. (itotal.gt.0) ) then
!       if (id.eq.0) then
        write(50,302) time, icc1_t, icc2_t, icc3_t, colltype_t(1,1), colltype_t(1,2), colltype_t(2,2)
       endif

302    format ( 1pe15.5,7i10 )

!  SAVE DATA FOR KINEMATIC COLLISION KERNEL, RDF AND RELATIVE VELOCITY AT CONTACT
!--------------------------------------------------------------------------------

       CALL MPI_ALLREDUCE (nir, nir_t, nset*nset, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
       itotal = SUM ( nir_t )

       if ( (id.eq.0) .and. (itotal.gt.0) ) then
!       if (id.eq.0) then
        write(51,303) time, nir_t(1,1), nir_t(1,2), nir_t(2,2)
       endif

303    format ( 1pe15.5,3i10 )

!       DO n = 0, nproc-1

!        CALL MPI_BARRIER (MPI_COMM_WORLD, error)

!        IF ( (id.eq.n) .and. (ij.GT.0) ) THEN

!         do kj = 1, ij
!          index1 = int ( buffer(1,kj) )
!          index2 = int ( buffer(2,kj) )

!          indexfile = index1*100 + index2*10

!          write(indexfile,191) time, id, buffer(3:12,kj)
!         enddo

!        ENDIF

!        CALL MPI_BARRIER (MPI_COMM_WORLD, error)

!       ENDDO

       IF ( id.eq.0 ) THEN

         if ( ij.GT.0 ) then
          do kj = 1, ij

           index1 = int ( buffer(1,kj) )
           index2 = int ( buffer(2,kj) )

           indexfile = index1*100 + index2*10

           write(indexfile,191) time, buffer(3:12,kj)

          enddo
         endif

         DO n = 1, nproc-1

          CALL MPI_RECV (ijl,1,MPI_INTEGER,n,111,MPI_COMM_WORLD,status1,ierr)

          if ( ijl.gt.0 ) then
      
           ALLOCATE ( bufferrecv(12,ijl) )

           CALL MPI_RECV (bufferrecv,12*ijl,MPI_REAL8,n,112,MPI_COMM_WORLD,status1,ierr)

           do kj = 1, ijl

            index1 = int ( bufferrecv(1,kj) )
            index2 = int ( bufferrecv(2,kj) )

            indexfile  = index1*100 + index2*10

            write(indexfile,191) time, bufferrecv(3:12,kj)

           enddo

           DEALLOCATE ( bufferrecv )

          endif

         ENDDO

       ELSE

         CALL MPI_SEND (ij,1,MPI_INTEGER,0,111,MPI_COMM_WORLD,status1,ierr)
         if ( ij.gt.0 ) then
          CALL MPI_SEND (buffer,12*ij,MPI_REAL8,0,112,MPI_COMM_WORLD,status1,ierr)
         endif 

       ENDIF

191    format ( f10.6, 10(1x,f11.5) )
!191    format ( f10.6, i10, 10(1x,f11.5) )

!----------------------- SAVING DATA - END ----------------------------------

       DEALLOCATE ( nir, nir_t        )
       DEALLOCATE ( thetabin_t, thetabin0 )
       DEALLOCATE ( urbin_t, urbin0   )
       DEALLOCATE ( wrbin_t, wrbin0   )
       DEALLOCATE ( wxbin_t, wxbin0   )
       DEALLOCATE ( wybin_t, wybin0   )
       DEALLOCATE ( wzbin_t, wzbin0   )
       DEALLOCATE ( nirbin_t, nirbin0 )
       DEALLOCATE ( colltype          )
       DEALLOCATE ( colltype_t        )
       DEALLOCATE ( buffer            )
       DEALLOCATE ( ncollpart_t       )
       DEALLOCATE ( removepart_t      )

       RETURN
       END SUBROUTINE COLLISION_DET

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE FINDINDEX (rad1,rad2,index1,index2)

!     This subroutine gives the set number of the colliding or overlapping particles. 

      USE params
      INCLUDE 'mpif.h'

      diff1 = 1000000.
      diff2 = 1000000.

      index1 = 1
      index2 = 1

      DO il = 1, nset
         diff1n = abs ( radius(il) - rad1 )
         diff2n = abs ( radius(il) - rad2 )

         if ( diff1n.LE.diff1 ) then
           diff1  = diff1n
           index1 = il
         endif

         if ( diff2n.LE.diff2 ) then
           diff2  = diff2n
           index2 = il
         endif
      ENDDO

      if ( index2.lt.index1 ) then
         indextemp = index2
         index2    = index1
         index1    = indextemp
      endif

      RETURN 
      END SUBROUTINE FINDINDEX

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE MEAN_VEL (id)

!     Postprocessing: this subroutine computes the mean values...
!     of given variable (here particle velocity)

      USE params
      INCLUDE 'mpif.h'

      INTEGER                           :: error, ierr
      REAL                              :: vphoriz(nset), vpxyz(nset)
      REAL                              :: vphoriz_tot(nset), vpxyz_tot(nset)

      REAL, ALLOCATABLE, DIMENSION(:,:) :: vpmean, vpvar, vpabs
      REAL, ALLOCATABLE, DIMENSION(:,:) :: vpmean_tot, vpvar_tot, vpabs_tot

!      stime = mpi_wtime()

      ALLOCATE ( vpmean     (3,nset) )
      ALLOCATE ( vpvar      (3,nset) )
      ALLOCATE ( vpabs      (3,nset) )
      ALLOCATE ( vpmean_tot (3,nset) )
      ALLOCATE ( vpvar_tot  (3,nset) )
      ALLOCATE ( vpabs_tot  (3,nset) )

      vpmean(:,:) = 0.
      vpvar(:,:)  = 0.

      vpabs(:,:) = 0.
      vphoriz(:) = 0.
      vpxyz(:)   = 0.

      do ih = 1, nps
        dradius_min = 100.

        do ij = 1, nset
          dradius = abs ( yp(7,ih) - radius(ij) )

          if ( dradius.lt.dradius_min ) then
            dradius_min = dradius
            id_part     = ij
          endif

        enddo 

        vpmean(:,id_part) = vpmean(:,id_part) + yp(4:6,ih)
        vpvar(:,id_part)  = vpvar(:,id_part) + yp(4:6,ih)**2

        vpabs(:,id_part) = vpabs(:,id_part) + abs ( yp(4:6,ih) )
        vphoriz(id_part) = vphoriz(id_part) + sqrt ( yp(5,ih)**2 + yp(6,ih)**2 )
        vpxyz(id_part)   = vpxyz(id_part) + sqrt ( yp(4,ih)**2 + yp(5,ih)**2 + yp(6,ih)**2 )

      enddo

      vpmean(:,:) = vpmean(:,:) / npart * nset
      vpvar(:,:)  = vpvar(:,:) / npart * nset

      vpabs(:,:) = vpabs(:,:) / npart * nset
      vphoriz(:) = vphoriz(:) / npart * nset 
      vpxyz(:)   = vpxyz(:) / npart * nset 

!      etime      = mpi_wtime()
!      time_mv_cp = time_mv_cp + etime - stime
!      stime      = mpi_wtime()

      CALL MPI_ALLREDUCE (vpmean,vpmean_tot,nset*3,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE (vpvar,vpvar_tot,nset*3,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

      CALL MPI_ALLREDUCE (vpabs,vpabs_tot,nset*3,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE (vphoriz,vphoriz_tot,nset,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE (vpxyz,vpxyz_tot,nset,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

!      etime      = mpi_wtime()
!      time_mv_cm = time_mv_cm + etime - stime
!      stime      = mpi_wtime()

      if ( id.eq.0 ) then

       write(106,108) time, vpmean_tot(:,1), vpvar_tot(:,1),vpmean_tot(:,2), vpvar_tot(:,2)
       write(107,109) time, vpabs_tot(:,1), vphoriz_tot(1),vpxyz_tot(1), vpabs_tot(:,2), vphoriz_tot(2),vpxyz_tot(2) 

108    format ( 13(1pe12.4) )
109    format ( 11(1pe12.4) )

      endif 

!      etime      = mpi_wtime()
!      time_mv_sv = time_mv_sv + etime - stime

      DEALLOCATE ( vpmean, vpvar, vpabs, vpmean_tot, vpvar_tot, vpabs_tot )

      RETURN
      END SUBROUTINE MEAN_VEL

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE RELOCATE (id)
                                                                                                                                                                      
!     This subroutine relocates the collided particles provided "Non-overlap" simulation.
!     The list of collided particles is inside "removepart".
                                                                                                                                                                      
      USE params
      INCLUDE 'mpif.h'

      stime = mpi_wtime()

      DO ih = 1, ncoll
                                                                                                                                                                      
        ind = removepart(ih,0)

        r1  = ranff2(iseedr,ivr,iyr,idum2)
        r2  = ranff2(iseedr,ivr,iyr,idum2)
        r3  = ranff2(iseedr,ivr,iyr,idum2)
        
        yp(1,ind) = r1 * hx
        yp(2,ind) = r2 * hy + real(indy) * hy
        yp(3,ind) = r3 * hz + real(indz) * hz
            
      ENDDO

      etime           = mpi_wtime()
      time_colldet_cp = time_colldet_cp + etime - stime

      RETURN
      END SUBROUTINE RELOCATE 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE RELOCATE_WHOLE (id)

!     This subroutine relocates all particles in the domain
!     This scenario happens when simulating no-turbulence and gravitational cases only.
!     "Non-overlap" simulation.

      USE params
      INCLUDE 'mpif.h'

      stime = mpi_wtime()

      DO ih = 1, nps

!       r1  = ranff2(iseedr,ivr,iyr,idum2)
        r2  = ranff2(iseedr,ivr,iyr,idum2)
        r3  = ranff2(iseedr,ivr,iyr,idum2)

!        yp(1,ih) = r1 * hx
        yp(2,ih) = r2 * hy + real(indy) * hy
        yp(3,ih) = r3 * hz + real(indz) * hz

      ENDDO

      etime           = mpi_wtime()
      time_colldet_cp = time_colldet_cp + etime - stime

      RETURN
      END SUBROUTINE RELOCATE_WHOLE

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE OVERLAP_DET (id)

!     This subroutine detects particle overlaping. If particles are not allowed to overlap, this subroutine...
!     relocates the overlapping particles. 

      USE params 
      INCLUDE 'mpif.h'

      INTEGER           :: req(2), error, ierr
      INTEGER   	:: status1(MPI_STATUS_SIZE)
      INTEGER           :: status2(MPI_STATUS_SIZE,2)
   
      REAL              :: xb(3), xd(3), dnij, rcoll

      IF ( nps.eq.0 ) RETURN

      d_overl = 2.0 * radmax 
      imax    = 1 + int ( d_overl * real(ncdx) /  hx )
      jmax    = 1 + int ( d_overl * real(ncdy) /  hy )
      kmax    = 1 + int ( d_overl * real(ncdz) /  hz )

1222  CONTINUE

      CALL HEADLIST (id, d_overl)

!-----For each particle ih1 inside the pencil:

      noverl = 0

      DO ih1 = 1, nps

        ix = 1 + int ( ypn(1,ih1) / hx * real(ncdx) )
        iy = 1 + int ( ( ypn(2,ih1) - real(indy) * hy ) / hy * real(ncdy) )
        iz = 1 + int ( ( ypn(3,ih1) - real(indz) * hz ) / hz * real(ncdz) )
                                                                                                                                                                      
        if ( ix .gt. ncdx ) ix = ncdx
        if ( iy .gt. ncdy ) iy = ncdy
        if ( iz .gt. ncdz ) iz = ncdz

        DO i = -imax, imax
         DO j = -jmax, jmax
          DO k = -kmax, kmax 

!          restrict counting to neighboring cells.
           IF ( ( abs(i) .GT. 1 ) .OR. ( abs(j) .GT. 1 ) .OR. ( abs(k) .GT. 1 ) ) GOTO 201

!-------   Find the second particle ih2 in neighbourhood of particle ih1.

           iy2 = iy + j
           IF ( ( iy2 .LT. (1-jmax) ) .OR. ( iy2 .GT. (ncdy+jmax) ) ) GOTO 201

           iz2 = iz + k
           IF ( ( iz2 .LT. (1-kmax) ) .OR. ( iz2 .GT. (ncdz+kmax) ) ) GOTO 201

!          eliminate double counting.
           IF ( iz2 .LT. 1 ) GOTO 201
           IF ( ( iz2 .LE. ncdz ) .AND. ( iy2 .LT. 1 ) ) GOTO 201

           ix2 = ix + i
           sx  = 0.0

           if ( ix2.lt.1 ) then
            ix2 = ncdx + ix2
            sx  = -1.
           elseif ( ix2.gt.ncdx ) then
            ix2 = mod (ix2, ncdx) 
            sx  = 1.
           endif

!          2nd particle's index
           ih2 = head (ix2,iy2,iz2)
           IF ( ( ih2 .LE. ih1 ) .OR. ( ih2 .EQ. 0 ) ) GOTO 201

199        CONTINUE

!          location of the second droplet
           xb(1) = ypn(1,ih2) + sx*hx
           xb(2) = ypn(2,ih2)
           xb(3) = ypn(3,ih2)

           xd   = ypn(1:3,ih1) - xb
           dnij = sqrt ( xd(1)**2 + xd(2)**2 + xd(3)**2 )

           rcoll = ypn(7,ih1) + ypn(7,ih2)

           if  ( dnij .LE. rcoll ) then 

!           relocate the 1st particle (with index ih1)
            noverl = noverl + 1

            ar1 = ranff2(iseedr,ivr,iyr,idum2)
            ar2 = ranff2(iseedr,ivr,iyr,idum2)
            ar3 = ranff2(iseedr,ivr,iyr,idum2)
                                                                                                                                                                      
            yp(1,ih1) = ar1 * hx
            yp(2,ih1) = ar2 * hy + real(indy) * hy
            yp(3,ih1) = ar3 * hz + real(indz) * hz
        
           endif

           ih2 = list(ih2)
           IF ( (ih2.LE.ih1) .OR. (ih2.EQ.0) ) GOTO 201
           GOTO 199

201        CONTINUE

          ENDDO
         ENDDO
        ENDDO

       ENDDO

      CALL MPI_ALLREDUCE (noverl,noverl_tot,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
         
!     repeat until no particle overlaps inside the slab.
      IF ( noverl_tot .gt. 0 ) GOTO 1222

!     if ( id.eq.0 ) write(*,*) 'Finished overlapping detection'

      RETURN
      END SUBROUTINE OVERLAP_DET

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE HEADLIST (id, dist)

!     This subroutine collects the boundary droplets and builds head and list.
!     "dist" is the thickness of the halo zone around the rectangular subdomain.

      USE params
      INCLUDE 'mpif.h'

      REAL      :: dist, dist_r, dist_b

      INTEGER   :: ny_dist, nz_dist
      INTEGER   :: req(2), error, ierr
      INTEGER   :: status2(MPI_STATUS_SIZE,2)
      INTEGER   :: status1(MPI_STATUS_SIZE)

      REAL, ALLOCATABLE, DIMENSION (:,:,:)  :: temp_p, temp_b

      nsize = int ( 20.0 * nuniform ) 

      ALLOCATE ( temp_p  (18,nsize,8)  )
      ALLOCATE ( temp_b  (18,nsize,8)  )

      stime = mpi_wtime()

      ih_p  = 0
      ih_pt = 0

      ny_dist = 1 + int ( dist * real(ncdy) /  hy )
      nz_dist = 1 + int ( dist * real(ncdz) /  hz )

!     ypn is list of droplets (like yp) extended to include neighbor droplets.
      ypn (1:8, 1:nps)   = yp(1:8, 1:nps) 
      ypn (9:11, 1:nps)  = dmove(:, 1:nps) 
      ypn (12:14, 1:nps) = fvp(:, 0, 1:nps) 
      ypn (15, 1:nps)    = yp(8, 1:nps) 
      vpn (:, 1:nps)     = vp(:, 1:nps) 

      CALL MPI_BARRIER (MPI_COMM_WORLD, error)

!---- Determine the boundary droplets and gather required data

      DO ih = 1, nps

        d_e = abs ( real(indy+1) * hy - yp(2,ih) ) 
        d_w = abs ( real(indy) * hy   - yp(2,ih) ) 
        d_n = abs ( real(indz+1) * hz - yp(3,ih) ) 
        d_s = abs ( real(indz) * hz   - yp(3,ih) ) 
                                                                                                                                                                     
        if ( d_e .le. dist ) then
                                                                                                                                 
           nid       = 1
           ih_p(nid) = ih_p(nid) + 1
           pindex (ih_p(nid), nid) = ih  
                         
           temp_p (1:8, ih_p(nid), nid)   = yp(1:8,ih)
           temp_p (9:11, ih_p(nid), nid)  = dmove(:,ih)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,0,ih)
           temp_p (15, ih_p(nid), nid)    = real (ih)
           temp_p (16:18, ih_p(nid), nid) = vp(1:3,ih) 
            
           if ( (indy+1) .eq. nprocY ) temp_p (2, ih_p(nid), nid) = yp(2,ih) - pi2 

        endif

        if ( sqrt(d_e**2+d_n**2) .le. dist ) then
                                                                                                                                                                     
           nid       = 2 
           ih_p(nid) = ih_p(nid) + 1
           pindex (ih_p(nid), nid) = ih  
                                                                                                                                                                     
           temp_p (1:8, ih_p(nid), nid)   = yp(1:8,ih)
           temp_p (9:11, ih_p(nid), nid)  = dmove(:,ih)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,0,ih)
           temp_p (15, ih_p(nid), nid)    = real (ih)
           temp_p (16:18, ih_p(nid), nid) = vp(1:3,ih) 

           if ((indy+1) .eq. nprocY) temp_p (2, ih_p(nid), nid) = yp(2,ih) - pi2 
           if ((indz+1) .eq. nprocZ) temp_p (3, ih_p(nid), nid) = yp(3,ih) - pi2 

        endif

        if ( d_n .le. dist ) then
                                                                                                                                                                     
           nid       = 3 
           ih_p(nid) = ih_p(nid) + 1
           pindex (ih_p(nid), nid) = ih  
                                                                                                                                                                     
           temp_p (1:8, ih_p(nid), nid)   = yp(1:8,ih)
           temp_p (9:11, ih_p(nid), nid)  = dmove(:,ih)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,0,ih)
           temp_p (15, ih_p(nid), nid)    = real (ih)
           temp_p (16:18, ih_p(nid), nid) = vp(1:3,ih) 

           if ( (indz+1) .eq. nprocZ ) temp_p (3, ih_p(nid), nid) = yp(3,ih) - pi2 

        endif

        if ( sqrt(d_w**2+d_n**2) .le. dist ) then
                                                                                                                                                                      
           nid       = 4 
           ih_p(nid) = ih_p(nid) + 1
           pindex (ih_p(nid), nid) = ih  
                                                                                                                                                                      
           temp_p (1:8, ih_p(nid), nid)   = yp(1:8,ih)
           temp_p (9:11, ih_p(nid), nid)  = dmove(:,ih)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,0,ih)
           temp_p (15, ih_p(nid), nid)    = real (ih)
           temp_p (16:18, ih_p(nid), nid) = vp(1:3,ih) 

           if (indy .eq. 0)          temp_p (2, ih_p(nid), nid) = yp(2,ih) + pi2 
           if ((indz+1) .eq. nprocZ) temp_p (3, ih_p(nid), nid) = yp(3,ih) - pi2 

        endif

        if ( d_w .le. dist ) then
                                                                                                                                                                      
           nid       = 5 
           ih_p(nid) = ih_p(nid) + 1
           pindex (ih_p(nid), nid) = ih  
                                                                                                                                                                      
           temp_p (1:8, ih_p(nid), nid)   = yp(1:8,ih)
           temp_p (9:11, ih_p(nid), nid)  = dmove(:,ih)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,0,ih)
           temp_p (15, ih_p(nid), nid)    = real (ih)
           temp_p (16:18, ih_p(nid), nid) = vp(1:3,ih) 

           if ( indy .eq. 0 ) temp_p (2, ih_p(nid), nid) = yp(2,ih) + pi2 

        endif

        if ( sqrt(d_w**2+d_s**2) .le. dist ) then
                                                                                                                                                                      
           nid       = 6 
           ih_p(nid) = ih_p(nid) + 1
           pindex (ih_p(nid), nid) = ih  
                                                                                                                                                                      
           temp_p (1:8, ih_p(nid), nid)   = yp(1:8,ih)
           temp_p (9:11, ih_p(nid), nid)  = dmove(:,ih)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,0,ih)
           temp_p (15, ih_p(nid), nid)    = real (ih)
           temp_p (16:18, ih_p(nid), nid) = vp(1:3,ih) 

           if (indy .eq. 0) temp_p (2, ih_p(nid), nid) = yp(2,ih) + pi2 
           if (indz .eq. 0) temp_p (3, ih_p(nid), nid) = yp(3,ih) + pi2 

        endif

        if ( d_s .le. dist ) then
                                                                                                                                                                      
           nid       = 7 
           ih_p(nid) = ih_p(nid) + 1
           pindex (ih_p(nid), nid) = ih  
                                                                                                                                                                      
           temp_p (1:8, ih_p(nid), nid)   = yp(1:8,ih)
           temp_p (9:11, ih_p(nid), nid)  = dmove(:,ih)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,0,ih)
           temp_p (15, ih_p(nid), nid)    = real (ih)
           temp_p (16:18, ih_p(nid), nid) = vp(1:3,ih) 

           if ( indz .eq. 0 ) temp_p (3, ih_p(nid), nid) = yp(3,ih) + pi2 

        endif

        if ( sqrt(d_e**2+d_s**2) .le. dist ) then
                                                                                                                                                                      
           nid       = 8 
           ih_p(nid) = ih_p(nid) + 1
           pindex (ih_p(nid), nid) = ih  
                                                                                                                                                                      
           temp_p (1:8, ih_p(nid), nid)   = yp(1:8,ih)
           temp_p (9:11, ih_p(nid), nid)  = dmove(:,ih)
           temp_p (12:14, ih_p(nid), nid) = fvp(:,0,ih)
           temp_p (15, ih_p(nid), nid)    = real (ih)
           temp_p (16:18, ih_p(nid), nid) = vp(1:3,ih) 

           if ((indy+1) .eq. nprocY) temp_p (2, ih_p(nid), nid) = yp(2,ih) - pi2 
           if (indz .eq. 0)          temp_p (3, ih_p(nid), nid) = yp(3,ih) + pi2 

        endif

      END DO

      etime           = mpi_wtime()
      time_colldet_cp = time_colldet_cp + etime - stime
      stime           = mpi_wtime()

!----- communicate the buffer zone droplets

       DO jid = 1, 8

        ids = nei_id ( jid )
        idr = nei_id ( mod (jid+3, 8) + 1 )

        CALL MPI_SENDRECV (ih_p(jid),1,MPI_INTEGER,ids,212,ih_pt(jid),1,MPI_INTEGER,idr,212,&
                       MPI_COMM_WORLD,status1,ierr)

        if ( ih_p(jid) .gt. 0 ) then
         CALL MPI_SEND (temp_p(:,1:ih_p(jid),jid),18*ih_p(jid),MPI_REAL8,ids,213,MPI_COMM_WORLD,ierr)
        endif

        if ( ih_pt(jid) .gt. 0 ) then
         CALL MPI_RECV (temp_b(:,1:ih_pt(jid),jid),18*ih_pt(jid),MPI_REAL8,idr,213,MPI_COMM_WORLD,status1,ierr)
        endif

!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAITALL (2,req,status2,ierr)
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).eq.0 ) CALL MPI_WAIT (req(1),status1,ierr)
!        if ( ih_p(jid).eq.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAIT (req(2),status1,ierr)

        if ( ih_pt(jid) .gt. 0 ) then

          nps1 = nps + SUM ( ih_pt ( 1 : (jid-1) ) )
          nps2 = nps + SUM ( ih_pt ( 1 : jid ) ) 

          ypn (1:15, (nps1+1):nps2) = temp_b(1:15,1:ih_pt(jid),jid)
          vpn (1:3, (nps1+1):nps2)  = temp_b(16:18,1:ih_pt(jid),jid)

         endif

        ENDDO

       if ( SUM ( ih_pt(1:8) ) .eq. 0 ) then
        nps1 = nps
        nps2 = nps
       endif

       etime           = mpi_wtime()
       time_colldet_cm = time_colldet_cm + etime - stime
       stime           = mpi_wtime()

!-----Find the heads and fill the list 

       list = 0
       head = 0

       DO ih = 1, nps2
                                                                                                                                                                    
        ix = 1 + int ( ypn(1,ih) / hx * real(ncdx) )

        dist_r = ypn(2,ih) - real(indy) * hy
        if ( dist_r .gt. 0 ) then
         iy = 1 + int ( dist_r / hy * real(ncdy) )
        else
         iy = - int ( abs (dist_r) / hy * real(ncdy) )
        endif

        dist_b = ypn(3,ih) - real(indz) * hz 
        if ( dist_b .gt. 0 ) then
         iz = 1 + int ( dist_b / hz * real(ncdz) )
        else
         iz = - int ( abs (dist_b) / hz * real(ncdz) )
        endif

        if ( ix .gt. ncdx ) ix = ncdx

        if ( iy .gt. (ncdy+ny_max) ) iy = ncdy + ny_max
        if ( iy .lt. (1-ny_max) )    iy = 1 - ny_max

        if ( iz .gt. (ncdz+nz_max) ) iz = ncdz + nz_max
        if ( iz .lt. (1-nz_max) )    iz = 1 - nz_max

        list (ih)       = head (ix,iy,iz)
        head (ix,iy,iz) = ih

       END DO

       CALL MPI_BARRIER (MPI_COMM_WORLD, error)

       DEALLOCATE ( temp_p, temp_b )

       etime           = mpi_wtime()
       time_colldet_cp = time_colldet_cp + etime - stime

       RETURN
       END SUBROUTINE HEADLIST

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE AFUN (x, id)

!     For any x, this subroutine computes Ax in which A is the matrix obtained by...
!     reforming linear HDI system into the form Ax = b. This is needed to build the ...
!     Krylov subspace. 
!     This subroutine is called from HDI_GMRES subroutine.
!     Written by Hossein Parishani, University of Delaware.
!     Questions, please contact hparish@udel.edu

      USE params 
      INCLUDE 'mpif.h'

      REAL      :: xb(3), xd(3), relv(3) ,relvb(3)
      INTEGER   :: req(2), hdipart, error, ierr
      INTEGER   :: status1(MPI_STATUS_SIZE), status2(MPI_STATUS_SIZE,2)

      REAL,    DIMENSION (3*nps,1)               :: x

      REAL,    ALLOCATABLE, DIMENSION (:,:)      :: Ax, xreshaped
      REAL,    ALLOCATABLE, DIMENSION (:,:)      :: RHSof6, bof6 
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)    :: temp_p, temp_b 

      stime   = mpi_wtime()

      hdipart = int ( 30 * nuniform )

      ALLOCATE ( xreshaped   (3,nps)   )
      ALLOCATE ( Ax          (3,nps)   )
      ALLOCATE ( RHSof6      (3,nps)   )
      ALLOCATE ( bof6        (3,nps)   )
      ALLOCATE ( temp_p  (6,hdipart,8) )
      ALLOCATE ( temp_b  (6,hdipart,8) )

      xreshaped(1, 1:nps) = x(1:nps,1) 
      xreshaped(2, 1:nps) = x(nps+1:2*nps,1) 
      xreshaped(3, 1:nps) = x(2*nps+1:3*nps,1) 

      ypn (16:18, 1:nps) = xreshaped(1:3, 1:nps) 
      ypn (19:21, 1:nps) = yp(4:6, 1:nps) - vp(1:3, 1:nps)  

      DO jid = 1, 8 
       DO jik = 1, ih_p(jid) 

        partind              = pindex (jik, jid)
        temp_p (1:3,jik,jid) = xreshaped(1:3, partind)  
        temp_p (4:6,jik,jid) = yp(4:6, partind) - vp(1:3, partind)  

       ENDDO
      ENDDO

      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime
      stime       = mpi_wtime()

!-----Communicate the particles close to right & left subdomain boundaries

      DO jid = 1, 8
                                                                                                                                                                    
        ids = nei_id ( jid )
        idr = nei_id ( mod (jid+3, 8) + 1 )
                                                                                                                                                                    
        if ( ih_p(jid) .gt. 0 ) then
         CALL MPI_SEND (temp_p(:,1:ih_p(jid),jid),6*ih_p(jid),MPI_REAL8,ids,213,MPI_COMM_WORLD,ierr)
        endif
                                                                                                                                                                    
        if ( ih_pt(jid) .gt. 0 ) then
         CALL MPI_RECV (temp_b(:,1:ih_pt(jid),jid),6*ih_pt(jid),MPI_REAL8,idr,213,MPI_COMM_WORLD,status1,ierr)
        endif
                                                                                                                                                                    
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAITALL (2,req,status2,ierr)
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).eq.0 ) CALL MPI_WAIT (req(1),status1,ierr)
!        if ( ih_p(jid).eq.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAIT (req(2),status1,ierr)
                                                                                                                                                                    
        if ( ih_pt(jid) .gt. 0 ) then
         
          nps1 = nps + SUM ( ih_pt ( 1 : (jid-1) ) )
          nps2 = nps + SUM ( ih_pt ( 1 : jid ) )
         
          ypn (16:21, (nps1+1):nps2) = temp_b(1:6,1:ih_pt(jid),jid)
         
        endif
         
      ENDDO

      etime       = mpi_wtime()
      time_hdi_cm = time_hdi_cm + etime - stime
      stime       = mpi_wtime()

      RHSof6 = 0.0
      bof6   = 0.0

      npstest = nps + SUM ( ih_pt (1:8) )

!-----afun for particles inside the current task, including boundary droplets.
!-----when completed, update RHSof6.

      DO ih1 = 1, nps

!       find the cells in the current slab which are falling in a sphere with radius d_hdi

        ix = 1 + int ( ypn(1,ih1) / hx*real(ncdx) )
        iy = 1 + int ( (ypn(2,ih1) - real(indy)*hy) / hy*real(ncdy) )
        iz = 1 + int ( (ypn(3,ih1) - real(indz)*hz) / hz*real(ncdz) )
       
        if ( ix.gt.ncdx ) ix = ncdx
        if ( iy.gt.ncdy ) iy = ncdy
        if ( iz.gt.ncdz ) iz = ncdz
       
        ix_max = ix + nx_max
        ix_min = ix - nx_max
        iy_max = iy + ny_max
        iy_min = iy - ny_max
        iz_max = iz + nz_max
        iz_min = iz - nz_max
       
        if ( iy_max .gt. (ncdy+ny_max) ) iy_max = ncdy + ny_max
        if ( iy_min .lt. (1-ny_max) )    iy_min = 1 - ny_max
       
        if ( iz_max .gt. (ncdz+nz_max) ) iz_max = ncdz + nz_max
        if ( iz_min .lt. (1-nz_max) )    iz_min = 1 - nz_max

!        dist_x = ypn(1,ih1) + d_hdi
!        if ( dist_x .gt. 0 ) then
!         ix_max = 1 + int ( dist_x / hx * real(ncdx) )
!        else
!         ix_max = - int ( abs (dist_x) / hx * real(ncdx) )
!        endif

!        dist_x = ypn(1,ih1) - d_hdi
!        if ( dist_x .gt. 0 ) then
!         ix_min = 1 + int ( dist_x / hx * real(ncdx) )
!        else
!         ix_min = - int ( abs (dist_x) / hx * real(ncdx) )
!        endif

!        dist_r = ypn(2,ih1) + d_hdi - real(indy) * hy
!        if ( dist_r .gt. 0 ) then
!         iy_max = 1 + int ( dist_r / hy * real(ncdy) )
!        else
!         iy_max = - int ( abs (dist_r) / hy * real(ncdy) )
!        endif

!        dist_r = ypn(2,ih1) - d_hdi - real(indy) * hy
!        if ( dist_r .gt. 0 ) then
!         iy_min = 1 + int ( dist_r / hy * real(ncdy) )
!        else
!         iy_min = - int ( abs (dist_r) / hy * real(ncdy) )
!        endif

!        if ( iy_max .gt. (ncdy+ny_max) ) iy_max = ncdy + ny_max
!        if ( iy_max .lt. (1-ny_max) )    iy_max = 1 - ny_max
!        if ( iy_min .gt. (ncdy+ny_max) ) iy_min = ncdy + ny_max
!        if ( iy_min .lt. (1-ny_max) )    iy_min = 1 - ny_max

!        dist_b = ypn(3,ih1) + d_hdi - real(indz) * hz
!        if ( dist_b .gt. 0 ) then
!         iz_max = 1 + int ( dist_b / hz * real(ncdz) )
!        else
!         iz_max = - int ( abs (dist_b) / hz * real(ncdz) )
!        endif

!        dist_b = ypn(3,ih1) - d_hdi - real(indz) * hz
!        if ( dist_b .gt. 0 ) then
!         iz_min = 1 + int ( dist_b / hz * real(ncdz) )
!        else
!         iz_min = - int ( abs (dist_b) / hz * real(ncdz) )
!        endif

!        if ( iz_max .gt. (ncdz+nz_max) ) iz_max = ncdz + nz_max
!        if ( iz_max .lt. (1-nz_max) )    iz_max = 1 - nz_max
!        if ( iz_min .gt. (ncdz+nz_max) ) iz_min = ncdz + nz_max
!        if ( iz_min .lt. (1-nz_max) )    iz_min = 1 - nz_max

        do i = ix_min, ix_max
         do j = iy_min, iy_max
          do k = iz_min, iz_max

           ii = i

           sx = 0.0
           if ( i .gt. ncdx ) then
            ii = mod (i, ncdx)
            sx = 1.0
           elseif ( i .lt. 1 ) then
            ii = ncdx + i 
            sx = -1.0
           endif

           ih2 = head (ii,j,k)
           IF ( ih2 .eq. 0 ) GOTO 2201

1199       CONTINUE

!          for all particles with greater indices          
           IF ( ih2 .le. ih1 ) GOTO 2201
  
!           location of the second droplet
            xb(1) = ypn(1,ih2) + sx*hx
            xb(2) = ypn(2,ih2)
            xb(3) = ypn(3,ih2)

!           distance between the 2 droplets
            xd(:) = ypn(1:3,ih1) - xb(:)
            dnij  = sqrt ( xd(1)**2 + xd(2)**2 + xd(3)**2 )

            if ( dnij .lt. (ypn(7,ih1)+ypn(7,ih2)) ) then
             write(*,*)'inside AFUN: overlapping is not allowed'
             write(*,*)'nps,npstest,ih1,ih2,dnij,ypn(1:8,ih1),ypn(1:8,ih2) =',nps,npstest,ih1,ih2,dnij,ypn(1:8,ih1),ypn(1:8,ih2)   
             CALL mpi_finalize(ierr)
             STOP
            endif

!           for particle with index ih1          

            hdirad2 = HDI_trunc * ypn(7,ih2) 
!            hdirad2 = HDI_trunc * radmax 

            if ( dnij .lt. hdirad2 ) then

             fact1 = ypn(7,ih2) / dnij
             fact2 = fact1**3
             fact1 = 0.75 * fact1
 
             a1 = fact1 - 0.75 * fact2
             a2 = fact1 + 0.25 * fact2

             relv(:)       = ypn(16:18,ih2) 
             aloop         = a1 * DOT_PRODUCT ( relv,xd ) / dnij**2
             RHSof6(:,ih1) = RHSof6(:,ih1) + aloop*xd(:) + a2*relv(:)

             relvb(:)    = ypn(19:21,ih2)
             aloopb      = a1 * DOT_PRODUCT ( relvb,xd ) / dnij**2
             bof6(:,ih1) = bof6(:,ih1) + aloopb*xd(:) + a2*relvb(:)

            endif

!           for particle with index ih2
                                                                                                                                                       
            IF ( ih2 .le. nps ) THEN

             hdirad2 = HDI_trunc * ypn(7,ih1)
!             hdirad2 = HDI_trunc * radmax 
                                                                                                                                                       
             if ( dnij .lt. hdirad2 ) then
                                                                                                                                                     
              xd    = - xd 
 
              fact1 = ypn(7,ih1) / dnij
              fact2 = fact1**3
              fact1 = 0.75 * fact1
                                                                                                                                                       
              a1 = fact1 - 0.75 * fact2
              a2 = fact1 + 0.25 * fact2
                                                                                                                                                       
              relv(:)       = ypn(16:18,ih1)
              aloop         = a1 * DOT_PRODUCT ( relv,xd ) / dnij**2
              RHSof6(:,ih2) = RHSof6(:,ih2) + aloop*xd(:) + a2*relv(:)
                                                                                                                                                       
              relvb(:)    = ypn(19:21,ih1)
              aloopb      = a1 * DOT_PRODUCT ( relvb,xd ) / dnij**2
              bof6(:,ih2) = bof6(:,ih2) + aloopb*xd(:) + a2*relvb(:)
                                                                                                                                                       
             endif

            ENDIF

           ih2 = list(ih2)
           IF ( ih2 .le. ih1 ) GOTO 2201

           GOTO 1199

2201       CONTINUE

          enddo 
         enddo
        enddo

      ENDDO

      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime
      stime       = mpi_wtime()

      Ax(:,1:nps) = RHSof6(:,1:nps) + xreshaped(:,1:nps)

      Axreshaped(1:nps,1)         = Ax(1, 1:nps)
      Axreshaped(nps+1:2*nps,1)   = Ax(2, 1:nps)
      Axreshaped(2*nps+1:3*nps,1) = Ax(3, 1:nps)

      breshaped(1:nps,1)          = bof6 (1, 1:nps)
      breshaped(nps+1:2*nps,1)    = bof6 (2, 1:nps)
      breshaped(2*nps+1:3*nps,1)  = bof6 (3, 1:nps)

      DEALLOCATE ( RHSof6     )
      DEALLOCATE ( temp_p     )
      DEALLOCATE ( temp_b     )
      DEALLOCATE ( bof6       )
      DEALLOCATE ( Ax         )
      DEALLOCATE ( xreshaped  )

      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime

      RETURN
      END SUBROUTINE AFUN

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE EUC_NORM (matx, dim1, dim2, id, iglob)

!     This subroutine computes Euclidean norm of matrix matx and updates NORM.
!     dim1 & dim2 are dimensions of matx.
!     If iglob=1, NORM will be global norm.
!     This subroutine is called from HDI_GMRES subroutine.

      USE params 
      INCLUDE 'mpif.h'

      INTEGER                        ::  ierr, dim1, dim2, i, j, iglob
      REAL                           ::  x_c, x_sum, x_abs
      REAL, DIMENSION  (dim1,dim2)   ::  matx 
                                                                                                      
      x_c   = 0.0
      x_abs = 0.0
      NORM  = 0.0 
                                                                                                                                                             
      DO j = 1, dim2 
       DO i = 1, dim1 
                                                                                                                                                             
        x_abs = matx(i,j)**2
        x_c   = x_c + x_abs 
                                                                                                                                                             
       ENDDO
      ENDDO

      NORM = sqrt ( x_c )

      if (iglob .eq. 1) then
       CALL MPI_ALLREDUCE (x_c,x_sum,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       NORM = sqrt ( x_sum )
      endif

      RETURN
      END SUBROUTINE EUC_NORM 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE HDI_GMRES (id)

!     This subroutine solves the HDI system of linear equations via GMRES and...
!     stores the disturbance velocity field in global matrix "pertvel".
!     Written by Hossein Parishani, University of Delaware.
!     Questions, please contact hparish@udel.edu

      USE params 
      INCLUDE 'mpif.h'

      REAL      :: norm_rhs, norm_b, aa, bb, cc, dd, ss, maxerr
      REAL      :: q1(2,1), q2(2,1), res, dotprod, globdot
      INTEGER   :: iter, maxiter, psize

      REAL,    ALLOCATABLE, DIMENSION (:)        :: v1, v2
      REAL,    ALLOCATABLE, DIMENSION (:,:)      :: bn, H, Rn_h, Qn, v
      REAL,    ALLOCATABLE, DIMENSION (:,:)      :: x0, x, Q, y, resvector

      stime = mpi_wtime()

      maxerr   = 1.0e-6
      maxiter  = 50 
      psize    = 3 * nps

      ALLOCATE ( bn       (maxiter+1,1)         )
      ALLOCATE ( H        (maxiter+1,maxiter)   )
      ALLOCATE ( Rn_h     (maxiter+1,maxiter+1) )
      ALLOCATE ( Qn       (maxiter+1,maxiter+1) )

      ALLOCATE ( x0  (psize,1)   )
      ALLOCATE ( x   (psize,1)   )
      ALLOCATE ( y   (maxiter,1) )
      ALLOCATE ( v1  (maxiter+1) )
      ALLOCATE ( v2  (maxiter+1) )
      ALLOCATE ( v   (psize,1)   )

      ALLOCATE ( Q   (psize,maxiter+1)   )
      ALLOCATE ( resvector (maxiter+1,1) )

      H    = 0.0
      Rn_h = 0.0
      Qn   = 0.0
      Q    = 0.0
      bn   = 0.0
      v    = 0.0
      y    = 0.0

!     initial guess, use pertvel from previous timestep

      x0(1:nps,1)         = pertvel(1,1:nps) 
      x0(nps+1:2*nps,1)   = pertvel(2,1:nps) 
      x0(2*nps+1:psize,1) = pertvel(3,1:nps) 

      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime

!     CALL AFUN (x0, id)
      CALL AFUN_PRECOMP (x0, id, 1)

      stime = mpi_wtime()

      Q (1:psize,1) = breshaped (1:psize,1) - Axreshaped (1:psize,1) 

      CALL EUC_NORM ( breshaped(1:psize,1), psize, 1, id, 1 ) 
      norm_b = NORM

      CALL EUC_NORM ( Q(1:psize,1), psize, 1, id, 1 ) 
      norm_rhs = NORM

!      if (id .eq. 10) write (*,*) 'norm_rhs =', norm_rhs

      IF ( norm_rhs .lt. 1e-5 ) then
       RETURN
      ENDIF

!     build basis for Krylov subspace, first vector
      Q(:,1)  = Q(:,1) / norm_rhs 

      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime

      DO iter = 1, maxiter 

!     CALL AFUN ( Q(1:psize,iter), id )
      CALL AFUN_PRECOMP ( Q(1:psize,iter), id, 2 )

       stime = mpi_wtime()

       v (1:psize,1) = Axreshaped (1:psize,1)

!      build the Hessenberg matrix H
       do j=1, iter

        dotprod = DOT_PRODUCT ( Q(1:psize,j), v(1:psize,1) )
        CALL MPI_ALLREDUCE (dotprod,globdot,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

        H(j,iter)    = globdot 
        v(1:psize,1) = v(1:psize,1) - H(j,iter) * Q(1:psize,j)

       enddo

       CALL EUC_NORM ( v(1:psize,1), psize, 1, id, 1 ) 
       H(iter+1,iter) = NORM

!      build basis for Krylov subspace, next vector
       Q(:,iter+1) = v(:,1) / H(iter+1,iter)

!================== SOLVING LEAST SQUARES PROBLEM =========================

       if ( iter .eq. 1 ) then

        q1(1:2,1) = H(1:2,1) / sqrt ( H(1,1)**2 + H(2,1)**2 )
        Qn(1:2,1) = q1(1:2,1)
        q2 (1,1)  = q1(2,1)
        q2 (2,1)  = -q1(1,1)
        Qn(1:2,2) = q2(1:2,1)

        Rn_h (1,1) = sqrt ( H(1,1)**2 + H(2,1)**2 ) 
        bn (1,1)   = q1(1,1) * norm_rhs
        y(1,1)     = bn (1,1) / Rn_h (1,1) 
                      
       else
          
!       Update Rn, Given's rotation
        Rn_h (1:iter, iter) = MATMUL ( Qn(1:iter, 1:iter), H(1:iter, iter) )
        aa = Rn_h (iter, iter)
        bb = H (iter+1, iter)
        dd = sqrt ( aa**2 + bb**2 )
        cc = aa / dd
        ss = bb / dd
        Rn_h (iter, iter) = dd
                                                                                                                                                            
!       Update Qn
        v1(1:iter) = cc * Qn (iter, 1:iter)
        v1(iter+1) = ss
        v2(1:iter) = -ss * Qn (iter, 1:iter)
        v2(iter+1) = cc

        Qn (iter, 1:iter+1)   = v1(1:iter+1)
        Qn (iter+1, 1:iter+1) = v2(1:iter+1)
                                                                                                                                                             
!       Update bnn
        bn(iter-1:iter,1) = Qn (iter-1:iter, 1) * norm_rhs
                                                                                                                                                            
!       Finding y, back substitution
        y(iter,1) = bn(iter,1) / Rn_h(iter, iter)
 
        do i = iter-1, 1, -1
         y(i,1) =  ( bn(i,1) - DOT_PRODUCT ( y(i+1:iter,1), Rn_h(i,i+1:iter) ) ) / Rn_h(i,i)
        enddo

       endif

!================== LEAST SQUARES PROBLEM - END =========================

!      Find the solution x = x0 + Qn * yn
       x(1:psize, 1) = x0(1:psize, 1) + MATMUL ( Q(1:psize,1:iter), y(1:iter,1) )

!      Update pertvel
       pertvel(1,1:nps) = x(1:nps,1) 
       pertvel(2,1:nps) = x(nps+1:2*nps,1) 
       pertvel(3,1:nps) = x(2*nps+1:psize,1)  

!      Compute residual
       resvector(1:iter+1, 1) = MATMUL ( H(1:iter+1, 1:iter), y(1:iter,1) )
       resvector(1, 1)        = resvector(1, 1) - norm_rhs

       CALL EUC_NORM ( resvector(1:iter+1, 1), iter+1, 1, id, 0 ) 
       res = NORM / norm_b

       if ( id .eq. 0 ) write (*,*) 'istep, iter, residual from GMRES code = ', istep, iter, res

!       CALL RESIDUAL (id,iter)

       if ( res .LT. maxerr) GOTO 973  
       if ( iter .EQ. (maxiter / 2) ) maxerr = maxerr * 10
       if ( iter .GE. (maxiter - 2) ) maxerr = maxerr * 10

       etime       = mpi_wtime()
       time_hdi_cp = time_hdi_cp + etime - stime

      ENDDO

      stime = mpi_wtime()

      if ( (res .GT. maxerr) .AND. (id .EQ. 0) ) then
       write (*,*) 'WARNING: HDI did not converge completely, iter =', iter
      endif

973   CONTINUE

      DEALLOCATE ( bn   )
      DEALLOCATE ( H    )
      DEALLOCATE ( Rn_h )
      DEALLOCATE ( Qn   )
      DEALLOCATE ( x0   )
      DEALLOCATE ( x    )
      DEALLOCATE ( y    )
      DEALLOCATE ( v1   )
      DEALLOCATE ( v2   )
      DEALLOCATE ( v    )
      DEALLOCATE ( Q    )
      DEALLOCATE ( resvector )
      DEALLOCATE ( nei_num   )
      DEALLOCATE ( nei_data  )

      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime

      RETURN
      END SUBROUTINE HDI_GMRES

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE RESIDUAL (id, iter)

!     This subroutine computes norm of the residulal for HDI solver. 

      USE params 
      INCLUDE 'mpif.h'

      REAL      :: xb_res(3), xd_res(3), relv_res(3), relvb(3), RHSminusLHS(3)
      REAL      :: error_norm

      INTEGER   :: req(2), ipart
      INTEGER   :: status2(MPI_STATUS_SIZE,2)
      INTEGER   :: status1(MPI_STATUS_SIZE)

      REAL,    ALLOCATABLE, DIMENSION (:,:)      :: RHSof6, RHSof6b 
      REAL,    ALLOCATABLE, DIMENSION (:,:)      :: ypres 
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)    :: temp_p, temp_b 

      maxpart = int ( 30 * nuniform )

      ALLOCATE ( temp_p  (6,maxpart,8) )
      ALLOCATE ( temp_b  (6,maxpart,8) )

      ALLOCATE ( ypres    (6,2*maxpart)   )
      ALLOCATE ( RHSof6   (3,nps)   )
      ALLOCATE ( RHSof6b  (3,nps)   )

      n_inter  = 0 

      ypres (1:3, 1:nps) = yp(4:6,1:nps) - vp(:,1:nps) - pertvel(:,1:nps)  
      ypres (4:6, 1:nps) = pertvel(:,1:nps)
                                                                                                                                                                     
      DO jid = 1, 8
       DO jik = 1, ih_p(jid)
                                                                                                                                                                     
        partind              = pindex (jik, jid)
        temp_p (1:3,jik,jid) = yp(4:6,partind) - vp(:,partind) - pertvel(:,partind) 
        temp_p (4:6,jik,jid) = pertvel(:,partind) 

       ENDDO
      ENDDO

!-----Communicate the particles close to subdomain boundaries

      DO jid = 1, 8
                                                                                                                                                                     
        ids = nei_id ( jid )
        idr = nei_id ( mod (jid+3, 8) + 1 )
                                                                                                                                                                     
        if ( ih_p(jid) .gt. 0 ) then
         CALL MPI_SEND (temp_p(:,1:ih_p(jid),jid),6*ih_p(jid),MPI_REAL8,ids,213,MPI_COMM_WORLD,ierr)
        endif
                                                                                                                                                                     
        if ( ih_pt(jid) .gt. 0 ) then
         CALL MPI_RECV (temp_b(:,1:ih_pt(jid),jid),6*ih_pt(jid),MPI_REAL8,idr,213,MPI_COMM_WORLD,status1,ierr)
        endif
                                                                                                                                                                     
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAITALL (2,req,status2,ierr)
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).eq.0 ) CALL MPI_WAIT (req(1),status1,ierr)
!        if ( ih_p(jid).eq.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAIT (req(2),status1,ierr)
                                                                                                                                                                     
        if ( ih_pt(jid) .gt. 0 ) then
                                                                                                                                                                     
          nps1 = nps + SUM ( ih_pt ( 1 : (jid-1) ) )
          nps2 = nps + SUM ( ih_pt ( 1 : jid ) )
                                                                                                                                                                     
          ypres (1:6, (nps1+1):nps2) = temp_b(1:6,1:ih_pt(jid),jid)
                                                                                                                                                                     
        endif
                                                                                                                                                                     
      ENDDO

      RHSof6  = 0.0
      RHSof6b = 0.0

!-----RESIDUAL for particles inside the current task. when completed, update RHSof6.

      DO ih1 = 1, nps

!       find the cells in the current slab which are falling in a sphere with radius d_hdi

        ix = 1 + int ( ypn(1,ih1) / hx*real(ncdx) )
        iy = 1 + int ( (ypn(2,ih1) - real(indy)*hy) / hy*real(ncdy) )
        iz = 1 + int ( (ypn(3,ih1) - real(indz)*hz) / hz*real(ncdz) )
                                                                                                                                                                     
        if ( ix.gt.ncdx ) ix = ncdx
        if ( iy.gt.ncdy ) iy = ncdy
        if ( iz.gt.ncdz ) iz = ncdz
                                                                                                                                                                     
        ix_max = ix + ny_max
        ix_min = ix - ny_max
        iy_max = iy + ny_max
        iy_min = iy - ny_max
        iz_max = iz + nz_max
        iz_min = iz - nz_max

        if ( iy_max .gt. (ncdy+ny_max) ) iy_max = ncdy + ny_max
        if ( iy_min .lt. (1-ny_max) )    iy_min = 1 - ny_max

        if ( iz_max .gt. (ncdz+nz_max) ) iz_max = ncdz + nz_max
        if ( iz_min .lt. (1-nz_max) )    iz_min = 1 - nz_max

        do i = ix_min, ix_max
         do j = iy_min, iy_max
          do k = iz_min, iz_max
                                                                                                                                                                     
           ii = i
                                                                                                                                                                     
           sx = 0.0
           if ( i .gt. ncdx ) then
            ii = mod (i, ncdx)
            sx = 1.0
           elseif ( i .lt. 1 ) then
            ii = ncdx + i
            sx = -1.0
           endif
                                                                                                                                                                     
           ih2 = head (ii,j,k)
           IF ( ih2 .eq. 0 ) GOTO 2201

1199       CONTINUE

!          for all particles other than the target particle itself          
           IF ( ih2 .ne. ih1 ) THEN 
  
!           location of the second droplet
            xb_res(1) = ypn(1,ih2) + sx*hx
            xb_res(2) = ypn(2,ih2)
            xb_res(3) = ypn(3,ih2)

!           distance between the 2 droplets
            xd_res(:) = ypn(1:3,ih1) - xb_res(:)
            dnij  = sqrt ( xd_res(1)**2 + xd_res(2)**2 + xd_res(3)**2 )

            hdirad2 = HDI_trunc * ypn(7,ih2) 

            if ( dnij .lt. hdirad2 ) then

             n_inter = n_inter + 1

             fact1 = ypn(7,ih2) / dnij
             fact2 = fact1**3
             fact1 = 0.75 * fact1
 
             a1 = fact1 - 0.75 * fact2
             a2 = fact1 + 0.25 * fact2

             relv_res(:)    = ypres(1:3,ih2)
             relvb(:)       = ypres(1:3,ih2) + ypres(4:6,ih2)
             aloop          = a1 * DOT_PRODUCT ( relv_res,xd_res ) / dnij**2
             aloopb         = a1 * DOT_PRODUCT ( relvb,xd_res ) / dnij**2
             RHSof6(:,ih1)  = RHSof6(:,ih1) + aloop*xd_res(:) + a2*relv_res(:)
             RHSof6b(:,ih1) = RHSof6b(:,ih1) + aloopb*xd_res(:) + a2*relvb(:)

            endif

           ENDIF

           ih2 = list(ih2)
           IF ( ih2 .eq. 0 ) GOTO 2201

           GOTO 1199

2201       CONTINUE

          enddo 
         enddo
        enddo

      ENDDO

!-----normalized error, final stage.

      RHSminusLHS_c   = 0.0
      RHSb            = 0.0
      errtest         = 0.90

      DO ipart = 1, nps

       RHSminusLHS(:)  = RHSof6(:,ipart) - pertvel(:,ipart)
       RHSminusLHS_abs = RHSminusLHS(1)**2 + RHSminusLHS(2)**2 + RHSminusLHS(3)**2 
       RHSminusLHS_c   = RHSminusLHS_c + RHSminusLHS_abs

       RHSb            = RHSb + RHSof6b(1,ipart)**2 + RHSof6b(2,ipart)**2 + RHSof6b(3,ipart)**2 

!       if ( (id .eq. 10 ) .and. (ipart .lt. 50) ) write (*,*) 'iter, RHSminusLHS_abs = ', iter, RHSminusLHS_abs
!       if ( (id .eq. 10 ) .and. (ipart .gt. nps-50) ) write (*,*) 'iter, RHSminusLHS_abs = ', iter, RHSminusLHS_abs
!       if ( ( RHSminusLHS_abs .gt. errtest ) .and. (iter .gt. 2) .and. (istep .gt. 1) ) write (*,*) 'istep, iter, ipart, RHSminusLHS_abs, yp(1:3,ipart) = ', istep, iter, ipart, RHSminusLHS_abs, yp(1:3,ipart) 

      ENDDO

      CALL MPI_ALLREDUCE (RHSminusLHS_c,RHSminusLHS_sum,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE (RHSb,RHSb_sum,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE (n_inter,nsum_inter,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)

      error_norm  = sqrt ( RHSminusLHS_sum / RHSb_sum ) 

      if ( id .eq. 0 ) write (*,*) 'istep, iter, error_norm, n_inter/npart = ', istep, iter, error_norm, nsum_inter/npart

      DEALLOCATE ( RHSof6         )
      DEALLOCATE ( RHSof6b        )
      DEALLOCATE ( ypres          ) 
      DEALLOCATE ( temp_p, temp_b ) 

      RETURN
      END SUBROUTINE RESIDUAL 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE NEI_LIST_PRECOMP (id)
                                                                                                                                                         
!     This subroutine gathres the data of neighboring particles affecting a particle via HDI.
!     nei_num(ih1) is the number of neighboring particles affecting particle ih1.                 
!     nei_data(:,nei_num(ih1),ih1) are the hdi coefficients. See GMRES paper for details.

      USE params
      INCLUDE 'mpif.h'
                 
      REAL                     :: xb(3), xd(3), fact1, fact2, gamma, beta

      ALLOCATE ( nei_num(nps)       )
      ALLOCATE ( nei_data(7,100,nps) )

      stime    = mpi_wtime()
                                                                                                                                                         
      nei_num(1:nps)      = 0 
      nei_data(:,:,1:nps) = 0.0 

      DO ih1 = 1, nps
                                                                                                                                                         
!       find the cells in the current slab which are falling in a sphere with radius d_hdi
                                                                                                                                                         
        ix = 1 + int ( ypn(1,ih1) / hx*real(ncdx) )
        iy = 1 + int ( (ypn(2,ih1) - real(indy)*hy) / hy*real(ncdy) )
        iz = 1 + int ( (ypn(3,ih1) - real(indz)*hz) / hz*real(ncdz) )
                                                                                                                                                         
        if ( ix.gt.ncdx ) ix = ncdx
        if ( iy.gt.ncdy ) iy = ncdy
        if ( iz.gt.ncdz ) iz = ncdz
                                                                                                                                                         
        ix_max = ix + nx_max
        ix_min = ix - nx_max
        iy_max = iy + ny_max
        iy_min = iy - ny_max
        iz_max = iz + nz_max
        iz_min = iz - nz_max
                                                                                                                                                         
        if ( iy_max .gt. (ncdy+ny_max) ) iy_max = ncdy + ny_max
        if ( iy_min .lt. (1-ny_max) )    iy_min = 1 - ny_max
                                                                                                                                                         
        if ( iz_max .gt. (ncdz+nz_max) ) iz_max = ncdz + nz_max
        if ( iz_min .lt. (1-nz_max) )    iz_min = 1 - nz_max

        do i = ix_min, ix_max
         do j = iy_min, iy_max
          do k = iz_min, iz_max
                                                                                                                                                         
           ii = i
                                                                                                                                                         
           sx = 0.0
           if ( i .gt. ncdx ) then
            ii = mod (i, ncdx)
            sx = 1.0
           elseif ( i .lt. 1 ) then
            ii = ncdx + i
            sx = -1.0
           endif
                                                                                                                                                         
           ih2 = head (ii,j,k)
           IF ( ih2 .eq. 0 ) GOTO 2201
                                                                                                                                                         
1199       CONTINUE
                  
!          for all particles with greater indices
           IF ( ih2 .le. ih1 ) GOTO 2201
                  
!           location of the second droplet
            xb(1) = ypn(1,ih2) + sx*hx
            xb(2) = ypn(2,ih2)
            xb(3) = ypn(3,ih2)
                  
!           distance between the 2 droplets
            xd(:) = ypn(1:3,ih1) - xb(:)
            dnij  = sqrt ( xd(1)**2 + xd(2)**2 + xd(3)**2 )
                  
            if ( dnij .lt. (ypn(7,ih1)+ypn(7,ih2)) ) then
             write(*,*)'inside NEI_LIST_PRECOMP: overlapping is not allowed'
             write(*,*)'nps,npstest,ih1,ih2,dnij,ypn(1:8,ih1),ypn(1:8,ih2) =',nps,npstest,ih1,ih2,dnij,ypn(1:8,ih1),ypn(1:8,ih2)
             CALL mpi_finalize(ierr)
             STOP
            endif

!           for particle with index ih1
            hdirad = HDI_trunc * ypn(7,ih2)
!            hdirad = HDI_trunc * radmax 
              
            if ( dnij .lt. hdirad ) then
             nei_num(ih1) = nei_num(ih1) + 1 
             inum         = nei_num(ih1)
             if ( inum .gt. 100 ) write(*,*) 'error: size of nei_data is small' 

             nei_data(1,inum,ih1) = real(ih2) 

             fact1 = ypn(7,ih2) / dnij
             fact2 = fact1**3
             fact1 = 0.75 * fact1
                                                                                                                                                         
             gamma = ( fact1 - 0.75 * fact2 ) / dnij**2
             beta  = fact1 + 0.25 * fact2

             nei_data(2,inum,ih1) = gamma * xd(1)**2 + beta 
             nei_data(3,inum,ih1) = gamma * xd(1) * xd(2)
             nei_data(4,inum,ih1) = gamma * xd(1) * xd(3) 
             nei_data(5,inum,ih1) = gamma * xd(2)**2 + beta 
             nei_data(6,inum,ih1) = gamma * xd(2) * xd(3) 
             nei_data(7,inum,ih1) = gamma * xd(3)**2 + beta 
            endif

!           for particle with index ih2
            hdirad = HDI_trunc * ypn(7,ih1)
!            hdirad = HDI_trunc * radmax 
            xd     = -xd 

            if ( ( ih2 .le. nps) .and. ( dnij .lt. hdirad ) ) then
             nei_num(ih2) = nei_num(ih2) + 1 
             inum         = nei_num(ih2)
             if ( inum .gt. 100 ) write(*,*) 'error: size of nei_data is small' 

             nei_data(1,inum,ih2) = real(ih1) 

             fact1 = ypn(7,ih1) / dnij
             fact2 = fact1**3
             fact1 = 0.75 * fact1
                                                                                                                                                         
             gamma = ( fact1 - 0.75 * fact2 ) / dnij**2
             beta  = fact1 + 0.25 * fact2
                                                                                                                                                         
             nei_data(2,inum,ih2) = gamma * xd(1)**2 + beta
             nei_data(3,inum,ih2) = gamma * xd(1) * xd(2)
             nei_data(4,inum,ih2) = gamma * xd(1) * xd(3)
             nei_data(5,inum,ih2) = gamma * xd(2)**2 + beta
             nei_data(6,inum,ih2) = gamma * xd(2) * xd(3)
             nei_data(7,inum,ih2) = gamma * xd(3)**2 + beta
            endif

           ih2 = list(ih2)
           IF ( ih2 .le. ih1 ) GOTO 2201
                                                                                                                                                         
           GOTO 1199
                                                                                                                                                         
2201       CONTINUE

          enddo
         enddo
        enddo
 
      ENDDO

      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime

      RETURN
      END SUBROUTINE NEI_LIST_PRECOMP

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE AFUN_PRECOMP (x, id, icall)
                                                                                                                                                         
!     For any x, this subroutine computes Ax in which A is the matrix obtained by...
!     reforming linear HDI system into the form Ax = b. This is needed to build the Krylov subspace.
!     This subroutine is called from HDI_GMRES subroutine.
!     Written by Hossein Parishani, University of Delaware.
!     Questions, please contact hparish@udel.edu
                                                                                                                                                         
      USE params
      INCLUDE 'mpif.h'
                                                                                                                                                         
      REAL      :: xb(3), xd(3), relv(3) ,relvb(3)
      INTEGER   :: req(2), hdipart, error, ierr, icall
      INTEGER   :: status1(MPI_STATUS_SIZE), status2(MPI_STATUS_SIZE,2)
                                                                                                                                                         
      REAL,    DIMENSION (3*nps,1)               :: x
                                                                                                                                                         
      REAL,    ALLOCATABLE, DIMENSION (:,:)      :: Ax, xreshaped
      REAL,    ALLOCATABLE, DIMENSION (:,:)      :: RHSof6, bof6
      REAL,    ALLOCATABLE, DIMENSION (:,:,:)    :: temp_p, temp_b
                                                                                                                                                         
      stime   = mpi_wtime()
       
      hdipart = int ( 30 * nuniform )
       
      ALLOCATE ( xreshaped   (3,nps)   )
      ALLOCATE ( Ax          (3,nps)   )
      ALLOCATE ( RHSof6      (3,nps)   )
      ALLOCATE ( bof6        (3,nps)   )
      ALLOCATE ( temp_p  (3,hdipart,8) )
      ALLOCATE ( temp_b  (3,hdipart,8) )
       
      xreshaped(1, 1:nps) = x(1:nps,1)
      xreshaped(2, 1:nps) = x(nps+1:2*nps,1)
      xreshaped(3, 1:nps) = x(2*nps+1:3*nps,1)
       
      ypn (16:18, 1:nps) = xreshaped(1:3, 1:nps)

      DO jid = 1, 8
       DO jik = 1, ih_p(jid)
                                                                                                                                                         
        partind              = pindex (jik, jid)
        temp_p (1:3,jik,jid) = xreshaped(1:3, partind)
                  
       ENDDO
      ENDDO
            
      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime
      stime       = mpi_wtime()
                  
!-----Communicate the particles close to right & left subdomain boundaries
                  
      DO jid = 1, 8
                  
                  
        ids = nei_id ( jid )
        idr = nei_id ( mod (jid+3, 8) + 1 )
                  
                  
        if ( ih_p(jid) .gt. 0 ) then
         CALL MPI_SEND (temp_p(:,1:ih_p(jid),jid),3*ih_p(jid),MPI_REAL8,ids,213,MPI_COMM_WORLD,ierr)
        endif
                  
                  
        if ( ih_pt(jid) .gt. 0 ) then
         CALL MPI_RECV (temp_b(:,1:ih_pt(jid),jid),3*ih_pt(jid),MPI_REAL8,idr,213,MPI_COMM_WORLD,status1,ierr)
        endif
                  
                  
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAITALL (2,req,status2,ierr)
!        if ( ih_p(jid).gt.0 .and. ih_pt(jid).eq.0 ) CALL MPI_WAIT (req(1),status1,ierr)
!        if ( ih_p(jid).eq.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAIT (req(2),status1,ierr)

        if ( ih_pt(jid) .gt. 0 ) then
                                                                                                                                                         
          nps1 = nps + SUM ( ih_pt ( 1 : (jid-1) ) )
          nps2 = nps + SUM ( ih_pt ( 1 : jid ) )
                                                                                                                                                         
          ypn (16:18, (nps1+1):nps2) = temp_b(1:3,1:ih_pt(jid),jid)
       
        endif
       
      ENDDO
       
      etime       = mpi_wtime()
      time_hdi_cm = time_hdi_cm + etime - stime
      stime       = mpi_wtime()
       
!-----Compute Ax and b with precomputed coefficients. See GMRES paper for details.
       
      RHSof6 = 0.0
      bof6   = 0.0
       
      DO ih1 = 1, nps
       DO ijk = 1, nei_num(ih1)

         ind = int ( nei_data(1,ijk,ih1) )

         relv(:)       = ypn(16:18,ind)
         RHSof6(1,ih1) = RHSof6(1,ih1) + nei_data(2,ijk,ih1)*relv(1) + nei_data(3,ijk,ih1)*relv(2) + nei_data(4,ijk,ih1)*relv(3) 
         RHSof6(2,ih1) = RHSof6(2,ih1) + nei_data(3,ijk,ih1)*relv(1) + nei_data(5,ijk,ih1)*relv(2) + nei_data(6,ijk,ih1)*relv(3) 
         RHSof6(3,ih1) = RHSof6(3,ih1) + nei_data(4,ijk,ih1)*relv(1) + nei_data(6,ijk,ih1)*relv(2) + nei_data(7,ijk,ih1)*relv(3) 

         if ( icall .eq. 1 ) then
          relvb(:)    = ypn(4:6,ind) - vpn(1:3,ind)
          bof6(1,ih1) = bof6(1,ih1) + nei_data(2,ijk,ih1)*relvb(1) + nei_data(3,ijk,ih1)*relvb(2) + nei_data(4,ijk,ih1)*relvb(3) 
          bof6(2,ih1) = bof6(2,ih1) + nei_data(3,ijk,ih1)*relvb(1) + nei_data(5,ijk,ih1)*relvb(2) + nei_data(6,ijk,ih1)*relvb(3) 
          bof6(3,ih1) = bof6(3,ih1) + nei_data(4,ijk,ih1)*relvb(1) + nei_data(6,ijk,ih1)*relvb(2) + nei_data(7,ijk,ih1)*relvb(3) 
         endif

       ENDDO
      ENDDO
                                                                                                                                                         
      Ax(:,1:nps) = RHSof6(:,1:nps) + xreshaped(:,1:nps)
                                                                                                                                                         
      Axreshaped(1:nps,1)         = Ax (1, 1:nps)
      Axreshaped(nps+1:2*nps,1)   = Ax (2, 1:nps)
      Axreshaped(2*nps+1:3*nps,1) = Ax (3, 1:nps)
                                                                                                                                                         
      breshaped(1:nps,1)          = bof6 (1, 1:nps)
      breshaped(nps+1:2*nps,1)    = bof6 (2, 1:nps)
      breshaped(2*nps+1:3*nps,1)  = bof6 (3, 1:nps)
                                                                                                                                                         
      DEALLOCATE ( RHSof6     )
      DEALLOCATE ( temp_p     )
      DEALLOCATE ( temp_b     )
      DEALLOCATE ( bof6       )
      DEALLOCATE ( Ax         )
      DEALLOCATE ( xreshaped  )
                                                                                                                                                         
      etime       = mpi_wtime()
      time_hdi_cp = time_hdi_cp + etime - stime
                                                                                                                                                         
      RETURN
      END SUBROUTINE AFUN_PRECOMP

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE PART_ACCELERATION (id, idump, destination)
                                                                                                                                                                
!     This subroutine computes and dumps particle acceleration.
                                                                                                                                                                
      USE params
      INCLUDE 'mpif.h'
                                                                                                                                                                
      INTEGER   :: req(2), nsize, error, ierr, partind
      INTEGER   :: status1(MPI_STATUS_SIZE), status2(MPI_STATUS_SIZE,2)
      INTEGER   :: npsgather(0:nproc-1)
      REAL      :: taup
                                                                                                                                                                
      REAL, ALLOCATABLE, DIMENSION    (:,:)     :: acc
      REAL, ALLOCATABLE, DIMENSION    (:,:)     :: acc_data, wxyzrt
                                                                                                                                                                
      CHARACTER*50 destination
      CHARACTER*70 fnm1, fnm2
                                                                                                                                                                
      ALLOCATE ( acc  (4,nps) )
      ALLOCATE ( acc_data (1:3, npart) )
                                                                                                                                                                
      acc_data (:,:) = 0.
                                                                                                                                                                
      DO ih = 1, nps
                                                                                                                                                                
       taup = tau_CLOUD_const2 * yp(7,ih)**2
!      if (id.eq.0) write(*,*) 'taup = ', taup
                                                                                                                                                                
       IF ( HDI_included ) THEN
        acc (1, ih) = - ( yp(4, ih) - vp(1, ih) - pertvel(1, ih) ) / taup + gravity_dns
        acc (2, ih) = - ( yp(5, ih) - vp(2, ih) - pertvel(2, ih) ) / taup
        acc (3, ih) = - ( yp(6, ih) - vp(3, ih) - pertvel(3, ih) ) / taup
        acc (4, ih) = yp(7, ih)
       ELSE
        acc (1, ih) = - ( yp(4, ih) - vp(1, ih) ) / taup + gravity_dns
        acc (2, ih) = - ( yp(5, ih) - vp(2, ih) ) / taup
        acc (3, ih) = - ( yp(6, ih) - vp(3, ih) ) / taup
        acc (4, ih) = yp(7, ih)
       ENDIF

      ENDDO
                                                                                                                                                                
!     if (id.eq.0) write(*,*) 'acc (1,1:5) = ', acc (1,1:5)
                                                                                                                                                                
      CALL MPI_GATHER (nps,1,MPI_INTEGER,npsgather,1,MPI_INTEGER,0,mpi_comm_world,ierr)

      if (id.eq.0) write(*,*) '==== gathering and dumping particle accelaration data ===='
!     if (id.eq.0) write(*,*) 'npart, SUM (npsgather), gravity_dns = ', npart, SUM (npsgather), gravity_dns                                                                                                                                                                
      ipart1 = 0
      ipart2 = 0
                                                                                                                                                                
      IF (id.eq.0) THEN
                                                                                                                                                                
       do i = 1, nps
                                                                                                                                                                
        if ( yp(7,i) .eq. radius(1) ) then
         ipart1 = ipart1 + 1
         acc_data (1:3,ipart1) = acc (1:3,i)
        elseif ( yp(7,i) .eq. radius(2) ) then
         ipart2 = ipart2 + 1
         acc_data (1:3,npart/2+ipart2) = acc (1:3,i)
        endif
                                                                                                                                                                
       enddo
                                                                                                                                                                
       do n=1, nproc-1
                                                                                                                                                                
        ALLOCATE ( wxyzrt(4,npsgather(n)) )
                                                                                                                                                                
        CALL MPI_RECV(wxyzrt,4*npsgather(n),MPI_REAL8,n,111,MPI_COMM_WORLD,status1,ierr)
                                                                                                                                                                
        do i = 1, npsgather(n)
                                                                                                                                                                
         if ( wxyzrt(4,i) .eq. radius(1) ) then
          ipart1 = ipart1 + 1
          acc_data (1:3,ipart1) = wxyzrt(1:3,i)
         elseif ( wxyzrt(4,i) .eq. radius(2) ) then
          ipart2 = ipart2 + 1
          acc_data (1:3,npart/2+ipart2) = wxyzrt(1:3,i)
         endif
                                                                                                                                                                
        enddo
                                                                                                                                                                
         DEALLOCATE ( wxyzrt )
                                                                                                                                                                
       enddo

      ELSE
       CALL MPI_SEND (acc(1:4,1:nps),4*nps,MPI_REAL8,0,111,MPI_COMM_WORLD,status1,ierr)
      ENDIF
                                                                                                                                                                
      CALL MPI_BARRIER ( MPI_COMM_WORLD, error )
                                                                                                                                                                
      IF (id.eq.0) THEN
                                                                                                                                                                
       i1step=int(istep/100000)
       i2step=int((istep-i1step*100000)/10000)
       i3step=int((istep-i1step*100000-i2step*10000)/1000)
       i4step=int((istep-i1step*100000-i2step*10000-i3step*1000)/100)
       i5step=int((istep-i1step*100000-i2step*10000-i3step*1000-i4step*100)/10)
       i6step=istep-i1step*100000-i2step*10000-i3step*1000-i4step*100-i5step*10
                                                                                                                                                                
       if(idump.lt.10) then
        fnm1=TRIM(destination)//'part_acc'//char(idump+48)//'.'//char(i1step+48) &
             //char(i2step+48)//char(i3step+48)//char(i4step+48)//char(i5step+48)//char(i6step+48)
       else
        i4dump=int(idump/10)
        i5dump=mod(idump,10)
        fnm1=TRIM(destination)//'part_acc'//char(i4dump+48)//char(i5dump+48)//'.'//char(i1step+48) &
             //char(i2step+48)//char(i3step+48)//char(i4step+48)//char(i5step+48)//char(i6step+48)
       endif
                                                                                                                                                                
       open(17,file=fnm1,status='unknown')
       write(17,6158) acc_data (1:3,1:npart)
       close(17)
                                                                                                                                                                
6158   format ( 3(1pe14.6) )
                                                                                                                                                                
      ENDIF
                                                                                                                                                                
      DEALLOCATE ( acc_data )
      DEALLOCATE ( acc )
                                                                                                                                                                
      RETURN
      END SUBROUTINE PART_ACCELERATION

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE COALESCE (id)

!     This subroutine substitutes the colliding droplets with a new droplet...
!     ... and computes the new droplet rasius, velocity and location.

      USE params
      INCLUDE 'mpif.h'

      INTEGER           :: sindex, eindex, kktot, sremove, p1ind, p2ind
      INTEGER           :: p1id, p2id, pindex1, pindex2, part1id, part2id 
      INTEGER           :: addindex, collmax, buffsize, kk1, error, ierr

!      REAL, ALLOCATABLE, DIMENSION (:,:)      :: clist_t
!      REAL, ALLOCATABLE, DIMENSION (:)        :: collist_t

      if ( id .eq. 0 ) write (*,*) 'id=0, coalescence subroutine flag 0'

!      collmax  = ncollmax * nproc 
!      buffsize = 20 * ncollmax * nproc 

!      ALLOCATE ( clist_t  (20, collmax) )
!      ALLOCATE ( collist_t  (buffsize) )

!      CALL MPI_ALLGATHER (ncollpair,1,MPI_INTEGER,ncollpair_t,1,MPI_INTEGER,mpi_comm_world)
!      if ( id .eq. 0 ) write (*,*) 'collist = ', collist
!      CALL MPI_ALLGATHER (collist,20*ncollmax,MPI_REAL8,collist_t,20*ncollmax,MPI_REAL8,mpi_comm_world)

!      CALL MPI_GATHER (collist,20*ncollmax,MPI_REAL8,collist_t,20*ncollmax,MPI_REAL8,0,mpi_comm_world,err)
!      CALL MPI_BCAST (collist_t,buffsize,MPI_REAL8,0,mpi_comm_world,ierr)

!      -----------------------------------------------------------

      jj = 0

      DO i = 1, ntotpair 

         p1ind = int ( collist (9, i) )
         p2ind = int ( collist (19, i) )

!        replace the 1st colliding particle with the daughter particle
         if ( p1ind .le. nps ) then

           yp (7, p1ind)   = ( collist (7, i)**3.0 + collist (17, i)**3.0 ) ** (1.0/3.0)

!          to conserve total momentum 
!          yp (4:6, p1ind) = ( collist(4:6,i)*collist(7,i)**3.0 + collist(14:16,i)*collist(17,i)**3.0 ) / yp(7,p1ind)**3  

!          falling with terminal velocity. no momentum conservation
           yp (4, p1ind) = Sv_const * yp(7,p1ind)**2
           yp (5, p1ind) = 0.0
           yp (6, p1ind) = 0.0

         endif

!        remove the 2nd particle in the host domain, then replace it
!        with the last particle in the list
         if ( p2ind .le. nps ) then

           jj = jj + 1

           yp (:, p2ind)    = yp (:, nps-jj+1)
           vp(:, p2ind)     = vp(:, nps-jj+1)
           dmove(:, p2ind)  = dmove(:, nps-jj+1) 
           fvp(:, 0, p2ind) = fvp(:, 0, nps-jj+1)
           fvp(:, 1, p2ind) = fvp(:, 1, nps-jj+1)
           fvp(:, 2, p2ind) = fvp(:, 2, nps-jj+1)
           fdp(:, 0, p2ind) = fdp(:, 0, nps-jj+1)
           fdp(:, 1, p2ind) = fdp(:, 1, nps-jj+1)
           fdp(:, 2, p2ind) = fdp(:, 2, nps-jj+1)
           fdp(:, 3, p2ind) = fdp(:, 3, nps-jj+1)

           yp (:, nps-jj+1)    = 0.0
           vp(:, nps-jj+1)     = 0.0
           dmove(:, nps-jj+1)  = 0.0
           fvp(:, 0, nps-jj+1) = 0.0
           fvp(:, 1, nps-jj+1) = 0.0
           fvp(:, 2, nps-jj+1) = 0.0
           fdp(:, 0, nps-jj+1) = 0.0
           fdp(:, 1, nps-jj+1) = 0.0
           fdp(:, 2, nps-jj+1) = 0.0
           fdp(:, 3, nps-jj+1) = 0.0

         endif

      ENDDO

      nps = nps - jj

      if (id .eq. 0) write (*,*)'id=0, Leaving coalescence subroutine flag 1'

!      -----------------------------------------------------------

!      DO i = 1, nproc 

!       sindex = SUM ( npair_t (1:i-1) )
!       eindex = SUM ( npair_t (1:i) )
!       kk1 = (i-1) * ncollmax * 20

!       DO j = sindex+1, eindex 

!        clist_t (1:20, j) = collist_t ( kk1+(j-1)*20+1 : kk1+j*20)

!       ENDDO
!      ENDDO

!      npair_total = SUM ( npair_t (1:nproc) )
!      kktot = 0

!      DO i = 1, npair_total-1
!        DO j = i+1, npair_total
                                                                                                                                                
!         gind1 = int ( clist_t (8, i) ) 
!         gind2 = int ( clist_t (18, i) )
!         gind3 = int ( clist_t (8, j) )
!         gind4 = int ( clist_t (18, j) )

!         sremove = 0
!         if ( ( gind1 .eq. gind3) .OR. ( gind1 .eq. gind4) ) sremove = 1 
!         if ( ( gind2 .eq. gind3) .OR. ( gind2 .eq. gind4) ) sremove = 1 

!         if ( sremove .eq. 1) then
!          clist_t (1:20, j) = clist_t (1:20, npair_total-kktot)
!          kktot = kktot + 1
!         endif

!        ENDDO
!      ENDDO

!      npair_new = npair_total - kktot

!      jj = 0

!      DO  i = 1, npair_new 

!          p1id = int ( clist_t (10, i) )
!          p2id = int ( clist_t (20, i) )

!          IF ( p1id .eq. p2id ) then
           
!            if ( p1id .eq. id ) then

!             jj = jj + 1
!             p1ind = int ( clist_t (9, i) )
!             p2ind = int ( clist_t (19, i) )

!             pindex1 = min ( p1ind, p2ind )
!             pindex2 = max ( p1ind, p2ind )

!             yp (7, pindex1) = ( clist_t (7, i)**3.0 + clist_t (17, i)**3.0 ) ** (1.0/3.0) 
!             yp (:, pindex2) = yp (:, nps-jj+1) 
!             yp (:, nps-jj+1) = 0.0 

!            endif

!          ELSE
        
!             p1ind = int ( clist_t (9, i) )
!             p2ind = int ( clist_t (19, i) )
                                                                                                                                                
!             pindex1 = min ( p1ind, p2ind )
!             pindex2 = max ( p1ind, p2ind )
                                                                                                                                                
!             part1id = int ( yp (10, pindex1) )
!             part2id = int ( yp (10, pindex2) )

!             if ( part1id .eq. id ) yp (7, pindex1) = ( clist_t (7, i)**3.0 + clist_t (17, i)**3.0 ) ** (1.0/3.0) 

!             if ( part2id .eq. id ) then
!              jj = jj + 1
!              yp (:, pindex2)  = yp (:, nps-jj+1) 
!              yp (:, nps-jj+1) = 0.0 
!             endif

!          ENDIF
         
!      ENDDO

!      nps = nps - jj

!      DEALLOCATE ( clist_t   )
!      DEALLOCATE ( collist_t )

      RETURN
      END SUBROUTINE COALESCE 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE SIZE_DIST_LINEAR (id)

!     This subroutine computes the size distribution of the droplets for growing case.
!     ... the size range is divided into bins. then the distribution is gathered, summed up and printed.

      USE params
      INCLUDE 'mpif.h'

      INTEGER                                 ::  nnbin, ierr, ir 
      INTEGER, ALLOCATABLE, DIMENSION  (:)    ::  sizebin_0, sizebin_t

      nnbin = 120

      ALLOCATE ( sizebin_0 (1:nnbin) )
      ALLOCATE ( sizebin_t (1:nnbin) )

      sizebin_0 = 0
      sizebin_t = 0

      DO  i = 1, nps

       ir = int ( ( yp (7,i) / radius(1) )**3.0 + 0.50 )

       if ( ir .lt. 1 )     ir = 1 
       if ( ir .gt. nnbin ) ir = nnbin

       sizebin_0 (ir) = sizebin_0 (ir) + 1

      ENDDO

      CALL MPI_ALLREDUCE (sizebin_0, sizebin_t, nnbin, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)

      if ( id.eq.0 ) write (1001, 515) istep, npart, npart_t, sizebin_t(1:nnbin)

515   format ( 2x, I10, I10, I10, 120(I10) )

      DEALLOCATE ( sizebin_0 )
      DEALLOCATE ( sizebin_t )

      RETURN
      END SUBROUTINE SIZE_DIST_LINEAR 

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE SIZE_DIST_LOGARITHMIC (id)

!     This subroutine computes the size distribution of the droplets for
!     growing case.
!     ... the size range is divided into bins. then the distribution is
!     gathered, summed up and printed.

      USE params
      INCLUDE 'mpif.h'

      REAL                                    ::  m_ratio, scal, rho_bin
      INTEGER                                 ::  nnbin, ierr, ir, error
      INTEGER, ALLOCATABLE, DIMENSION  (:)    ::  sizebin_0, sizebin_t
      REAL, ALLOCATABLE, DIMENSION  (:)       ::  mdensbin_0, mdensbin_t

      nnbin   = 360
      scal    = 32.0
      rho_bin = 2.0 ** ( 1.0 / scal ) 

      ALLOCATE ( sizebin_0 (1:nnbin) )
      ALLOCATE ( sizebin_t (1:nnbin) )
      ALLOCATE ( mdensbin_0 (1:nnbin) )
      ALLOCATE ( mdensbin_t (1:nnbin) )

      sizebin_0 = 0
      sizebin_t = 0

      mdensbin_0 = 0.0
      mdensbin_t = 0.0

      DO  i = 1, nps

       m_ratio = ( yp (7,i) / radius(1) )**3.0
       ir = int ( ( log (m_ratio) / log (rho_bin) ) + 1.5 )

       if ( ir .lt. 1 )     ir = 1
       if ( ir .gt. nnbin ) ir = nnbin

       sizebin_0 (ir) = sizebin_0 (ir) + 1
       mdensbin_0 (ir) = mdensbin_0 (ir) + (4.0/3.0) * rho_water * pi * yp (7,i)**3 

      ENDDO

      CALL MPI_ALLREDUCE (sizebin_0, sizebin_t, nnbin, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
      CALL MPI_ALLREDUCE (mdensbin_0, mdensbin_t, nnbin, MPI_REAL8, MPI_SUM, mpi_comm_world, ierr)

      mdensbin_t = 3.0 * scal * mdensbin_t / log (2.0)  

      if ( id.eq.0 ) write (1001, 515) istep, npart, npart_t, sizebin_t(1:nnbin)
      if ( id.eq.0 ) write (1002, 516) istep, npart, npart_t, mdensbin_t(1:nnbin)

515   format ( 2x, I10, I10, I10, 360(I10) )
516   format ( 2x, I10, I10, I10, 360(1pe14.6) )

      DEALLOCATE ( sizebin_0 )
      DEALLOCATE ( sizebin_t )
      DEALLOCATE ( mdensbin_0 )
      DEALLOCATE ( mdensbin_t )

      RETURN
      END SUBROUTINE SIZE_DIST_LOGARITHMIC

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE SIZE_AVG (id)

!     This subroutine computes the average radius of the droplets.

      USE params
      INCLUDE 'mpif.h'

      REAL                                    ::  r_sum_0, r_max_0
      REAL                                    ::  r_sum_t, r_max_t
      REAL                                    ::  r_avg, r_max_avg
      INTEGER                                 ::  ierr, error

      r_sum_0 = 0.0 
      r_max_0 = 0.0 

      DO  i = 1, nps
       r_sum_0 = r_sum_0 + yp (7,i)**3.0 
      ENDDO

      r_max_0 = ( MAXVAL ( yp (7,1:nps) ) )**3.0

      CALL MPI_ALLREDUCE (r_sum_0, r_sum_t, 1, MPI_REAL8, MPI_SUM, mpi_comm_world, ierr)
      CALL MPI_ALLREDUCE (r_max_0, r_max_t, 1, MPI_REAL8, MPI_SUM, mpi_comm_world, ierr)

      r_avg     = ( r_sum_t / real(npart_t) ) ** (1.0/3.0)
      r_avg     = r_avg / radius(1)

      r_max_avg = ( r_max_t / real(nproc) ) ** (1.0/3.0)
      r_max_avg = r_max_avg / radius(1)

      if ( id.eq.0 ) write (1003, 517) istep, npart, npart_t, r_avg, r_max_avg 

517   format ( 2x, I10, I10, I10, 2(1pe14.6) )

      RETURN
      END SUBROUTINE SIZE_AVG

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE DISSIP(ux,uy,uz,kkx,kky,kkz,rnul,edissip,id)

!     This subroutine computes the local dissipation rate when called.

      USE params
      INCLUDE 'mpif.h'

      REAL,DIMENSION      (llx,lly,lz)  ::  ux,uy,uz
      REAL,DIMENSION      (llx,lly,lz)  ::  wx1,wy1
      REAL,DIMENSION      (llx,lly,lz)  ::  kkx,kky,kkz
      REAL,DIMENSION      (llx,lly,lz)  ::  tmph
      REAL,DIMENSION      (llx,lly,lz)  ::  edissip

! e12--------------------------------------------------------------

      wx1 = kky * ux   !du/dy
      wy1 = kkx * uy   !dv/dx
      tmph = wx1
      wx1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wy1
      wy1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wy1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
       
      tmph=0.5*(wx1+wy1)

!     Parallel FFT to real space
      CALL mpifft3DCR (tmph,id)
 
      edissip=2.0*tmph*tmph

! e23--------------------------------------------------------------
      wx1 = kkz * uy   !dv/dz
      wy1 = kky * uz   !dw/dy
      tmph = wx1
      wx1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wy1
      wy1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wy1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)

      tmph=0.5*(wx1+wy1)

!     Parallel FFT
      CALL mpifft3DCR (tmph,id)

      edissip=edissip+2.0*tmph*tmph

! e31--------------------------------------------------------------
      wx1 = kkx * uz   !dw/dx
      wy1 = kkz * ux   !du/dz
      tmph = wx1
      wx1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wy1
      wy1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wy1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)

      tmph=0.5*(wx1+wy1)

!     Parallel FFT
      CALL mpifft3DCR (tmph,id)

      edissip=edissip+2.0*tmph*tmph

! e11--------------------------------------------------------------
      wx1 = kkx * ux   !du/dx
      tmph = wx1
      wx1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)

      tmph=wx1

!     Parallel FFT to real space
      CALL mpifft3DCR (tmph,id)

      edissip=edissip+tmph*tmph

! e22--------------------------------------------------------------
      wx1 = kky * uy   !dv/dy
      tmph = wx1
      wx1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)

      tmph=wx1

!     Parallel FFT to real space
      CALL mpifft3DCR (tmph,id)

      edissip=edissip+tmph*tmph

! e33--------------------------------------------------------------
      wx1 = kkz * uz   !dw/dz
      tmph = wx1
      wx1(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx1(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)

      tmph=wx1

!     Parallel FFT to real space
      CALL mpifft3DCR (tmph,id)

      edissip=edissip+tmph*tmph

      edissip=2.0*rnul*edissip

      END SUBROUTINE DISSIP

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE COMP_ENSTROPHY(wx,wy,wz,enstrophy)

      USE params
      INCLUDE 'mpif.h'

      REAL,DIMENSION      (llx,lly,lz)  ::  wx,wy,wz
      REAL,DIMENSION      (llx,lly,lz)  ::  enstrophy

      enstrophy=SQRT(wx**2+wy**2+wz**2)


      END SUBROUTINE COMP_ENSTROPHY

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE INTERP_SCALAR (v, jjj, id)

!     This subroutine does the 6-pt Lagrange interpolation of...
!     the data field u using interpolation factors bg.


      USE params
      INCLUDE 'mpif.h'


      INTEGER                              :: id,jjj
      INTEGER, DIMENSION (6)               :: ix, iy, iz
      INTEGER                              :: status1(MPI_STATUS_SIZE), error, ierr

      REAL, DIMENSION  (llx,lly,lz)        :: v
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: vext
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: v1s, v1r, v2s, v2r, v3s, v3r, v4s, v4r
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: v5s, v5r, v6s, v6r, v7s, v7r, v8s, v8r

      stime = mpi_wtime()

      ALLOCATE ( v1s  (lx,3,lz)  )
      ALLOCATE ( v1r  (lx,3,lz)  )
      ALLOCATE ( v2s  (lx,3,3)   )
      ALLOCATE ( v2r  (lx,3,3)   )
      ALLOCATE ( v3s  (lx,ly,3)  )
      ALLOCATE ( v3r  (lx,ly,3)  )
      ALLOCATE ( v4s  (lx,2,3)   )
      ALLOCATE ( v4r  (lx,2,3)   )
      ALLOCATE ( v5s  (lx,2,lz)  )
      ALLOCATE ( v5r  (lx,2,lz)  )
      ALLOCATE ( v6s  (lx,2,2)   )
      ALLOCATE ( v6r  (lx,2,2)   )
      ALLOCATE ( v7s  (lx,ly,2)  )
      ALLOCATE ( v7r  (lx,ly,2)  )
      ALLOCATE ( v8s  (lx,3,2)   )
      ALLOCATE ( v8r  (lx,3,2)   )
      ALLOCATE ( vext (lx,ly+5,lz+5) )

!-------------------------------------------------------------------

      etime          = mpi_wtime()
      time_interp_cp = time_interp_cp + etime - stime
      stime          = mpi_wtime()

      ls1 = 3 * lz * lx
      ls2 = 3 * 3  * lx
      ls3 = 3 * ly * lx
      ls4 = 2 * 3  * lx
      ls5 = 2 * lz * lx
      ls6 = 2 * 2  * lx
      ls7 = 2 * ly * lx
      ls8 = 2 * 3  * lx

      v1s ( 1:lx, 1:3, 1:lz ) = v (1:lx, 1:3, 1:lz)
      v2s ( 1:lx, 1:3, 1:3 )  = v (1:lx, 1:3, 1:3)
      v3s ( 1:lx, 1:ly, 1:3 ) = v (1:lx, 1:ly, 1:3)
      v4s ( 1:lx, 1:2, 1:3 )  = v (1:lx, (ly-1):ly, 1:3)
      v5s ( 1:lx, 1:2, 1:lz ) = v (1:lx, (ly-1):ly, 1:lz)
      v6s ( 1:lx, 1:2, 1:2 )  = v (1:lx, (ly-1):ly, (lz-1):lz)
      v7s ( 1:lx, 1:ly, 1:2 ) = v (1:lx, 1:ly, (lz-1):lz)
      v8s ( 1:lx, 1:3, 1:2 )  = v (1:lx, 1:3, (lz-1):lz)

! right
      CALL MPI_SENDRECV (v1s,ls1,MPI_REAL8,id_w,311,v1r,ls1,&
             MPI_REAL8,id_e,311,MPI_COMM_WORLD,status1,ierr1)

! left
      CALL MPI_SENDRECV (v5s,ls5,MPI_REAL8,id_e,315,v5r,ls5,&
             MPI_REAL8,id_w,315,MPI_COMM_WORLD,status1,ierr5)

! top
      CALL MPI_SENDRECV (v3s,ls3,MPI_REAL8,id_s,313,v3r,ls3,&
             MPI_REAL8,id_n,313,MPI_COMM_WORLD,status1,ierr3)

! bottom
      CALL MPI_SENDRECV (v7s,ls7,MPI_REAL8,id_n,317,v7r,ls7,&
             MPI_REAL8,id_s,317,MPI_COMM_WORLD,status1,ierr7)

! se corner
      CALL MPI_SENDRECV (v8s,ls8,MPI_REAL8,id_nw,318,v8r,ls8,&
             MPI_REAL8,id_se,318,MPI_COMM_WORLD,status1,ierr8)

! ne corner
      CALL MPI_SENDRECV (v2s,ls2,MPI_REAL8,id_sw,312,v2r,ls2,&
             MPI_REAL8,id_ne,312,MPI_COMM_WORLD,status1,ierr2)

! nw corner
      CALL MPI_SENDRECV (v4s,ls4,MPI_REAL8,id_se,314,v4r,ls4,&
             MPI_REAL8,id_nw,314,MPI_COMM_WORLD,status1,ierr4)

! sw corner
      CALL MPI_SENDRECV (v6s,ls6,MPI_REAL8,id_ne,316,v6r,ls6,&
             MPI_REAL8,id_sw,316,MPI_COMM_WORLD,status1,ierr6)

      vext ( 1:lx, 3:(ly+2), 3:(lz+2) )           = v   (1:lx, 1:ly, 1:lz)
      vext ( 1:lx, (ly+3):(ly+5), 3:(lz+2) )      = v1r (1:lx, 1:3, 1:lz)
      vext ( 1:lx, (ly+3):(ly+5), (lz+3):(lz+5) ) = v2r (1:lx, 1:3, 1:3)
      vext ( 1:lx, 3:(ly+2), (lz+3):(lz+5) )      = v3r (1:lx, 1:ly, 1:3)
      vext ( 1:lx, 1:2, (lz+3):(lz+5) )           = v4r (1:lx, 1:2, 1:3)
      vext ( 1:lx, 1:2, 3:(lz+2) )                = v5r (1:lx, 1:2, 1:lz)
      vext ( 1:lx, 1:2, 1:2 )                     = v6r (1:lx, 1:2, 1:2)
      vext ( 1:lx, 3:(ly+2), 1:2 )                = v7r (1:lx, 1:ly, 1:2)
      vext ( 1:lx, (ly+3):(ly+5), 1:2 )           = v8r (1:lx, 1:3, 1:2)

      etime          = mpi_wtime()
      time_interp_cm = time_interp_cm + etime - stime
      stime          = mpi_wtime()

!-------------------------------------------------------------------

      mn1 = lx

      if (jjj.eq.1) enst_p(1:nps) = 0.0
      if (jjj.eq.2) disp_p(1:nps) = 0.0

      DO ih0 = 1, nps

        jx    = lhnode(ih0,1) - 2
        ix(1) = ( mod (jx, mn1) + mn1*(1-isign(1,jx))/2 ) + 1

        if ( ix(1) .le. (mn1-5) ) then
          ix(2) = ix(1) + 1
          ix(3) = ix(2) + 1
          ix(4) = ix(3) + 1
          ix(5) = ix(4) + 1
          ix(6) = ix(5) + 1
        else
          ix(2) = mod (ix(1), mn1) + 1
          ix(3) = mod (ix(2), mn1) + 1
          ix(4) = mod (ix(3), mn1) + 1
          ix(5) = mod (ix(4), mn1) + 1
          ix(6) = mod (ix(5), mn1) + 1
        endif

        iy(1) = lhnode(ih0,2) - indy * ly + 1
        iy(2) = iy(1) + 1
        iy(3) = iy(2) + 1
        iy(4) = iy(3) + 1
        iy(5) = iy(4) + 1
        iy(6) = iy(5) + 1

        iz(1) = lhnode(ih0,3) - indz * lz + 1
        iz(2) = iz(1) + 1
        iz(3) = iz(2) + 1
        iz(4) = iz(3) + 1
        iz(5) = iz(4) + 1
        iz(6) = iz(5) + 1

        do  ndx = 1, 6
         ix1 = ix(ndx)
         do  ndy = 1, 6
          iy1 = iy(ndy)
          do  ndz = 1, 6
           iz1 = iz(ndz)

            if (jjj.eq.1) then 
               enst_p(ih0) = enst_p(ih0) + vext(ix1,iy1,iz1) * bg(ih0,ndx,1) * bg(ih0,ndy,2) * bg(ih0,ndz,3)
            endif

            if (jjj.eq.2) then
               disp_p(ih0) = disp_p(ih0) + vext(ix1,iy1,iz1) * bg(ih0,ndx,1) * bg(ih0,ndy,2) * bg(ih0,ndz,3)
            endif

          enddo
         enddo
        enddo


        if ((jjj.eq.1).AND.(enst_p(ih0).lt.0.)) enst_p(ih0) = 0.00001
!      if ((jjj.eq.1).AND.(enst_p(ih0).lt.0.))then
!        write(*,*)'id,enst_p lt than 0',id,enst_p(ih0),ih0
!       endif


        if ((jjj.eq.2).AND.(disp_p(ih0).lt.0.)) disp_p(ih0) = 0.00001
!       if ((jjj.eq.2).AND.(disp_p(ih0).lt.0.))then
!        write(*,*)'id,disp_p lt than 0',id,disp_p(ih0),ih0
!       endif

      END DO

      DEALLOCATE ( v1s, v1r, v2s, v2r, v3s, v3r, v4s, v4r, vext )
      DEALLOCATE ( v5s, v5r, v6s, v6r, v7s, v7r, v8s, v8r )

      etime          = mpi_wtime()
      time_interp_cp = time_inter_cp + etime - stime

      RETURN
      END SUBROUTINE INTERP_SCALAR

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      SUBROUTINE LOCAL_St_AVG (id)

!     This subroutine computes the average local St of the droplets.

      USE params
      INCLUDE 'mpif.h'

      REAL                      ::  st_sum_0, st_max_0, st_1
      REAL                      ::  st_sum_t, st_max_t
      REAL                      ::  st_avg, st_max_avg
      INTEGER                   ::  ierr, error

      st_sum_0 = 0.0
      st_max_0 = 0.0

      DO  i = 1, nps
       taup = tau_CLOUD_const2 * yp(7,i)**2
       st_sum_0 = st_sum_0 + taup * (disp_p(i)/rnu)**0.50
      ENDDO

      ind_max_r = MAXLOC ( yp (7,1:nps), DIM=1 )
      st_max_0 = tau_CLOUD_const2 * yp (7,ind_max_r)**2.0 * (disp_p(ind_max_r)/rnu)**0.50

      CALL MPI_ALLREDUCE (st_sum_0, st_sum_t, 1, MPI_REAL8, MPI_SUM, mpi_comm_world, ierr)
      CALL MPI_ALLREDUCE (st_max_0, st_max_t, 1, MPI_REAL8, MPI_SUM, mpi_comm_world, ierr)

      st_1       = (tau_CLOUD_const2 * radius(1)**2) / tauk_DNS 
      st_avg     = ( st_sum_t / real(npart_t) ) / st_1
      st_max_avg = ( st_max_t / real(nproc) ) / st_1

!      if ( id.eq.0 ) write (*,*) istep,npart,npart_t,ind_max_r,yp (7,ind_max_r),radius(1:2),disp_p(ind_max_r)
      if ( id.eq.0 ) write (1004, 518) istep, npart, npart_t, st_avg, st_max_avg

518   format ( 2x, I10, I10, I10, 2(1pe14.6) )

      RETURN
      END SUBROUTINE LOCAL_St_AVG

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
