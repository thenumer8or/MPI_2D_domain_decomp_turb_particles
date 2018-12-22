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
         CALL MPI_ISEND (temp_p(:,1:ih_p(jid),jid),18*ih_p(jid),MPI_REAL8,ids,213,MPI_COMM_WORLD,req(1),ierr)
        endif
                                                                                                                                                                     
        if ( ih_pt(jid) .gt. 0 ) then
         CALL MPI_IRECV (temp_b(:,1:ih_pt(jid),jid),18*ih_pt(jid),MPI_REAL8,idr,213,MPI_COMM_WORLD,req(2),ierr)
        endif
                                                                                                                                                                     
        if ( ih_p(jid).gt.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAITALL (2,req,status2,ierr)
        if ( ih_p(jid).gt.0 .and. ih_pt(jid).eq.0 ) CALL MPI_WAIT (req(1),status1,ierr)
        if ( ih_p(jid).eq.0 .and. ih_pt(jid).gt.0 ) CALL MPI_WAIT (req(2),status1,ierr)
                                                                                                                                                                     
        if ( ih_pt(jid) .gt. 0 ) then
                                                                                                                                                                     
          nps1 = nps + SUM ( ih_pt ( 1 : (jid-1) ) )
          nps2 = nps + SUM ( ih_pt ( 1 : jid ) ) 

          ypn (1:15, (nps1+1):nps2) = temp_b(1:15,1:ih_pt(jid),jid)
          vpn (1:3, (nps1+1):nps2)  = temp_b(16:18,1:ih_pt(jid),jid)
                                                                                                                                                                     
         endif
                                                                                                                                                                     
        ENDDO

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
