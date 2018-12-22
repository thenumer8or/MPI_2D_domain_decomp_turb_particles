! ----------------------------------------------------------------
!       Spectral code for three-D isotropic turbulence using message-passing 
!       2D domain decomposition implemented by hparish and oayala. 
! 
!       The parallel FFT is described in the code mpifft.f90
!       It is used a 2D decomposition strategy with cyclic permitation
!       scheme for the communications among processors.
! ----------------------------------------------------------------- 
!       nx*2 = ny = nz = 2**n   ( n = 4,5,...)
!       nx = x-diretion mesh size in each node
!       ny = y-diretion mesh size in each node
!       nz = z-direction mesh size in each node
!       nproc = number of process (nodes)
!       lx=nx,ly=ny/nproxY,lz=nz/nprocZ   (local dimensions asigned to each node)
!       Total number of processors = nproxY * nprocZ
! ----------------------------------------------------------------- 

      PROGRAM TWOD_DNS_FLOW

      USE params
!      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'plans.h'
      INCLUDE 'fftw_include.h'

      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  vx, vy, vz
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  edissip, enstrophy
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  vxc, vyc, vzc
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  wx, wy, wz
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  ox, oy, oz
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  kx, ky, kz, k2
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  tmp
      INTEGER, ALLOCATABLE, DIMENSION (:,:,:)   ::  ik2
      REAL, ALLOCATABLE, DIMENSION    (:)       ::  enum
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  a1r,a2r,a3r
      REAL, ALLOCATABLE, DIMENSION    (:,:,:)   ::  b1r,b2r,b3r
      REAL, ALLOCATABLE, DIMENSION    (:)       ::  kxr,kyr,kzr
      REAL                                      ::  force(2)
      REAL                                      ::  co1,co2,co1_t,co2_t,vol
      INTEGER                                   ::  iseedf,iyf,ivf(32)
      INTEGER                                   ::  idp,nstep,forcing
      INTEGER                                   ::  i,j,k
      INTEGER                                   ::  id, ierr, requ(2), iseed_m0
      INTEGER                                   ::  resultlen
      INTEGER, DIMENSION (MPI_STATUS_SIZE)      ::  status
      LOGICAL                                   ::  new_flow
      REAL                                      ::  start_time,end_time,dtime,timeend,timestart,wtime
      REAL                                      ::  time_v(12)
      REAL, ALLOCATABLE, DIMENSION (:)          ::  time_full

      CHARACTER*(MPI_MAX_PROCESSOR_NAME) name
      CHARACTER*60 directvelo
      CHARACTER*60 directstat
      CHARACTER*60 velo_omp 
      CHARACTER*60 str
      CHARACTER*60 part_in
      CHARACTER*60 part_out
      CHARACTER*60 part_acc

      time_flow_cp      = 0.0
      time_flow_cm      = 0.0
      time_partadv_cp   = 0.0
      time_partadv_cm   = 0.0
      time_periodic_cp  = 0.0
      time_periodic_cm  = 0.0
      time_interp_cp    = 0.0
      time_interp_cm    = 0.0
      time_colldet_cp   = 0.0
      time_colldet_cm   = 0.0
      time_hdi_cp       = 0.0
      time_hdi_cm       = 0.0
      time_v            = 0.0   

!     Setup MPI environment
      CALL mpi_init (ierr)
      CALL mpi_comm_size (MPI_COMM_WORLD, nproc, ierr)
      CALL mpi_comm_rank (MPI_COMM_WORLD, id, ierr)
!      print *, "MYID = ", id

      CALL mpi_get_processor_name (name,resultlen,ierr)
!      print *, "Name = ", name

      start_time = mpi_wtime()

      CALL MPI_BCAST(start_time,1,MPI_REAL8,0,mpi_comm_world,ierr)

! ----------------------------------------------------------------- 
!     path to input, output files

      directstat = '/glade/scratch/hparish/grow_turb/stat/'
      directvelo = '/glade/scratch/hparish/grow_turb/velo/'
      velo_omp   = '/glade/scratch/hparish/grow_turb/velo_omp/'
      part_in    = '/glade/scratch/hparish/grow_turb/particle_in/'
      part_out   = '/glade/scratch/hparish/grow_turb/particle_out/'
      part_acc   = '/glade/scratch/hparish/grow_turb/particle_acc/'

! ----------------------------------------------------------------- 
!     id = 0 reads parameters from input file and writes them to file 

      IF ( id.eq.0 ) THEN
        open(90,file='parameter.d',status='unknown')
        read(90,*) nx, ny, nz
        read(90,*) str
        read(90,*) nprocY
        read(90,*) str
        read(90,*) nstep
        read(90,*) str
        read(90,*) dt
        read(90,*) str
        read(90,*) rnu
        read(90,*) str
        read(90,*) time
        read(90,*) str
        read(90,*) new_flow
        read(90,*) str
        read(90,*) idp
        read(90,*) str
        read(90,*) iseed_m0
        read(90,*) str
        read(90,*) forcing
        read(90,*) str
        read(90,*) force
        read(90,*) str
        read(90,*) iseedf
        read(90,*) str
        read(90,*) ieout
        read(90,*) str
        read(90,*) ispec
        read(90,*) str
        read(90,*) divstat
        read(90,*) str
        read(90,*) rho_water
        read(90,*) str
        read(90,*) rho_air
        read(90,*) str
        read(90,*) etk_DNS, tauk_DNS
        read(90,*) str
        read(90,*) st
        read(90,*) str
        read(90,*) ediss_CLOUD
        read(90,*) str
        read(90,*) viscosity_CLOUD
        read(90,*) str
        read(90,*) gravity
        read(90,*) str
        read(90,*) new_particle
        read(90,*) str
        read(90,*) nonoverlap
        read(90,*) str
        read(90,*) HDI_included
        read(90,*) str
        read(90,*) HDI_trunc
        read(90,*) str
        read(90,*) npart 
        read(90,*) str
        read(90,*) ndrag
        read(90,*) str
        read(90,*) nset
        read(90,*) str

        do il = 1, nset
         read(90,*) radius(il)
        enddo

        close(90)
      ENDIF 

      IF ( id.eq.0 ) THEN

        if ( HDI_included ) nonoverlap = .TRUE.

        open (20,file=TRIM(directstat)//'parameter.dat',status='unknown') 

        write(20,*) 'nx, ny, nz'
        write(20,*) nx, ny, nz
        write(20,*) 'Number of processors along Y-array and Z-array direction'
        write(20,*) nprocY,nproc/nprocY
        write(20,*) 'Number of time steps'
        write(20,*) nstep
        write(20,*) 'dt'
        write(20,*) dt
        write(20,*) 'rnu'
        write(20,*) rnu
        write(20,*) 'Initial time'
        write(20,*) time
        write(20,*) 'New flow ?'
        write(20,*) new_flow
        write(20,*) 'Seed for initial flow in case of deterministic forcing'
        write(20,*) iseed_m0
        write(20,*) 'Velocity start point'
        write(20,*) idp
        write(20,*) 'Forcing scheme: 1 deterministic, 2 stochastic'
        write(20,*) forcing
        write(20,*) 'Force for deterministic forcing scheme'
        write(20,*) force
        write(20,*) 'Seeding for stochastic forcing scheme'
        write(20,*) iseedf
        write(20,*) 'Save turbulence statistics every Y time steps, Y='
        write(20,*) ieout
        write(20,*) 'Save spectrum every X time steps, X='
        write(20,*) ispec
        write(20,*) 'Want divergence velocity and other statistics?'
        write(20,*) divstat
        write(20,*) 'water density [g/cm3]'
        write(20,*) rho_water
        write(20,*) 'air density [g/cm3]'
        write(20,*) rho_air
        write(20,*) 'DNS kolmogorov parameters: eta,time'
        write(20,*) etk_DNS, tauk_DNS
        write(20,*) 'shell thickness'
        write(20,*) st
        write(20,*) 'energy dissipation rate [cm2/s3]'
        write(20,*) ediss_CLOUD
        write(20,*) 'air viscosity [cm2/s]'
        write(20,*) viscosity_CLOUD
        write(20,*) 'gravity [cm/s2]'
        write(20,*) gravity
        write(20,*) 'New particles?'
        write(20,*) new_particle 
        write(20,*) 'nonoverlapping particles ?'
        write(20,*) nonoverlap
        write(20,*) 'HDI included ?'
        write(20,*) HDI_included
        write(20,*) 'HDI truncation radius ?'
        write(20,*) HDI_trunc
        write(20,*) 'Total number of particles?'
        write(20,*) npart 
        write(20,*) 'drag 0-Stokes 1-nonlinear'
        write(20,*) ndrag
        write(20,*) 'number of particle sets'
        write(20,*) nset
        write(20,*) 'particle radii in [um]'
        write(20,*) radius(1:nset)

        close(20)
      ENDIF 
    
      if ( id.eq.0 ) CALL WALLCLOCK_LIMIT (wtime)
!      if ( id.eq.0 ) wtime = 99999999.0

! ----------------------------------------------------------------- 
!     Broadcast parameters from id = 0 to rest of the world

      CALL MPI_BCAST (nx,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (ny,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (nz,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (nprocY,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (nstep,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (dt,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (rnu,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (time,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (new_flow,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (iseed_m0,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (idp,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (forcing,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (force,2,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (iseedf,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (ieout,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (ispec,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (divstat,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (wtime,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (rho_water,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (rho_air,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (tauk_DNS,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (etk_DNS,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (st,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (ediss_CLOUD,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (viscosity_CLOUD,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (gravity,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (new_particle,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (nonoverlap,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (HDI_included,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (HDI_trunc,1,MPI_REAL8,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (npart,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (ndrag,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (nset,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_BCAST (radius,nset,MPI_REAL8,0,mpi_comm_world,ierr)

      if ( id.eq.0 ) write(*,*) 'SAFE EXIT after', wtime-3.,'minutes'

      if((mod(id,nprocY)+1).eq.nprocY) then
        lyext=2
      else
        lyext=0
      endif

      nprocZ = nproc/nprocY
      lx  = nx
      llx = lx + 2
      ly  = ny/nprocY
      lly = ly + lyext
      lz  = nz/nprocZ
      idlocalY = id-int(id/nprocY)*nprocY
      idlocalZ = int(id/nprocY)
      pi  = 4.*atan(1.)
      pi2 = pi*2.0

      hx  = pi2
      hy  = REAL ( pi2 / nprocY)
      hz  = REAL ( pi2 / nprocZ)

      akmax = nx/2 - 1.5
      nek   = int(akmax)

! ----------------------------------------------------------------- 
!     Create FFTW plans, store in plan1DF, plan1DB, plan2DRC, plan2DCR...
!     (passed through common in included file plans.h).

      CALL set_plan1D_RC(nx)
      CALL set_plan1D_CC(ny)

      ALLOCATE ( vxc(llx,lly,lz)  )
      ALLOCATE ( vyc(llx,lly,lz)  )
      ALLOCATE ( vzc(llx,lly,lz)  )
      ALLOCATE (  vx(llx,lly,lz)  )
      ALLOCATE (  vy(llx,lly,lz)  )
      ALLOCATE (  vz(llx,lly,lz)  )
      ALLOCATE (  wx(llx,lly,lz)  )
      ALLOCATE (  wy(llx,lly,lz)  )
      ALLOCATE (  wz(llx,lly,lz)  )
      ALLOCATE (  ox(llx,lly,lz)  )
      ALLOCATE (  oy(llx,lly,lz)  )
      ALLOCATE (  oz(llx,lly,lz)  )
      ALLOCATE (  kx(llx,lly,lz)  )
      ALLOCATE (  ky(llx,lly,lz)  )
      ALLOCATE (  kz(llx,lly,lz)  )
      ALLOCATE (  k2(llx,lly,lz)  )
      ALLOCATE ( ik2(llx,lly,lz)  )
      ALLOCATE ( tmp(llx,lly,lz)  )
      ALLOCATE ( edissip(llx,lly,lz)   )
      ALLOCATE ( enstrophy(llx,lly,lz) )
      ALLOCATE ( enum(nek)        )
      ALLOCATE ( time_full(15*nproc) )
!      ALLOCATE ( removepart(m2send)  )

      if(forcing.eq.2) then
        ALLOCATE ( kxr(nx+2) )
        ALLOCATE ( kyr(ny) )
        ALLOCATE ( kzr(nz) )
        ALLOCATE ( a1r(6,5,5) )
        ALLOCATE ( a2r(6,5,5) )
        ALLOCATE ( a3r(6,5,5) )
        ALLOCATE ( b1r(6,5,5) )
        ALLOCATE ( b2r(6,5,5) )
        ALLOCATE ( b3r(6,5,5) )
      endif

      if ( id.eq.0 ) print*, ' memory allocated'

      open(110 ,file=TRIM(directstat)//'Vp11.dat',position='append',status='unknown')
      open(120 ,file=TRIM(directstat)//'Vp12.dat',position='append',status='unknown')
      open(220 ,file=TRIM(directstat)//'Vp22.dat',position='append',status='unknown')

      IF ( id.eq.0 ) THEN
         open(20  ,file=TRIM(directstat)//'spectrum.dat',status='unknown') 
         open(22  ,file=TRIM(directstat)//'transfer.dat',status='unknown') 
         open(40  ,file=TRIM(directstat)//'monitora.dat',status='unknown')
         open(41  ,file=TRIM(directstat)//'monitorb.dat',status='unknown')
         open(42  ,file=TRIM(directstat)//'monitorc.dat',status='unknown')
!         open(43  ,file=TRIM(directstat)//'DuDT.dat',status='unknown')
!         open(44  ,file=TRIM(directstat)//'DvDT.dat',status='unknown')
!         open(45  ,file=TRIM(directstat)//'DwDT.dat',status='unknown')
         open(50  ,file=TRIM(directstat)//'collision.dat',status='unknown')
         open(51  ,file=TRIM(directstat)//'detectpairs.dat',status='unknown')
!         open(110 ,file=TRIM(directstat)//'Vp11.dat',status='unknown')
!         open(120 ,file=TRIM(directstat)//'Vp12.dat',status='unknown')
!         open(220 ,file=TRIM(directstat)//'Vp22.dat',status='unknown')
         open(106 ,file=TRIM(directstat)//'meanv_part_vel.dat',status='unknown')
         open(107 ,file=TRIM(directstat)//'vel_mag_part.dat',status='unknown')
!         open(61  ,file=TRIM(directstat)//'concentration.dat',status='unknown')
         open(1011,file=TRIM(directstat)//'nirbin11.dat',status='unknown')
         open(1012,file=TRIM(directstat)//'nirbin12.dat',status='unknown')
         open(1022,file=TRIM(directstat)//'nirbin22.dat',status='unknown')
         open(5110,file=TRIM(directstat)//'wxbin11.dat',status='unknown')
         open(5120,file=TRIM(directstat)//'wxbin12.dat',status='unknown')
         open(5220,file=TRIM(directstat)//'wxbin22.dat',status='unknown')
         open(6110,file=TRIM(directstat)//'wybin11.dat',status='unknown')
         open(6120,file=TRIM(directstat)//'wybin12.dat',status='unknown')
         open(6220,file=TRIM(directstat)//'wybin22.dat',status='unknown')
         open(7110,file=TRIM(directstat)//'wzbin11.dat',status='unknown')
         open(7120,file=TRIM(directstat)//'wzbin12.dat',status='unknown')
         open(7220,file=TRIM(directstat)//'wzbin22.dat',status='unknown')
         open(1110,file=TRIM(directstat)//'wrbin11.dat',status='unknown')
         open(1120,file=TRIM(directstat)//'wrbin12.dat',status='unknown')
         open(1220,file=TRIM(directstat)//'wrbin22.dat',status='unknown')
         open(2110,file=TRIM(directstat)//'thetabin11.dat',status='unknown')
         open(2120,file=TRIM(directstat)//'thetabin12.dat',status='unknown')
         open(2220,file=TRIM(directstat)//'thetabin22.dat',status='unknown')
         open(3110,file=TRIM(directstat)//'urbin11.dat',status='unknown')
         open(3120,file=TRIM(directstat)//'urbin12.dat',status='unknown')
         open(3220,file=TRIM(directstat)//'urbin22.dat',status='unknown')
         open(4110,file=TRIM(directstat)//'wrtanbin11.dat',status='unknown')
         open(4120,file=TRIM(directstat)//'wrtanbin12.dat',status='unknown')
         open(4220,file=TRIM(directstat)//'wrtanbin22.dat',status='unknown')
         open(88  ,file=TRIM(directstat)//'timing.dat',status='unknown')
         open(1001,file=TRIM(directstat)//'size_dist.dat',status='unknown')
         open(1002,file=TRIM(directstat)//'mdens_dist.dat',status='unknown')
         open(1003,file=TRIM(directstat)//'size_avg.dat',status='unknown')
         open(1004,file=TRIM(directstat)//'st_avg.dat',status='unknown')
      END IF

      CALL wavenumber (kx,ky,kz,k2,id)

      ik2 = int ( sqrt(k2) + 0.5 )
      ik2(lx+1:llx,:,:) = 0

      DO i = 1, nek
        co1 = float ( count(ik2(1:lx,:,:).eq.i.and.kx(1:lx,:,:).gt.0.1) )
        co2 = float ( count(ik2(1:lx,:,:).eq.i.and.kx(1:lx,:,:).lt.0.1) )
        CALL MPI_ALLREDUCE (co1,co1_t,1,MPI_REAL8,MPI_SUM,&
                                     mpi_comm_world,ierr)
        CALL MPI_ALLREDUCE (co2,co2_t,1,MPI_REAL8,MPI_SUM,&
                                     mpi_comm_world,ierr)
        enum(i) = co1_t + co2_t/2.
        vol = 4./3. * pi * ( (float(i)+0.5)**3 - (float(i)-0.5)**3 )
        if ( id.eq.0 ) write(*,*) i,enum(i)/vol
      ENDDO

      if(forcing.eq.2) CALL SETWAVERAND (kxr,kyr,kzr)

! ------- initiate particles---------------------------------

      IF ( npart .gt. 0 ) then

       CALL PART_INIT (id, part_in)
       CALL PERIODIC (id)

       IF ( nonoverlap ) CALL OVERLAP_DET (id)

      ENDIF

      stime = mpi_wtime()

! =================================================================  
!     Generate initial condition

      IF ( new_flow ) THEN
        IF (forcing.eq.1) THEN
         CALL initial_flow(vx,vy,vz,iseed_m0,ik2,enum,kx,ky,kz,k2,force,id)
        ELSEIF (forcing.eq.2) THEN
          iseedf = -iseedf
          ivf(:) = 0.0
          iyf = 0
          b1r = 0.0
          b2r = 0.0
          b3r = 0.0
          vx = 0.0
          vy = 0.0
          vz = 0.0
        ENDIF
        ox = vx
        oy = vy
        oz = vz
      ELSE
! ---- Reading from DISK for the initial flow ----------------------

        CALL input_output (directvelo,vx,vy,vz,ox,oy,oz,idp,id,1)
        if(forcing.eq.2) CALL input_outputf (directvelo,b1r,b2r,b3r,iseedf,ivf,iyf,idp,id,1)

        if(forcing.eq.2) then
         CALL MPI_BCAST (iseedf,1,MPI_INTEGER,0,mpi_comm_world,ierr)
         CALL MPI_BCAST (ivf,NTAB,MPI_INTEGER,0,mpi_comm_world,ierr)
         CALL MPI_BCAST (iyf,1,MPI_INTEGER,0,mpi_comm_world,ierr)
         CALL MPI_BCAST (b1r,150,MPI_REAL8,0,mpi_comm_world,ierr)
         CALL MPI_BCAST (b2r,150,MPI_REAL8,0,mpi_comm_world,ierr)
         CALL MPI_BCAST (b3r,150,MPI_REAL8,0,mpi_comm_world,ierr)
        endif

        if ( id .eq. 0 ) write(*,*) 'flow field & stoch. forcing parametres are read from previos MPI_2d run' 

      ENDIF

! --- Save data in Open MP format ----------------------
!      CALL SAVEOMP (vx,vy,vz,ox,oy,oz,velo_omp,id)
!      if(forcing.eq.2) CALL input_outputf (velo_omp,b1r,b2r,b3r,iseedf,ivf,iyf,idp+1,id,2)

!     CALL MPI_FINALIZE (ierr)
!     STOP

! =================================================================

      dtr  = 1.5 * dt
      dt_h = 0.5 * dt
      istep = 0

      etime        = mpi_wtime()
      time_flow_cp = time_flow_cp + etime - stime

      timestart = mpi_wtime()

1     istep = istep + 1
      jstep = istep

      stime = mpi_wtime()

! -------Postprocessing section---------------------------------
! Write out Du/DT to learn its importance
!      if (istep.gt.1) then
!        CALL DUDTPOSTPROC(kx,ky,kz,vx,vy,vz,vxc,1,id)
!        CALL DUDTPOSTPROC(kx,ky,kz,vx,vy,vz,vyc,2,id)
!        CALL DUDTPOSTPROC(kx,ky,kz,vx,vy,vz,vzc,3,id)
!      endif
! ------END Postprocessing section------------------------------

      wx = ky*vz - kz*vy
      wy = kz*vx - kx*vz
      wz = kx*vy - ky*vx
      tmp = wx
      wx(1:lx , 1:lly-1:2 , :) = -tmp(1:lx , 2:lly:2   , :)
      wx(1:lx , 2:lly:2   , :) =  tmp(1:lx , 1:lly-1:2 , :)
      tmp = wy
      wy(1:lx , 1:lly-1:2 , :) = -tmp(1:lx , 2:lly:2   , :)
      wy(1:lx , 2:lly:2   , :) =  tmp(1:lx , 1:lly-1:2 , :)
      tmp = wz
      wz(1:lx , 1:lly-1:2 , :) = -tmp(1:lx , 2:lly:2   , :)
      wz(1:lx , 2:lly:2   , :) =  tmp(1:lx , 1:lly-1:2 , :)

      vxc = vx
      vyc = vy
      vzc = vz

!     Parallel FFT, vorticity to physical space
      CALL mpifft3DCR (wx,id)
      CALL mpifft3DCR (wy,id)
      CALL mpifft3DCR (wz,id)

      nstep_grow = 40000 

      IF ((istep .ge. nstep_grow) .and. (mod (istep,100) .eq. 0)) THEN 
!       CALL COMP_ENSTROPHY(wx,wy,wz,enstrophy)
       CALL DISSIP(vx,vy,vz,kx,ky,kz,rnu,edissip,id)
      ENDIF

!     Parallel FFT, velocity to physical space

      CALL mpifft3DCR (vx,id)
      CALL mpifft3DCR (vy,id)
      CALL mpifft3DCR (vz,id)

      etime        = mpi_wtime()
      time_flow_cp = time_flow_cp + etime - stime 

! --------------- particle part: interpolate velocity for particle locations -------------------

      if ( npart .gt. 0 ) then

!      if no turbulence, this is a funny way of turning the turbulence off for particles.
!      vx = 0
!      vy = 0
!      vz = 0

       CALL PREPARE_INTERP (id)
       CALL INTERPOLATE (vx, 1, id)
       CALL INTERPOLATE (vy, 2, id)
       CALL INTERPOLATE (vz, 3, id)

       IF ((istep .ge. nstep_grow) .and. (mod (istep,100) .eq. 0)) THEN 
!        CALL INTERP_SCALAR (enstrophy, 1, id)
        CALL INTERP_SCALAR (edissip, 2, id)
        CALL LOCAL_St_AVG (id)
       endif

       CALL PART_HISTORY_FVP (id)

       dist_headlist = d_coll
       if ( HDI_included ) dist_headlist = d_hdi

       CALL HEADLIST (id, dist_headlist)

       if ( HDI_included ) CALL NEI_LIST_PRECOMP (id)
       if ( HDI_included ) CALL HDI_GMRES (id)

       CALL PART_HISTORY_FDP (id)

       CALL PART_ADVANCE (id)

!       IF ( (jstep .gt. (nstep/2)) .and. (mod(jstep,ieout) .eq. 0) ) CALL PART_ACCELERATION (id, idp+1, part_acc)

       CALL COMPLETE_YPN (id)

       IF ( mod(jstep,ieout) .eq. 0 ) CALL MEAN_VEL (id)

       CALL COLLISION_DET (id)

       CALL UPDATE_PART_POSITION (id)

!       nstep_grow = 18000 

       IF ( (istep .lt. nstep_grow) .and. ( nonoverlap .and. (ncoll .gt. 0) ) ) CALL RELOCATE (id)
       IF ( (istep .ge. nstep_grow) .and. ( nonoverlap .and. (ncoll .gt. 0) ) ) CALL COALESCE (id)

       CALL MPI_BARRIER (MPI_COMM_WORLD, error)
       CALL MPI_ALLREDUCE (nps, npart_t, 1, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)

!       IF ((istep .ge. nstep_grow) .and. (mod (istep,100) .eq. 0)) CALL SIZE_DIST_LINEAR (id)
       IF ((istep .ge. nstep_grow) .and. (mod (istep,100) .eq. 0)) CALL SIZE_DIST_LOGARITHMIC (id)
       IF ((istep .ge. nstep_grow) .and. (mod (istep,100) .eq. 0)) CALL SIZE_AVG (id)

       IF ( npart_t .lt. (npart/2) ) THEN
        if ( id.eq.0 ) write (*,*) 'npart_t < npart/2 = ', npart_t, npart
        CALL mpi_finalize(ierr)
        STOP
       ENDIF

       if (istep .ge. nstep_grow) CALL HEADLIST (id, dist_headlist)

       CALL PERIODIC (id)

       IF ( nonoverlap ) CALL OVERLAP_DET (id)

117    format ( 3I6, 8(1pe18.10) )

! ------------------------------------------- particle part: end ----------------------------------

!       open (178786,file='fort.178786',status='unknown',position='append')
!       do ih = 1, nps
!          partind = int ( yp(8,ih) )
!          if ( partind .eq. 178786 ) write(178786,117) id, jstep, ih, yp(:,ih)
!       enddo
!       close (178786)

      endif
       
      stime = mpi_wtime()

! -------Postprocessing section---------------------------------
      if(mod(jstep,ieout).eq.0) then
         CALL POSTPROC(vxc,vyc,vzc,vx,vy,vz,wx,wy,wz, &
                        kx,ky,kz,k2,ik2,rnu,id)
      endif
! ------END Postprocessing section------------------------------

      tmp = vy*wz - vz*wy
      wy  = vx*wy - vy*wx
      wz  = vz*wx - vx*wz
      wx  = tmp

      vx  = vxc
      vy  = vyc
      vz  = vzc

!     Parallel FFT, (u x w) to spectral space
      CALL mpifft3DRC (wx,id)
      CALL mpifft3DRC (wz,id)
      CALL mpifft3DRC (wy,id)

      CALL dealiasing (wx,ik2)
      CALL dealiasing (wy,ik2)
      CALL dealiasing (wz,ik2)

! -------Postprocessing section---------------------------------
      if(mod(jstep,ispec).eq.0) then
!         CALL TRANSPOSTPROC(vx,vy,vz,wx,wz,wy,k2,id)
         CALL TRANSPOSTPROC(vx,vy,vz,wx,wz,wy,ik2,id)
      endif
! ------END Postprocessing section------------------------------

! --------------------------------------------------------------
!     Crank-Nicholson integration scheme

      vx = vx * (1.-rnu*dt_h*k2) + dtr*wx - dt_h*ox
      vy = vy * (1.-rnu*dt_h*k2) + dtr*wz - dt_h*oy
      vz = vz * (1.-rnu*dt_h*k2) + dtr*wy - dt_h*oz

      if(forcing.eq.2) CALL supforSTH (vx,vy,vz,iseedf,ivf,iyf,  &
                            kxr,kyr,kzr,a1r,a2r,a3r,b1r,b2r,b3r,id)

      tmp  = ( kx*vx + ky*vy + kz*vz ) / k2
      vx = vx - kx*tmp
      vy = vy - ky*tmp
      vz = vz - kz*tmp

      vx = vx / (1.+rnu*dt_h*k2)
      vy = vy / (1.+rnu*dt_h*k2)
      vz = vz / (1.+rnu*dt_h*k2)

      if(forcing.eq.1) CALL supforDET (vx,vy,vz,k2,ik2,force,id)


      CALL symmetrize (vx,id)
      CALL symmetrize (vy,id)
      CALL symmetrize (vz,id)

      CALL dealiasing (vx,ik2)
      CALL dealiasing (vy,ik2)
      CALL dealiasing (vz,ik2)

      ox = wx
      oy = wz
      oz = wy

      time  = time + dt

      end_time = mpi_wtime()
      dtime = ( end_time - start_time ) / 60.
      CALL MPI_BCAST (dtime,1,MPI_REAL8,0,mpi_comm_world,ierr)

      if ( id.eq.0 ) write(*,*) istep

!-------------------------------------------------------------------
!     for timing

      etime        = mpi_wtime()
      time_flow_cp = time_flow_cp + etime - stime 

      if ( (jstep.lt.nstep) .and. ((wtime-dtime).gt.3.) ) goto 1

      timeend = mpi_wtime()
      if ( ((wtime-dtime).lt.3.).and.(id.eq.0) ) write(*,*) 'SAFE EXIT'

      if (id.eq.0) write(*,*) 'total time              =', timeend-timestart
      if (id.eq.0) write(*,*) 'Next run starts at time =', time

! --- Save final velocity field & forcing data   ----------------------
      CALL input_output (directvelo,vx,vy,vz,ox,oy,oz,idp+1,id,2)
      if(forcing.eq.2) CALL input_outputf (directvelo,b1r,b2r,b3r,iseedf,ivf,iyf,idp+1,id,2)

! --- Save final velocity field in Open MP format ----------------------
!      CALL SAVEOMP (vx,vy,vz,ox,oy,oz,velo_omp,id)
!      if(forcing.eq.2) CALL input_outputf (velo_omp,b1r,b2r,b3r,iseedf,ivf,iyf,idp+1,id,2)

      if ( npart.gt.0 ) CALL SAVE_PARTICLE (id,part_out)
!      if ( npart.gt.0 ) CALL SAVE_PART_GLOBAL (id,part_out)                                                                                                                               
!       index_yp = 463 + id
!       write (index_yp,6060) yp
!6060   format ( 2x,7(1pe25.15) )
                                                                                                                                                                 
!      CALL mpi_finalize(ierr)
!      STOP

      if ( id.eq.0 ) write(*,*) 'Successfully completed'

      time_v(1)  = real ( time_flow_cp     )
      time_v(2)  = real ( time_flow_cm     )
      time_v(3)  = real ( time_partadv_cp  )
      time_v(4)  = real ( time_partadv_cm  )
      time_v(5)  = real ( time_periodic_cp )
      time_v(6)  = real ( time_periodic_cm )
      time_v(7)  = real ( time_interp_cp   )
      time_v(8)  = real ( time_interp_cm   )
      time_v(9)  = real ( time_colldet_cp  ) 
      time_v(10) = real ( time_colldet_cm  ) 
      time_v(11) = real ( time_hdi_cp      ) 
      time_v(12) = real ( time_hdi_cm      ) 

      CALL MPI_GATHER (time_v,12,MPI_REAL8,time_full,12,MPI_REAL8,0,mpi_comm_world,err)
                                                                                                                                           
      if ( id.eq.0 ) then
        do idm = 0, nproc-1
          write(88,505) idm, time_full((1+12*idm):(12+12*idm))
        enddo
      endif
       
505   format ( 2x, I4, 12(1pe10.2) )

      if ( id.eq.0 ) then
        close(20)
        close(22)
        close(40)
        close(41)
        close(42)
!        close(43)
!        close(44)
!        close(45)
        close(50)
        close(51)
        close(88)
        close(110)
        close(120)
        close(220)
        close(106)
        close(107)
!        close(61)
        close(1011)
        close(1012)
        close(1022)
        close(1110)
        close(1120)
        close(1220)
        close(2110)
        close(2120)
        close(2220)
        close(3110)
        close(3120)
        close(3220)
        close(1001)
        close(1002)
        close(1003)
        close(1004)
      end if

      DEALLOCATE ( vx, vy, vz, wx, wy, wz    )
      DEALLOCATE ( vxc, vyc, vzc, ox, oy, oz )
      DEALLOCATE ( kx, ky, kz, k2, ik2       )
      DEALLOCATE ( tmp, enum                 )
      DEALLOCATE ( time_full                 )
      DEALLOCATE ( edissip, enstrophy        )

      if ( forcing .eq. 2 ) then
        DEALLOCATE ( kxr, kyr, kzr )
        DEALLOCATE ( a1r, a2r, a3r )
        DEALLOCATE ( b1r, b2r, b3r )
      endif

      if ( npart .gt. 0) then

        DEALLOCATE ( ypn, vpn       )
        DEALLOCATE ( pindex         )
        DEALLOCATE ( removepart     )
        DEALLOCATE ( yp, head, list )
        DEALLOCATE ( wxbin, wybin, wzbin )
        DEALLOCATE ( nirbin, wrbin       )
        DEALLOCATE ( urbin, thetabin     )
        DEALLOCATE ( wrtanbin            )
        DEALLOCATE ( collist             )
        DEALLOCATE ( enst_p, disp_p      )

        if ( HDI_included ) then
         DEALLOCATE ( pertvel, Axreshaped, breshaped )
!        DEALLOCATE ( nei_num, nei_data              )
        endif

      endif

      CALL destroy_plan1D_RC
      CALL destroy_plan1D_CC

      CALL mpi_finalize(ierr)
      STOP
      END PROGRAM TWOD_DNS_FLOW
