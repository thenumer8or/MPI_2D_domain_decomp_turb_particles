       MODULE PARAMS

       PARAMETER            (NTAB=32)

       LOGICAL           :: divstat, new_particle, nonoverlap, HDI_included

       INTEGER           :: nx,ny,nz,lx,ly,lz,lyext,llx,lly
       INTEGER           :: nproc,nprocY,nprocZ,idlocalY,idlocalZ,nbin

       INTEGER           :: iseedr, iyr, ivr(NTAB), idum2
       INTEGER           :: istep, jstep, ieout, ispec, ncdx, ncdy, ncdz, nek
       INTEGER           :: npart, nps, nset, ndrag, msize, ncoll, ih_p(8), ih_pt(8), nei_id(8)
       INTEGER           :: ntotpair, npart_t
       INTEGER           :: nx_coll, ny_coll, nz_coll, nx_hdi, ny_hdi, nz_hdi
       INTEGER           :: nx_max, ny_max, nz_max, indy, indz, nuniform, ncollmax
       INTEGER           :: id_e, id_ne, id_n, id_nw, id_w, id_sw, id_s, id_se

!      INTEGER           :: irr,idp,iii,kkk,jjj
!      INTEGER           :: isign,ieout,itout,nek,istep

       REAL              :: dt, dt_h, dtr, dtc, dtau, time, d_coll, d_hdi, d_overl, HDI_trunc
       REAL              :: pi, pi2, hx, hy, hz, radius(10), radmax, dxyz(3)
       REAL              :: rnu, rho_water, rho_air, gravity, gravity_dns, ucharac, st, akmax 
       REAL              :: etk_CLOUD, viscosity_CLOUD, ediss_CLOUD, tauk_CLOUD, vk_CLOUD
       REAL              :: vk_DNS, etk_DNS, tauk_DNS, vk_ratio, etk_ratio, tauk_ratio
       REAL              :: tau_CLOUD_const, Sv_const, tau_CLOUD_const2, NORM

       DOUBLE PRECISION  :: etime, stime
       DOUBLE PRECISION  :: time_flow_cp, time_flow_cm
       DOUBLE PRECISION  :: time_interp_cp, time_interp_cm
       DOUBLE PRECISION  :: time_partadv_cp, time_partadv_cm
       DOUBLE PRECISION  :: time_periodic_cp, time_periodic_cm
       DOUBLE PRECISION  :: time_colldet_cp, time_colldet_cm
       DOUBLE PRECISION  :: time_hdi_cp, time_hdi_cm

       INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: head
       INTEGER, ALLOCATABLE, DIMENSION (:)     :: list
       INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: lhnode
       INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: pindex 
       INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: removepart
       INTEGER, ALLOCATABLE, DIMENSION (:)     :: ncollpart, npair_t 
       INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: nirbin
       INTEGER, ALLOCATABLE, DIMENSION (:)     :: nei_num 

!      INTEGER, ALLOCATABLE, DIMENSION (:)     :: lindex

       REAL, ALLOCATABLE, DIMENSION (:)        :: disp_p, enst_p
       REAL, ALLOCATABLE, DIMENSION (:,:)      :: yp, vp, dmove
       REAL, ALLOCATABLE, DIMENSION (:,:)      :: ypn, vpn, collist
       REAL, ALLOCATABLE, DIMENSION (:,:,:)    :: bg
       REAL, ALLOCATABLE, DIMENSION (:,:,:)    :: fvp, fdp
       REAL, ALLOCATABLE, DIMENSION (:,:,:)    :: wrbin, wxbin, wybin, wzbin
       REAL, ALLOCATABLE, DIMENSION (:,:,:)    :: urbin, thetabin, wrtanbin
       REAL, ALLOCATABLE, DIMENSION (:,:)      :: pertvel, breshaped, Axreshaped
       REAL, ALLOCATABLE, DIMENSION (:,:,:)    :: nei_data 

       END MODULE PARAMS

!====================================================

