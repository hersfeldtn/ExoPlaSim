!
!     **********************************************************
!     *  Parameters and subroutines for aerosol transport       *
!     **********************************************************

      module aeromod

!     **********************************************************
!     * This module contains the parameters, arrays and        *
!     * subroutines that are needed for transporting aerosols   *
!     * using the Flux-Form Semi-Lagrangian (FFSL) algorithm   *
!     * developed by S.-J. Lin (now at GFDL).                  *
!     **********************************************************
!     * The original transport code (in F77) was written by    *
!     * S.-J. Lin. Adaptation for the Planet Simulator was     *
!     * done by Hui Wan (MPI-M).
!     **********************************************************

      use pumamod, only: NLAT,NLEV,nud

      logical,parameter :: aero_debug  = .TRUE.
      logical,parameter :: aero_zcross = .TRUE.
      logical,parameter :: aero_deform = .FALSE.

      logical,parameter :: aero_fill = .FALSE.
      logical,parameter :: aero_mfct = .FALSE.

      integer,parameter :: aero_iord = 2
      integer,parameter :: aero_jord = 2
      integer,parameter :: aero_kord = 3
	  
	  integer,parameter :: l_source = 1 ! 1 = photochemical haze (source at top level)
	                                  ! 2 = dust (source at bottom level)

      integer,parameter :: aero_cnst = 1   ! 1 = constant preserving
                                           ! 2 = mass conserving

!      integer,parameter :: aero_j1  = 2  ! 1st lat. outside polar cap
!      integer,parameter :: aero_j2  = NLAT + 1 - aero_j1 
                                         ! last lat. outside polar cap
										 
	  real,parameter :: apart = 50e-06 ! Radius of aerosol particle in meters
	  real,parameter :: rhop = 1000 ! Density of aerosol particle in kg/m3
	  real,parameter :: fcoeff = 10e-13 ! Haze particle mass production rate in kg/m2s

      end module aeromod


!     ======================
!     SUBROUTINE AERO_MAIN
!     ======================

      subroutine aero_main

      use pumamod, only: du,dv,dp,du0,dv0,dp0,daeros, &
                         NLON,NLAT,NLEV,NAERO,       &
                         mypid,NROOT,sigmah,dt
      use tracermod
	  use aeromod
      implicit none

      real :: zu   (NLON,NLAT,NLEV)
      real :: zv   (NLON,NLAT,NLEV)
      real :: zps0 (NLON,NLAT)
      real :: zps1 (NLON,NLAT)

      real :: x (NLON+1,NLAT,NLEV,NAERO)  ! for GUI output
      real :: y (NLON+1,NLAT,NLEV)         ! for GUI output

      integer :: j,jc

      character(len=9) :: aero_name

!     --- 

      call prepare_uvps( zu,zv,zps0,zps1,      & ! output
                         du0,dv0,dp0,du,dv,dp)   ! input

      if (mypid == NROOT .and. aero_debug) then
         write(nud,'(a,f11.2)') '* max aero u   =',maxval(abs(zu))
         write(nud,'(a,f11.2)') '* max v   =',maxval(abs(zv))
         write(nud,'(a,f11.2)') '* max ps0 =',maxval(zps0)
         write(nud,'(a,f11.2)') '* max ps1 =',maxval(zps1)
       ! write(nud,*)
       ! write(nud,*) 'ps0 NP'
       ! write(nud,*)  dp0(1:NLON)
       ! write(nud,*) 'zps0 NP'
       ! write(nud,*) zps0(1:NLON,NLAT)
       ! write(nud,*) 'ps0 SP'
       ! write(nud,*)  dp0(NLON*(NLAT-1)+1:)
       ! write(nud,*) 'zps0 SP'
       ! write(nud,*) zps0(1:NLON,1)
      end if
    
      if (mypid == NROOT) then

         call aerocore(daeros,l_source,sigmah,dt,         &
                      zps0,zps1,zu,zv,                    &
                      dtoa,dtdx,apart,rhop,fcoeff,        &
                      aero_iord,aero_jord,aero_kord,      & 
                      NAERO,NLON,NLAT,NLEV,dap,dbk,       &
                      iml,ffsl_j1,ffsl_j2,js0,jn0,        &
                      colae,colad,rcolad,dlat,rcap,       &
                      aero_cnst,aero_deform,aero_zcross,  &
                      aero_fill,aero_mfct,aero_debug,nud)

!        preparation for the GUI output: 
!        invert the meridional direction and add the 360 deg. longitude

         do j=1,NLAT
            x(1:NLON,j,:,:) = daeros(:,NLAT+1-j,:,:)
            x(NLON+1,j,:,:) = daeros(1,NLAT+1-j,:,:)
         end do

!        send all tracer fields to output

         do jc=1,NAERO
            write(aero_name,'(a,i2.2)') 'DAEROS',jc
            call guiput(aero_name // char(0), x(1,1,1,jc), NLON+1,NLAT,NLEV)
         enddo

!        check the correlation between tracers 3 and 4

!        y(:,:,:) = x(:,:,:,3)+x(:,:,:,4)
!        call guiput('TRC03+04' // char(0), y, NLON+1,NLAT,NLEV)

      end if

      return
      end subroutine aero_main
