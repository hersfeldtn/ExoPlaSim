! Test functions for haze radiative transfer

PROGRAM generate_data
	IMPLICIT NONE

	INTEGER,PARAMETER :: im=3
	INTEGER,PARAMETER :: jm=3
	INTEGER,PARAMETER :: nl=6
	INTEGER,PARAMETER :: ts=10
    INTEGER,PARAMETER :: NAERO=1
	INTEGER :: i,j,k,t


	REAL :: ta(im,jm,nl)
	REAL :: sigma(im,jm,nl)
	REAL :: ps(im,jm,nl)
	REAL :: rhog(im,jm,nl)
	REAL :: dh(im,jm,nl-1)
    REAL :: mmr(im,jm,nl,NAERO)
    REAL :: numrho(im,jm,nl,NAERO)

	DO k=1,nl
		DO i=1,im
			DO j=1,jm
				ta(i,j,k) = 300.-i-j-k
				ps(i,j,k) = 101325
				sigma(i,j,k) = 101325-i-j-10*k
                mmr(i,j,k,NAERO) = 3e-13-i*0.15e-13+j*0.3e-13-k*1.2e-13
			END DO
		END DO
	END DO
	
    call density(im,jm,nl,ta,ps,sigma,rhog)
    call hthick(im,jm,nl,sigma,ps,rhog,dh)
	
	WRITE(*,*) rhog,dh
	

END PROGRAM



    SUBROUTINE density(im,jm,nl,temp,ps,sigmah,   &
                        rhog)

!   Calculate the gas density in kg/m3 (input for terminal velocity calculation)
!****6***0*********0*********0*********0*********0*********0**********72
    
    IMPLICIT NONE
    
! Dimensions of arrays

    INTEGER,INTENT(IN) :: im,jm,nl ! Length of x, y, and z dimensions
    
! Define arrays

    REAL,INTENT(IN) :: temp(im,jm,nl) ! Temperature
    REAL,INTENT(IN) :: sigmah(nl) ! Sigma
    REAL,INTENT(IN) :: ps(im,jm) ! Surface air pressure
    REAL,INTENT(OUT) :: rhog(im,jm,nl) ! Density of surrounding gas

! Local arrays	
    REAL :: airp(im,jm,nl) ! Air pressure 
    INTEGER :: k ! For vertical loop
    
! Constants

    REAL :: R_gas ! Gas constant for ideal gas law
    REAL :: M_mass ! Molar mass of nitrogen
        
    R_gas = 8.314 ! SI units : J/K mol    
    M_mass = 0.0280 ! kg/mol for N2
    
! Calculations
    
    DO k=1,nl
        airp(:,:,k) = ps(:,:)*sigmah(k) ! Air pressure at mid level in Pa
    ENDDO
        
    rhog = (airp*M_mass)/(R_gas*temp)
    
    RETURN
    END
    
    
!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE hthick(im,jm,nl,sigma,ps,rhog,   &
                      dh)
!****6***0*********0*********0*********0*********0*********0**********72
    
    IMPLICIT NONE
    
    ! Dimensions of arrays

    INTEGER,INTENT(IN) :: im,jm,nl ! Length of x, y, and z dimensions
    
    ! Define arrays
    REAL,INTENT(IN) :: sigma(nl) ! Sigma at levels
    REAL,INTENT(IN) :: ps(im,jm) ! Surface air pressure
    REAL,INTENT(IN) :: rhog(im,jm,nl) ! Density of surrounding gas
    REAL,INTENT(OUT) :: dh(im,jm,nl-1) ! Height thickness
    
! Local arrays	
    REAL :: airp(im,jm,nl) ! Air pressure
    REAL :: ph(im,jm,nl-1) ! Pressure thickness
    
    INTEGER :: k ! For vertical loop
    
    REAL :: grav ! Planet's gravitational constant
    
    grav = 9.81 ! Earth
    ! Calculations
    
    DO k=1,nl
        airp(:,:,k) = ps(:,:)*sigma(k) ! Air pressure at mid level in Pa
    ENDDO
    
    DO k=1,nl-1
        ph(:,:,k) = airp(:,:,k) - airp(:,:,k+1)
    ENDDO
    
    dh = ph/(grav*rhog)
    
    RETURN
    END
    

!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE mmr2n(mmr,apart,rhop,rhog,NAERO,im,jm,nl,   &
                      numrho)
! Convert mass mixing ratio (kg/kg) from aerocore into number density (particles/m3)
!****6***0*********0*********0*********0*********0*********0**********72
    IMPLICIT NONE
    
    integer,intent(in) :: NAERO,im,jm,nl
    real,intent(in) ::   mmr(im,jm,nl,NAERO) ! Mixing ratio of aerosol
    real,intent(in)    :: apart ! Particle radius
    real,intent(in)    :: rhop ! Particle density
    real,intent(in)    :: rhog ! Bulk gas density
    real,intent(out) :: numrho(im,jm,nl,NAERO) ! Number density of aerosol
    
    REAL :: PI ! Value of pi
    PI = 3.14159
    rhop = 1262 ! kg/m3
    apart = 0.5e-06 ! m
    
    svol = (4/3)*PI*(apart**3) ! Sphere volume
    mpart = svol*rhop
    
    numrho = mmr*(1/mpart)*rhog
    RETURN
    END
    
    
!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE atrscat(im,jm,nl,numrho,qscat,apart,dz,   &
                      trscat)
! Convert mass mixing ratio (kg/kg) from aerocore into number density (particles/m3)
!****6***0*********0*********0*********0*********0*********0**********72
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: im,jm,nl ! Length of x, y, and z dimensions
    
    REAL,INTENT(IN) :: numrho(im,jm,nl)
    REAL,INTENT(IN) :: dz(im,jm,nl)
    REAL,INTENT(IN) :: qscat
    REAL,INTENT(IN) :: apart
    
    REAL,INTENT(OUT) :: trscat(im,jm,nl)
    
    qscat = 2.7881319182436437
    apart = 0.5e-06
    PI = 3.14159
    
    xsec = PI*(apart**2)*qscat
    trscat = exp(-numrho*dz)
    
    RETURN
    END
    