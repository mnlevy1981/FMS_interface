module constants_mod

use shr_kind_mod,   only : R8 => shr_kind_r8
use shr_const_mod,  only : SHR_CONST_PI, &
     SHR_CONST_CDAY, &
     SHR_CONST_OMEGA, &
     SHR_CONST_REARTH, &
     SHR_CONST_G, &
     SHR_CONST_STEBOL,&             !   = 5.67e-8_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
     SHR_CONST_AVOGAD,&             !   = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
     SHR_CONST_MWDAIR,&             !   = 28.966_R8       ! molecular weight dry air ~ kg/kmole
     SHR_CONST_MWWV,&               !   = 18.016_R8       ! molecular weight water vapor
     SHR_CONST_RDAIR,&              ! Dry air gas constant     ~ J/K/kg
     SHR_CONST_RWV,&                ! Water vapor gas constant ~ J/K/kg
     SHR_CONST_KARMAN,&             !   = 0.4_R8          ! Von Karman constant
     SHR_CONST_PSTD,&               !   = 101325.0_R8     ! standard pressure ~ pascals
     SHR_CONST_TKFRZ,&              !   = 273.15_R8       ! freezing T of fresh water          ~ K
     SHR_CONST_RHODAIR,&            ! density of dry air at STP  ~ kg/m^3
     SHR_CONST_RHOFW,&              !   = 1.000e3_R8      ! density of fresh water     ~ kg/m^3
     SHR_CONST_RHOSW,&              !   = 1.026e3_R8      ! density of sea water       ~ kg/m^3
     SHR_CONST_CPDAIR,&             !   = 1.00464e3_R8    ! specific heat of dry air   ~ J/kg/K
     SHR_CONST_CPWV,&               !   = 1.810e3_R8      ! specific heat of water vap ~ J/kg/K
     SHR_CONST_CPSW,&
     SHR_CONST_LATICE,&             !   = 3.337e5_R8      ! latent heat of fusion      ~ J/kg
     SHR_CONST_LATVAP               !   = 2.501e6_R8      ! latent heat of evaporation ~ J/kg

implicit none
private

! Include variable "version" to be written to log file.
#include<file_version.h>
!-----------------------------------------------------------------------
! version is public so that write_version_number can be called for constants_mod
! by fms_init
public :: version

real :: realnumber !< dummy variable to use in HUGE initializations

!-----------------------------------------------------------------------
! The following variables are overridden by CESM values in order to keep
! dynamics and physics values consistent and allow a more accurate comparison
! of energy between dynamics and physics.

! Below are the original GFDL FMS values and the CESM value

!! gfdl_val real,public, parameter :: RADIUS = 6371.0e+3_r8_kind * small_fac   !< Radius of the Earth [m]
real(R8),        public, parameter :: RADIUS = SHR_CONST_REARTH                         !< Radius of the Earth [m]

!! gfdl_val real(kind=8), public, parameter :: PI_8   = 3.14159265358979323846_r8_kind  !< Ratio of circle circumference to diameter [N/A]
real(R8),        public, parameter :: PI_8   = SHR_CONST_PI                             !< Ratio of circle circumference to diameter [N/A]

!! gfdl_val real, public, parameter :: PI     = 3.14159265358979323846_r8_kind  !< Ratio of circle circumference to diameter [N/A]
real(R8),         public, parameter :: PI     = SHR_CONST_PI                            !< Ratio of circle circumference to diameter [N/A]

!! gfdl_val real, public, parameter :: OMEGA  = 7.292e-5_r8_kind / small_fac    !< Rotation rate of the Earth [1/s]
real(R8),         public, parameter :: OMEGA  = SHR_CONST_OMEGA                         !< Rotation rate of the Earth [1/s]

!! gfdl_val real, public, parameter :: GRAV   = 9.80_r8_kind             !< Acceleration due to gravity [m/s^2]
real(R8),         public, parameter :: GRAV   = SHR_CONST_G              !< Acceleration due to gravity [m/s^2]

!! gfdl_val real, public, parameter :: RDGAS  = 287.04_r8_kind           !< Gas constant for dry air [J/kg/deg]
real(R8),         public, parameter :: RDGAS  = SHR_CONST_RDAIR          !< Gas constant for dry air [J/kg/deg]

!! gfdl_val real, public, parameter :: RVGAS  = 461.50_r8_kind           !< Gas constant for water vapor [J/kg/deg]
real(R8),         public, parameter :: RVGAS  = SHR_CONST_RWV            !< Gas constant for water vapor [J/kg/deg]

!! gfdl_val real, public, parameter :: HLV = 2.500e6_r8_kind             !< Latent heat of evaporation [J/kg]
real(R8),         public, parameter :: HLV = SHR_CONST_LATVAP            !< Latent heat of evaporation [J/kg]

!! gfdl_val real, public, parameter :: HLF = 3.34e5_r8_kind              !< Latent heat of fusion [J/kg]
real(R8),         public, parameter :: HLF = SHR_CONST_LATICE            !< Latent heat of fusion [J/kg]

!! gfdl_val real, public, parameter :: CP_AIR = RDGAS/KAPPA              !< Specific heat capacity of dry air at constant pressure [J/kg/deg]
real(R8),         public, parameter :: CP_AIR = SHR_CONST_CPDAIR         !< Specific heat capacity of dry air at constant pressure [J/kg/deg]

!! gfdl_val real, public, parameter :: KAPPA  = 2.0_r8_kind/7.0_r8_kind  !< RDGAS / CP_AIR [dimensionless]
!! gfdl_val real, public, parameter :: KAPPA  = RDGAS/CP_AIR             !< RDGAS / CP_AIR [dimensionless]
real(R8),         public, parameter :: KAPPA  = RDGAS/CP_AIR             !< RDGAS / CP_AIR [dimensionless]

!! gfdl_val real, public, parameter :: TFREEZE = 273.16_r8_kind          !< Freezing temperature of fresh water [K]
real(R8),         public, parameter :: TFREEZE = SHR_CONST_TKFRZ         !< Freezing temperature of fresh water [K]

!! gfdl_val real, public, parameter :: STEFAN  = 5.6734e-8_r8_kind !< Stefan-Boltzmann constant [W/m^2/deg^4]
real(R8),         public, parameter :: STEFAN  = SHR_CONST_STEBOL  !< Stefan-Boltzmann constant [W/m^2/deg^4]

!! gfdl_val real, public, parameter :: CP_VAPOR = 4.0_r8_kind*RVGAS      !< Specific heat capacity of water vapor at constant pressure [J/kg/deg]
real(R8),         public, parameter :: CP_VAPOR = SHR_CONST_CPWV         !< Specific heat capacity of water vapor at constant pressure [J/kg/deg]

!! gfdl_val real, public, parameter :: CP_OCEAN = 3989.24495292815_r8_kind !< Specific heat capacity taken from McDougall (2002) 
                                                               !! "Potential Enthalpy ..." [J/kg/deg]
real(R8),         public, parameter :: CP_OCEAN = SHR_CONST_CPSW      !     = 3.996e3_R8      ! specific heat of sea h2o   ~ J/kg/K

!! gfdl_val real, public, parameter :: RHO0    = 1.035e3_r8_kind  !< Average density of sea water [kg/m^3]
real(R8),         public, parameter :: RHO0    = SHR_CONST_RHOSW  !< Average density of sea water [kg/m^3]

!! gfdl_val real, public, parameter :: RHO0R   = 1.0_r8_kind/RHO0 !< Reciprocal of average density of sea water [m^3/kg]
real(R8),         public, parameter :: RHO0R   = 1.0_r8/RHO0      !< Reciprocal of average density of sea water [m^3/kg]

!! gfdl_val real, public, parameter :: RHO_CP  = RHO0*CP_OCEAN    !< (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C) [J/m^3/deg]
real(R8),         public, parameter :: RHO_CP  = RHO0*CP_OCEAN    !< (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C) [J/m^3/deg]

!! gfdl_val real, public, parameter :: ES0 = 1.0_r8_kind        !< Humidity factor. Controls the humidity content of the atmosphere through
                                                                !! the Saturation Vapour Pressure expression when using DO_SIMPLE. [dimensionless]
real(R8),         public, parameter :: ES0 = 1.0_r8             !< Humidity factor. Controls the humidity content of the atmosphere through

!! gfdl_val real, public, parameter :: DENS_H2O = 1000._r8_kind !< Density of liquid water [kg/m^3]
real(R8),         public, parameter :: DENS_H2O = SHR_CONST_RHOFW !< Density of liquid water [kg/m^3]

!! gfdl_val real, public, parameter :: HLS = HLV + HLF          !< Latent heat of sublimation [J/kg]
real(R8),         public, parameter :: HLS = HLV + HLF          !< Latent heat of sublimation [J/kg]

!! gfdl_val real, public, parameter :: WTMAIR   = 2.896440E+01_r8_kind   !< Molecular weight of air [AMU]
real(R8),         public, parameter :: WTMAIR   = SHR_CONST_MWDAIR       !< Molecular weight of air [AMU]

!! gfdl_val real, public, parameter :: WTMH2O   = WTMAIR*(RDGAS/RVGAS)   !< Molecular weight of water [AMU]
real(R8),         public, parameter :: WTMH2O   = SHR_CONST_MWWV         !< Molecular weight of water [AMU]

!! gfdl_val real, public, parameter :: WTMOZONE =  47.99820_r8_kind      !< Molecular weight of ozone [AMU]
real(R8),         public, parameter :: WTMOZONE =  47.99820_r8           !< Molecular weight of ozone [AMU]

!! gfdl_val real, public, parameter :: WTMC     =  12.00000_r8_kind      !< Molecular weight of carbon [AMU]
real(R8),         public, parameter :: WTMC     =  12.00000_r8           !< Molecular weight of carbon [AMU]

!! gfdl_val real, public, parameter :: WTMCO2   =  44.00995_r8_kind      !< Molecular weight of carbon dioxide [AMU]
real(R8),         public, parameter :: WTMCO2   =  44.00995_r8           !< Molecular weight of carbon dioxide [AMU]

!! gfdl_val real, public, parameter :: WTMCH4   =  16.0425_r8_kind       !< Molecular weight of methane [AMU]
real(R8),         public, parameter :: WTMCH4   =  16.0425_r8            !< Molecular weight of methane [AMU]

!! gfdl_val real, public, parameter :: WTMO2    =  31.9988_r8_kind       !< Molecular weight of molecular oxygen [AMU]
real(R8),         public, parameter :: WTMO2    =  31.9988_r8            !< Molecular weight of molecular oxygen [AMU]

!! gfdl_val real, public, parameter :: WTMCFC11 = 137.3681_r8_kind       !< Molecular weight of CFC-11 (CCl3F) [AMU]
real(R8),         public, parameter :: WTMCFC11 = 137.3681_r8            !< Molecular weight of CFC-11 (CCl3F) [AMU]

!! gfdl_val real, public, parameter :: WTMCFC12 = 120.9135_r8_kind       !< Molecular weight of CFC-21 (CCl2F2) [AMU]
real(R8),         public, parameter :: WTMCFC12 = 120.9135_r8            !< Molecular weight of CFC-21 (CCl2F2) [AMU]

!! gfdl_val real, public, parameter :: WTMN     =  14.0067_r8_kind       !< Molecular weight of Nitrogen [AMU]
real(R8),         public, parameter :: WTMN     =  14.0067_r8            !< Molecular weight of Nitrogen [AMU]

!! gfdl_val real, public, parameter :: DIFFAC   = 1.660000E+00_r8_kind   !< Diffusivity factor [dimensionless]
real(R8),         public, parameter :: DIFFAC   = 1.660000E+00_r8        !< Diffusivity factor [dimensionless]

!! gfdl_val real, public, parameter :: AVOGNO   = 6.023000E+23_r8_kind   !< Avogadro's number [atoms/mole]
real(R8),         public, parameter :: AVOGNO   = SHR_CONST_AVOGAD * 1.0e-3_r8 !< convert cesm atoms/kmole to Avogadro's number [atoms/mole]

!! gfdl_val real, public, parameter :: PSTD     = 1.013250E+06_r8_kind   !< Mean sea level pressure [dynes/cm^2]
real(R8),         public, parameter :: PSTD     = SHR_CONST_PSTD*10.0_r8 !< convert cesm units N/m^2 to dynes/cm^2

!! gfdl_val real, public, parameter :: PSTD_MKS = 101325.0_r8_kind       !< Mean sea level pressure [N/m^2]
real(R8),         public, parameter :: PSTD_MKS = SHR_CONST_PSTD         !< Mean sea level pressure [N/m^2]

!! gfdl_val real, public, parameter :: SECONDS_PER_DAY    = 8.640000E+04_r8_kind !< Seconds in a day [s]
real(R8),         public, parameter :: SECONDS_PER_DAY    = SHR_CONST_CDAY       !< Seconds in a day [s]

!! gfdl_val real, public, parameter :: SECONDS_PER_HOUR   = 3600._r8_kind        !< Seconds in an hour [s]
real(R8),         public, parameter :: SECONDS_PER_HOUR   = 3600._r8             !< Seconds in an hour [s]

!! gfdl_val real, public, parameter :: SECONDS_PER_MINUTE = 60._r8_kind          !< Seconds in a minute [s]
real(R8),         public, parameter :: SECONDS_PER_MINUTE = 60._r8               !< Seconds in a minute [s]

!! gfdl_val real, public, parameter :: RAD_TO_DEG         = 180._r8_kind/PI      !< Degrees per radian [deg/rad]
real(R8),         public, parameter :: RAD_TO_DEG         = 180._r8     /PI      !< Degrees per radian [deg/rad]

!! gfdl_val real, public, parameter :: DEG_TO_RAD         = PI/180._r8_kind      !< Radians per degree [rad/deg]
real(R8),         public, parameter :: DEG_TO_RAD         = PI/180._r8           !< Radians per degree [rad/deg]

!! gfdl_val real, public, parameter :: RADIAN             = RAD_TO_DEG           !< Equal to RAD_TO_DEG for backward compatability. [rad/deg]
real(R8),         public, parameter :: RADIAN             = RAD_TO_DEG           !< Equal to RAD_TO_DEG for backward compatability. [rad/deg]

!! gfdl_val real, public, parameter :: ALOGMIN            = -50.0_r8_kind        !< Minimum value allowed as argument to log function [N/A]
real(R8),         public, parameter :: ALOGMIN            = -50.0_r8             !< Minimum value allowed as argument to log function [N/A]

!! gfdl_val real, public, parameter :: EPSLN              = 1.0e-40_r8_kind      !< A small number to prevent divide by zero exceptions [N/A]
real(R8),         public, parameter :: EPSLN              = 1.0e-40_r8           !< A small number to prevent divide by zero exceptions [N/A]

!! gfdl_val real, public, parameter :: RADCON = ((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY !< Factor used to convert flux divergence to
real(R8),         public, parameter :: RADCON = ((1.0E+02_r8*GRAV)/(1.0E+04_r8*CP_AIR))*SECONDS_PER_DAY !< Factor used to convert flux divergence to
                                                                                      !! heating rate in degrees per day [deg sec/(cm day)]
!! gfdl_val real, public, parameter :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY !< Factor used to convert flux divergence to
real(R8),         public, parameter :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY !< Factor used to convert flux divergence to
                                                                       !! heating rate in degrees per day [deg sec/(m day)]
!! gfdl_val real, public, parameter :: O2MIXRAT    = 2.0953E-01_r8_kind !< Mixing ratio of molecular oxygen in air [dimensionless]
real(R8),         public, parameter :: O2MIXRAT    = 2.0953E-01_r8      !< Mixing ratio of molecular oxygen in air [dimensionless]

!! gfdl_val real, public, parameter :: RHOAIR      = 1.292269_r8_kind   !< Reference atmospheric density [kg/m^3]
real(R8),         public, parameter :: RHOAIR      = SHR_CONST_RHODAIR  !< Reference atmospheric density [kg/m^3]

!! gfdl_val real, public, parameter :: VONKARM     = 0.40_r8_kind       !< Von Karman constant [dimensionless]
real(R8),         public, parameter :: VONKARM     = SHR_CONST_KARMAN   !< Von Karman constant [dimensionless]

!! gfdl_val real, public, parameter :: C2DBARS     = 1.e-4_r8_kind      !< Converts rho*g*z (in mks) to dbars: 1dbar = 10^4 (kg/m^3)(m/s^2)m [dbars]
real(R8),         public, parameter :: C2DBARS     = 1.e-4_r8           !< Converts rho*g*z (in mks) to dbars: 1dbar = 10^4 (kg/m^3)(m/s^2)m [dbars]

!! gfdl_val real, public, parameter :: KELVIN      = 273.15_r8_kind     !< Degrees Kelvin at zero Celsius [K]
real(R8),         public, parameter :: KELVIN      = 273.15_r8          !< Degrees Kelvin at zero Celsius [K]

public :: constants_init

contains

!> \brief dummy routine.
subroutine constants_init

end subroutine constants_init

end module constants_mod
