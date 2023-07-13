!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!! \brief Defines useful constants for Earth. Constants are defined as real
!! parameters. Constants are accessed through the "use" statement.
!!
!! \author Bruce Wyman <Bruce.Wyman@noaa.gov>
module constants_mod

#ifdef USE_FMSCONST

!---variable for strong typing grid parameters
use platform_mod, only: r8_kind
implicit none
private

! Include variable "version" to be written to log file.
#include<file_version.h>
!-----------------------------------------------------------------------
! version is public so that write_version_number can be called for constants_mod
! by fms_init
public :: version

real :: realnumber !< dummy variable to use in HUGE initializations

!! The small_fac parameter is used to alter the radius of the earth to allow one to
!! examine non-hydrostatic effects without the need to run full-earth high-resolution
!! simulations (<13km) that will tax hardware resources.
#ifdef SMALL_EARTH
#if defined(DCMIP) || (defined(HIWPP) && defined(SUPER_K))
 real, public, parameter :: small_fac =  1._r8_kind / 120._r8_kind   ! only needed for supercell test
#elif defined(HIWPP)
 real, public, parameter :: small_fac = 1._r8_kind / 166.7_r8_kind
#else
 real, public, parameter :: small_fac = 1._r8_kind / 10._r8_kind
#endif
#else
 real, public, parameter :: small_fac = 1._r8_kind
#endif

#ifdef GFS_PHYS
! SJL: the following are from fv3_gfsphysics/gfs_physics/physics/physcons.f90
real,               public, parameter :: RADIUS = 6.3712e+6_r8_kind * small_fac !< Radius of the Earth [m]
real(kind=r8_kind), public, parameter :: PI_8   = 3.1415926535897931_r8_kind    !< Ratio of circle circumference to diameter [N/A]
real,               public, parameter :: PI     = 3.1415926535897931_r8_kind    !< Ratio of circle circumference to diameter [N/A]
real,               public, parameter :: OMEGA  = 7.2921e-5_r8_kind / small_fac !< Rotation rate of the Earth [1/s]
real,               public, parameter :: GRAV   = 9.80665_r8_kind     !< Acceleration due to gravity [m/s^2]
real(kind=r8_kind), public, parameter :: GRAV_8 = 9.80665_r8_kind     !< Acceleration due to gravity [m/s^2] (REAL(KIND=8))
real,               public, parameter :: RDGAS  = 287.05_r8_kind      !< Gas constant for dry air [J/kg/deg]
real,               public, parameter :: RVGAS  = 461.50_r8_kind      !< Gas constant for water vapor [J/kg/deg]
! Extra:
real,               public, parameter :: HLV      = 2.5e6_r8_kind     !< Latent heat of evaporation [J/kg]
real,               public, parameter :: HLF      = 3.3358e5_r8_kind  !< Latent heat of fusion [J/kg]
real,               public, parameter :: con_cliq = 4.1855e+3_r8_kind !< spec heat H2O liq [J/kg/K]
real,               public, parameter :: con_csol = 2.1060e+3_r8_kind !< spec heat H2O ice [J/kg/K]
real,               public, parameter :: CP_AIR = 1004.6_r8_kind      !< Specific heat capacity of dry air at constant pressure [J/kg/deg]
real,               public, parameter :: KAPPA  = RDGAS/CP_AIR        !< RDGAS / CP_AIR [dimensionless]
real,               public, parameter :: TFREEZE = 273.15_r8_kind     !< Freezing temperature of fresh water [K]

#else

real,         public, parameter :: RADIUS = 6371.0e+3_r8_kind * small_fac   !< Radius of the Earth [m]
real(kind=8), public, parameter :: PI_8   = 3.14159265358979323846_r8_kind  !< Ratio of circle circumference to diameter [N/A]
real,         public, parameter :: PI     = 3.14159265358979323846_r8_kind  !< Ratio of circle circumference to diameter [N/A]
real,         public, parameter :: OMEGA  = 7.292e-5_r8_kind / small_fac    !< Rotation rate of the Earth [1/s]
real,         public, parameter :: GRAV   = 9.80_r8_kind             !< Acceleration due to gravity [m/s^2]
real,         public, parameter :: RDGAS  = 287.04_r8_kind           !< Gas constant for dry air [J/kg/deg]
real,         public, parameter :: RVGAS  = 461.50_r8_kind           !< Gas constant for water vapor [J/kg/deg]
! Extra:
real,         public, parameter :: HLV = 2.500e6_r8_kind             !< Latent heat of evaporation [J/kg]
real,         public, parameter :: HLF = 3.34e5_r8_kind              !< Latent heat of fusion [J/kg]
real,         public, parameter :: KAPPA  = 2.0_r8_kind/7.0_r8_kind  !< RDGAS / CP_AIR [dimensionless]
real,         public, parameter :: CP_AIR = RDGAS/KAPPA              !< Specific heat capacity of dry air at constant pressure [J/kg/deg]
real,         public, parameter :: TFREEZE = 273.16_r8_kind          !< Freezing temperature of fresh water [K]
#endif

real, public, parameter :: STEFAN  = 5.6734e-8_r8_kind !< Stefan-Boltzmann constant [W/m^2/deg^4]

real, public, parameter :: CP_VAPOR = 4.0_r8_kind*RVGAS      !< Specific heat capacity of water vapor at constant pressure [J/kg/deg]
real, public, parameter :: CP_OCEAN = 3989.24495292815_r8_kind !< Specific heat capacity taken from McDougall (2002)
                                                               !! "Potential Enthalpy ..." [J/kg/deg]
real, public, parameter :: RHO0    = 1.035e3_r8_kind  !< Average density of sea water [kg/m^3]
real, public, parameter :: RHO0R   = 1.0_r8_kind/RHO0 !< Reciprocal of average density of sea water [m^3/kg]
real, public, parameter :: RHO_CP  = RHO0*CP_OCEAN    !< (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C) [J/m^3/deg]

real, public, parameter :: ES0 = 1.0_r8_kind        !< Humidity factor. Controls the humidity content of the atmosphere through
                                                    !! the Saturation Vapour Pressure expression when using DO_SIMPLE. [dimensionless]
real, public, parameter :: DENS_H2O = 1000._r8_kind !< Density of liquid water [kg/m^3]
real, public, parameter :: HLS = HLV + HLF          !< Latent heat of sublimation [J/kg]

real, public, parameter :: WTMAIR   = 2.896440E+01_r8_kind   !< Molecular weight of air [AMU]
real, public, parameter :: WTMH2O   = WTMAIR*(RDGAS/RVGAS)   !< Molecular weight of water [AMU]
real, public, parameter :: WTMOZONE =  47.99820_r8_kind      !< Molecular weight of ozone [AMU]
real, public, parameter :: WTMC     =  12.00000_r8_kind      !< Molecular weight of carbon [AMU]
real, public, parameter :: WTMCO2   =  44.00995_r8_kind      !< Molecular weight of carbon dioxide [AMU]
real, public, parameter :: WTMCH4   =  16.0425_r8_kind       !< Molecular weight of methane [AMU]
real, public, parameter :: WTMO2    =  31.9988_r8_kind       !< Molecular weight of molecular oxygen [AMU]
real, public, parameter :: WTMCFC11 = 137.3681_r8_kind       !< Molecular weight of CFC-11 (CCl3F) [AMU]
real, public, parameter :: WTMCFC12 = 120.9135_r8_kind       !< Molecular weight of CFC-21 (CCl2F2) [AMU]
real, public, parameter :: WTMN     =  14.0067_r8_kind       !< Molecular weight of Nitrogen [AMU]
real, public, parameter :: DIFFAC   = 1.660000E+00_r8_kind   !< Diffusivity factor [dimensionless]
real, public, parameter :: AVOGNO   = 6.023000E+23_r8_kind   !< Avogadro's number [atoms/mole]
real, public, parameter :: PSTD     = 1.013250E+06_r8_kind   !< Mean sea level pressure [dynes/cm^2]
real, public, parameter :: PSTD_MKS = 101325.0_r8_kind       !< Mean sea level pressure [N/m^2]

real, public, parameter :: SECONDS_PER_DAY    = 8.640000E+04_r8_kind !< Seconds in a day [s]
real, public, parameter :: SECONDS_PER_HOUR   = 3600._r8_kind        !< Seconds in an hour [s]
real, public, parameter :: SECONDS_PER_MINUTE = 60._r8_kind          !< Seconds in a minute [s]
real, public, parameter :: RAD_TO_DEG         = 180._r8_kind/PI      !< Degrees per radian [deg/rad]
real, public, parameter :: DEG_TO_RAD         = PI/180._r8_kind      !< Radians per degree [rad/deg]
real, public, parameter :: RADIAN             = RAD_TO_DEG           !< Equal to RAD_TO_DEG for backward compatability. [rad/deg]
real, public, parameter :: ALOGMIN            = -50.0_r8_kind        !< Minimum value allowed as argument to log function [N/A]
real, public, parameter :: EPSLN              = 1.0e-40_r8_kind      !< A small number to prevent divide by zero exceptions [N/A]

real, public, parameter :: RADCON = ((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY !< Factor used to convert flux divergence to
                                                                                      !! heating rate in degrees per day [deg sec/(cm day)]
real, public, parameter :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY !< Factor used to convert flux divergence to
                                                                       !! heating rate in degrees per day [deg sec/(m day)]
real, public, parameter :: O2MIXRAT    = 2.0953E-01_r8_kind !< Mixing ratio of molecular oxygen in air [dimensionless]
real, public, parameter :: RHOAIR      = 1.292269_r8_kind   !< Reference atmospheric density [kg/m^3]
real, public, parameter :: VONKARM     = 0.40_r8_kind       !< Von Karman constant [dimensionless]
real, public, parameter :: C2DBARS     = 1.e-4_r8_kind      !< Converts rho*g*z (in mks) to dbars: 1dbar = 10^4 (kg/m^3)(m/s^2)m [dbars]
real, public, parameter :: KELVIN      = 273.15_r8_kind     !< Degrees Kelvin at zero Celsius [K]

#else

use shr_kind_mod,   only : R8 => shr_kind_r8
use shr_const_mod,  only : &
     SHR_CONST_PSTD,                      & ! standard pressure ~ pascals
     PI              => SHR_CONST_PI,     &
     PI_8            => SHR_CONST_PI,     &
     SECONDS_PER_DAY => SHR_CONST_CDAY,   &
     OMEGA           => SHR_CONST_OMEGA,  &
     RADIUS          => SHR_CONST_REARTH, &
     GRAV            => SHR_CONST_G,      &
     STEFAN          => SHR_CONST_STEBOL, & ! Stefan-Boltzmann constant ~ W/m^2/K^4
     AVOGNO          => SHR_CONST_AVOGAD, & ! Avogadro's number ~ molecules/kmole
     WTMAIR          => SHR_CONST_MWDAIR, & ! molecular weight dry air ~ kg/kmole
     WTMH2O          => SHR_CONST_MWWV,   & ! molecular weight water vapor
     RDGAS           => SHR_CONST_RDAIR,  & ! Dry air gas constant     ~ J/K/kg
     RVGAS           => SHR_CONST_RWV,    & ! Water vapor gas constant ~ J/K/kg
     VONKARM         => SHR_CONST_KARMAN, & ! Von Karman constant
     TFREEZE         => SHR_CONST_TKFRZ,  & ! freezing T of fresh water          ~ K
     RHOAIR          => SHR_CONST_RHODAIR,& ! density of dry air at STP  ~ kg/m^3
     DENS_H2O        => SHR_CONST_RHOFW,  & ! density of fresh water     ~ kg/m^3
     RHO0            => SHR_CONST_RHOSW,  & ! density of sea water       ~ kg/m^3
     CP_AIR          => SHR_CONST_CPDAIR, & ! specific heat of dry air   ~ J/kg/K
     CP_VAPOR        => SHR_CONST_CPWV,   & ! specific heat of water vap ~ J/kg/K
     CP_OCEAN        => SHR_CONST_CPSW,   &
     HLF             => SHR_CONST_LATICE, & ! latent heat of fusion      ~ J/kg
     HLV             => SHR_CONST_LATVAP    ! latent heat of evaporation ~ J/kg

implicit none
private

public :: PI, PI_8, SECONDS_PER_DAY, OMEGA, RADIUS, GRAV, STEFAN, AVOGNO, &
          WTMAIR, WTMH2O, RDGAS, RVGAS, VONKARM, TFREEZE, RHOAIR, DENS_H2O, &
          RHO0, CP_AIR, CP_VAPOR, CP_OCEAN, HLF, HLV

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

! Constants below use CESM shr values

real(R8),   public, parameter :: KAPPA  = RDGAS/CP_AIR             !< RDGAS / CP_AIR [dimensionless]
real(R8),   public, parameter :: RHO0R   = 1.0_r8/RHO0        	   !< Reciprocal of average density of sea water [m^3/kg]
real(R8),   public, parameter :: RHO_CP  = RHO0*CP_OCEAN           !< (kg/m^3)*(cal/kg/deg C)(joules/cal) =
                                                                   !<(joules/m^3/deg C) [J/m^3/deg]
real(R8),   public, parameter :: ES0 = 1.0_r8                      !< Humidity factor. Controls the humidity content of the
                                                                   !< atmosphere through the Saturation Vapour Pressure
                                                                   !< expression when using DO_SIMPLE. [dimensionless]
real(R8),   public, parameter :: HLS = HLV + HLF                   !< Latent heat of sublimation [J/kg]
real(R8),   public, parameter :: WTMOZONE =  47.99820_r8           !< Molecular weight of ozone [AMU]
real(R8),   public, parameter :: WTMC     =  12.00000_r8           !< Molecular weight of carbon [AMU]
real(R8),   public, parameter :: WTMCO2   =  44.00995_r8           !< Molecular weight of carbon dioxide [AMU]
real(R8),   public, parameter :: WTMCH4   =  16.0425_r8            !< Molecular weight of methane [AMU]
real(R8),   public, parameter :: WTMO2    =  31.9988_r8            !< Molecular weight of molecular oxygen [AMU]
real(R8),   public, parameter :: WTMCFC11 = 137.3681_r8            !< Molecular weight of CFC-11 (CCl3F) [AMU]
real(R8),   public, parameter :: WTMCFC12 = 120.9135_r8            !< Molecular weight of CFC-21 (CCl2F2) [AMU]
real(R8),   public, parameter :: WTMN     =  14.0067_r8            !< Molecular weight of Nitrogen [AMU]
real(R8),   public, parameter :: DIFFAC   = 1.660000E+00_r8        !< Diffusivity factor [dimensionless]
real(R8),   public, parameter :: PSTD_MKS = SHR_CONST_PSTD         !< standard pressure ~ pascals
real(R8),   public, parameter :: PSTD     = SHR_CONST_PSTD*10.0_r8 !< convert cesm units N/m^2 to dynes/cm^2
real(R8),   public, parameter :: SECONDS_PER_HOUR   = 3600._r8     !< Seconds in an hour [s]
real(R8),   public, parameter :: SECONDS_PER_MINUTE = 60._r8       !< Seconds in a minute [s]
real(R8),   public, parameter :: RAD_TO_DEG         = 180._r8/PI   !< Degrees per radian [deg/rad]
real(R8),   public, parameter :: DEG_TO_RAD         = PI/180._r8   !< Radians per degree [rad/deg]
real(R8),   public, parameter :: RADIAN             = RAD_TO_DEG   !< Equal to RAD_TO_DEG for backward compatability. [rad/deg]
real(R8),   public, parameter :: ALOGMIN            = -50.0_r8     !< Minimum value allowed as argument to log function [N/A]
real(R8),   public, parameter :: EPSLN              = 1.0e-40_r8   !< A small number to prevent divide by zero exceptions [N/A]
real(R8),   public, parameter :: RADCON = ((1.0E+02_r8*GRAV)/(1.0E+04_r8*CP_AIR))*SECONDS_PER_DAY !< convert flux divergence
                                                                   !to heating rate in degrees per day [deg sec/(cm day)]
real(R8),   public, parameter :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY !< Factor used to convert flux divergence to
                                                                   !< heating rate in degrees per day [deg sec/(m day)]
real(R8),   public, parameter :: O2MIXRAT    = 2.0953E-01_r8       !< Mixing ratio of molecular oxygen in air [dimensionless]
real(R8),   public, parameter :: C2DBARS     = 1.e-4_r8            !< rho*g*z(mks) to dbars: 1dbar = 10^4 (kg/m^3)(m/s^2)m [dbars]
real(R8),   public, parameter :: KELVIN      = 273.15_r8           !< Degrees Kelvin at zero Celsius [K]

#ifdef SMALL_EARTH
#if defined(DCMIP) || (defined(HIWPP) && defined(SUPER_K))
 real, public, parameter :: small_fac =  1._r8 / 120._r8   ! only needed for supercell test
#elif defined(HIWPP)
 real, public, parameter :: small_fac = 1._r8 / 166.7_r8
#else
 real, public, parameter :: small_fac = 1._r8 / 10._r8
#endif
#else
 real, public, parameter :: small_fac = 1._r8
#endif

#endif    ! ifdef USE_FMSCONST

public :: constants_init

contains

!> \brief dummy routine.
subroutine constants_init

end subroutine constants_init

end module constants_mod

!   <FUTURE>
!   1.  Renaming of constants.
!   </FUTURE>
!   <FUTURE>
!   2.  Additional constants.
!   </FUTURE>
!   <NOTE>
!    Constants have been declared as type REAL, PARAMETER.
!
!    The value a constant can not be changed in a users program.
!    New constants can be defined in terms of values from the
!    constants module using a parameter statement.<br><br>
!
!    The name given to a particular constant may be changed.<br><br>
!
!    Constants can be used on the right side on an assignment statement
!    (their value can not be reassigned).
!
!
!<TESTPROGRAM NAME="EXAMPLE">
!<PRE>
!    use constants_mod, only:  TFREEZE, grav_new => GRAV
!    real, parameter :: grav_inv = 1.0 / grav_new
!    tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!    geopotential(:,:) = height(:,:) * grav_new
!</PRE>
!</TESTPROGRAM>
!   </NOTE>
