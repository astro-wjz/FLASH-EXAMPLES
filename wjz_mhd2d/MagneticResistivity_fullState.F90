!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/MagneticResistivity_fullState
!!
!! NAME
!!  MagneticResistivity_fullState
!!
!! SYNOPSIS
!!  call MagneticResistivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                                     real(out)   :: resPar,
!!                            OPTIONAL,real(out)   :: resPerp)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron Magnetic Resistivity for all materials,
!! including those with Z > 1. The expressions used here comes from NRL plasma form.
!!
!!  Returns Magnetic Resistivity, parallel and perpendicular components.
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   resPar   :   parallel component of Magnetic Resistivity
!!   resPerp :    perpendicular component of Magnetic Resistivity
!!
!!***

#include "Flash.h"
#include "constants.h"  


subroutine MagneticResistivity_fullState(solnVec,resPar, resPerp)
  use MagneticResistivity_interface, ONLY: MagneticResistivity
  use MagneticResistivity_data, ONLY: mag_useMagneticResistivity, &
       res_mele, res_qele, res_navo, res_speedlt, res_boltz, res_hbar
  use MagneticResistivity_data, ONLY: res_mUnit
  use MagneticResistivity_data, ONLY: res_coef
  use MagneticResistivity_data, ONLY: res_maxRes
  use Eos_interface, ONLY: Eos_getAbarZbar

  implicit none

  real, intent(in)  :: solnVec(NUNK_VARS)
  real, intent(out) :: resPar
  real, OPTIONAL, intent(out) :: resPerp

  real :: dens
  real :: tele, tele_eV
  real :: tion
  real :: nele
  real :: nion
  real :: abar
  real :: zbar
  real :: eqtime
  real :: resPerpLoc
  real :: ll

  call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)

  dens = solnVec(DENS_VAR)
  nion = dens * res_navo / abar
  nele = zbar * nion

#ifdef FLASH_3T
  tele = solnVec(TELE_VAR)
  tion = solnVec(TION_VAR)
#else
  tele = solnVec(TEMP_VAR)
  tion = solnVec(TEMP_VAR)
#endif

  tele_eV = tele/11604.5221
  tele_eV = max(tele_eV,0.1)
  call loglambda(tele, zbar*nion, zbar, ll)

  !resPerpLoc = 8.196e5*max(zbar,1.0)*ll/tele_eV**1.5  !! In CGS -- Here tele has to be in eV added max(zbar,1)
  resPerpLoc = tele_eV/(tele_eV**2.5/(8.196e5*max(zbar,1.0)*ll)+5e6/8e9*(dens/8.96) + 3e5*tele_eV/8e9*(dens/8.96))  !! In CGS -- Here tele has to be in eV added max(zbar,1),Shu: For high density
  !!Cu, I used the formula from J R Davis PRE 65 026407

  ! This formula is only valid when the magnetic field is strong (and
  ! only for hydrogen):
  resPar = resPerpLoc/1.96

  !! FOR CGS--> SI unit conversion
  if (res_mUnit == "SI" .or. res_mUnit == "si" ) then
     resPerpLoc = resPerpLoc*1.0e-4
     resPar  = resPar *1.0e-4
  end if


  if (present(resPerp)) then
    resPerp = resPerpLoc
    resPerp = min(resPerp, res_maxRes) ! cap resPerp to avoid unphysically cold regions

   !  write(*,*) resPerp,'Perp full state'

  end if

  !! cap resPar to avoid unphysically cold regions
  resPar = min(resPar, res_maxRes)
  
  !write(*,*) resPar,'Par full state'
contains

  subroutine loglambda(tele, nele, zbar, ll)
    implicit none

    ! This subroutine computes the Coulomb logarithm. The formula used
    ! comes from Atzeni.
    
    real, intent(in)  :: tele ! electron temperature [K]
    real, intent(in)  :: nele ! electron number density [cm^-3]
    real, intent(in)  :: zbar ! the average ionization [unitless]
    real, intent(out) :: ll   ! the coulomb logarithm [unitless]

    real :: bmax, bmin, bmin_classic, bmin_quantum
    real, parameter :: ll_floor = 1.0
    
    bmax = sqrt(res_boltz * tele / (4*PI * res_qele**2 * nele))
    
    bmin_classic = zbar * res_qele**2 / (3*res_boltz*tele)
    bmin_quantum = res_hbar / (2*sqrt(3*res_boltz*tele*res_mele))
    bmin = max(bmin_classic, bmin_quantum)

    ll = max(log(bmax/bmin), ll_floor)

  end subroutine loglambda

end subroutine MagneticResistivity_fullState
