!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  real :: xmin, xmax, ymax
  integer :: lrefine_max, nblockx, nblocky
  character(len=MAX_STRING_LENGTH) :: str

!!  call RuntimeParameters_get('ymin', ymin) !!added

  call RuntimeParameters_get('sim_targLength', sim_targLength)
  call RuntimeParameters_get('sim_targ2Length', sim_targ2Length)
  call RuntimeParameters_get('sim_targ3Length', sim_targ3Length)
  call RuntimeParameters_get('sim_targWidth', sim_targWidth)
  call RuntimeParameters_get('sim_targ2Width', sim_targ2Width)
  call RuntimeParameters_get('sim_targ3Width', sim_targ3Width)
  call RuntimeParameters_get('sim_targThick', sim_targThick)
  call RuntimeParameters_get('sim_targ2Thick', sim_targ2Thick)
  call RuntimeParameters_get('sim_targ3Thick', sim_targ3Thick)
  call RuntimeParameters_get('sim_targRadius', sim_targRadius)
  call RuntimeParameters_get('sim_targDist', sim_targDist)

  call RuntimeParameters_get('sim_targetRadius', sim_targetRadius)
  call RuntimeParameters_get('sim_targetHeight', sim_targetHeight)
  call RuntimeParameters_get('sim_vacuumHeight', sim_vacuumHeight)
  
  call RuntimeParameters_get('sim_rhoTarg', sim_rhoTarg)
  call RuntimeParameters_get('sim_teleTarg', sim_teleTarg)
  call RuntimeParameters_get('sim_tionTarg', sim_tionTarg)
  call RuntimeParameters_get('sim_tradTarg', sim_tradTarg)

  call RuntimeParameters_get('sim_rhoTar2', sim_rhoTar2)
  call RuntimeParameters_get('sim_teleTar2', sim_teleTar2)
  call RuntimeParameters_get('sim_tionTar2', sim_tionTar2)
  call RuntimeParameters_get('sim_tradTar2', sim_tradTar2)

  call RuntimeParameters_get('sim_rhoTar3', sim_rhoTar3)
  call RuntimeParameters_get('sim_teleTar3', sim_teleTar3)
  call RuntimeParameters_get('sim_tionTar3', sim_tionTar3)
  call RuntimeParameters_get('sim_tradTar3', sim_tradTar3)

!!  call RuntimeParameters_get('sim_magx', sim_magx)
!!  call RuntimeParameters_get('sim_magy', sim_magy)

  call RuntimeParameters_get('sim_rhoCham', sim_rhoCham)
  call RuntimeParameters_get('sim_teleCham', sim_teleCham)
  call RuntimeParameters_get('sim_tionCham', sim_tionCham)
  call RuntimeParameters_get('sim_tradCham', sim_tradCham)

  call RuntimeParameters_get('smallX', sim_smallX)

  call RuntimeParameters_get('sim_initGeom', sim_initGeom)

#ifdef FLASH_USM_MHD
  call RuntimeParameters_get('killdivb', sim_killdivb)
#endif
end subroutine Simulation_init
