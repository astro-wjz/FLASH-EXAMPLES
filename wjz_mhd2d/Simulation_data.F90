!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!!
!!***
module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!

  real, save :: sim_targLength
  real, save :: sim_targ2Length
  real, save :: sim_targ3Length
  real, save :: sim_targWidth
  real, save :: sim_targ2Width
  real, save :: sim_targ3Width
  real, save :: sim_targThick
  real, save :: sim_targ2Thick
  real, save :: sim_targ3Thick
  real, save :: sim_targRadius
  real, save :: sim_targDist

  real, save :: sim_targetRadius
  real, save :: sim_targetHeight
  real, save :: sim_vacuumHeight

  real,    save :: sim_rhoTarg  
  real,    save :: sim_teleTarg 
  real,    save :: sim_tionTarg 
  real,    save :: sim_tradTarg 
  real,    save :: sim_zminTarg
  integer, save :: sim_eosTarg

  real,    save :: sim_rhoTar2
  real,    save :: sim_teleTar2
  real,    save :: sim_tionTar2
  real,    save :: sim_tradTar2
  real,    save :: sim_zminTar2
  integer, save :: sim_eosTar2

  real,    save :: sim_rhoTar3
  real,    save :: sim_teleTar3
  real,    save :: sim_tionTar3
  real,    save :: sim_tradTar3
  real,    save :: sim_zminTar3
  integer, save :: sim_eosTar3

  real,    save :: sim_rhoCham  
  real,    save :: sim_teleCham 
  real,    save :: sim_tionCham 
  real,    save :: sim_tradCham 
  integer, save :: sim_eosCham  

  logical, save :: sim_killdivb = .FALSE.
  real, save :: sim_smallX
  character(len=MAX_STRING_LENGTH), save :: sim_initGeom


end module Simulation_data
