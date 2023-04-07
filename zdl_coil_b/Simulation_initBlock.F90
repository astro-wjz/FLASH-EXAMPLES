!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(IN) :: blockID)
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData,&
                           Grid_getBlkPtr,         &
                           Grid_releaseBlkPtr

  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction
  ! (including guardcells)

  integer, intent(in) :: blockId

  integer :: i, j, k, n
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: tradActual
  real :: loc_x, loc_y, loc_z, frac
  real :: rho, tele, trad, tion, zbar, abar,magx,magy,magz,magp
  integer :: species
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData
#if NDIM > 0
  real, pointer, dimension(:,:,:,:) :: facezData
#endif

#ifndef CHAM_SPEC
  integer :: CHAM_SPEC = 1, TARG_SPEC = 2, TAR2_SPEC=3
#endif


  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))
  allocate(zcent(blkLimitsGC(HIGH, KAXIS)))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., &
       zcent, blkLimitsGC(HIGH, KAXIS))

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM>2) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  !------------------------------------------------------------------------------

  ! Loop over cells and set the initial state
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           species = CHAM_SPEC
           if (sim_initGeom == "slab") then
              if(NDIM == 2) then
                 if ( xcent(i) <= -sim_targDist+sim_targThick .and. xcent(i) >= -sim_targDist .and. &
                      ycent(j) <= abs(sim_targLength)/2.0 ) then !! CH
                    species = TARG_SPEC
                 endif
                 if ( xcent(i) >= sim_targAmplitude*sin(2*PI*ycent(j)/sim_targPeriod) .and. &
                      xcent(i) <= sim_targAmplitude*sin(2*PI*ycent(j)/sim_targPeriod) &
                       + sim_targ2Thick .and. &
                      ycent(j) <= abs(sim_targ2Length)/2.0 ) then  !! shelf
                    species = TAR2_SPEC
                 end if 
!                  if ( xcent(i) >= sim_targDist .and. xcent(i)<=sim_targDist+sim_targ2Thick .and. &
!                       abs(ycent(j)) <= sim_targ3Length ) then
!                    species = TAR2_SPEC
!                  end if         
              endif
!               if(NDIM == 3) then
!                  if ( xcent(i) <= sim_targThick .and. xcent(i) >= 0.0 .and. &
!                       ycent(j) <= 0.0 .and. ycent(j) >= -sim_targLength .and. &
!                       abs(zcent(k)) <= sim_targWidth/2.0 ) then !! CH
!                     species = TARG_SPEC
!                  endif
!                  if ( xcent(i) <= sim_targ2Thick .and. xcent(i) >= 0.0 .and. &
!                       ycent(j) <= sim_targ2Length .and. ycent(j) >= 0.0 .and. &
!                       abs(zcent(k)) <= sim_targ2Width/2.0 ) then  !! shelf
!                     species = TAR2_SPEC
!                  end if
!                  if ( xcent(i) <= sim_targDist .and. ycent(j) <= 0.0 .and. &
!                       ycent(j) >= -sim_targ3Length .and. &
!                       (xcent(i)-sim_targDist)**2 + zcent(k)**2 <= sim_targRadius**2 ) then  !! Cu
!                     species = TAR3_SPEC
!                  end if
!               end if
           else
               if (sqrt(xcent(i)**2+ycent(j)**2+zcent(k)**2)<= sim_targetRadius) then
                  species = TARG_SPEC
              end if
           end if

           if(species == TARG_SPEC) then
              rho = sim_rhoTarg
              tele = sim_teleTarg
              tion = sim_tionTarg
              trad = sim_tradTarg
           else if(species==TAR2_SPEC) then
              rho = sim_rhoTar2
              tele = sim_teleTar2
              tion = sim_tionTar2
              trad = sim_tradTar2
!            else if(species==TAR3_SPEC) then
!               rho = sim_rhoTar3
!               tele = sim_teleTar3
!               tion = sim_tionTar3
!               trad = sim_tradTar3
           else
              rho = sim_rhoCham
              tele = sim_teleCham
              tion = sim_tionCham
              trad = sim_tradCham
           end if
           
!           if (NDIM == 2 .and. useDipole == .true. ) then
!               loc_x = xcent(i)
!               loc_y = ycent(j)
!               frac = (loc_x**2 + loc_y**2)**2.5
!               if (xcent(i)<=-sim_targAmplitude) then
!                  sim_magx = -1.77147/loc_x**2*3*loc_x*loc_y/frac
!                  sim_magy = -1.77147/loc_x**2*(2*loc_y**2-loc_x**2)/frac
!               else if (xcent(i)<=sim_targAmplitude ) then
!                  if (frac==0) then
!                     sim_magx = 0.0
!                     sim_magy = 150000.0
!                  else
!                     if (loc_x<=0) then
!                        sim_magx = -205761317*loc_x**6*3*loc_x*loc_y/frac
!                        sim_magy = -205761317*loc_x**6*(2*loc_y**2-loc_x**2)/frac + 150000.0
!                     else
!                        sim_magx = 205761317*loc_x**6*3*loc_x*loc_y/frac
!                        sim_magy = 205761317*loc_x**6*(2*loc_y**2-loc_x**2)/frac + 150000.0
!                     end if 
!                  end if
!               else
!                  sim_magx = 0.0
!                  sim_magy = 0.0
!               end if
!           end if

          if (NDIM == 2 .and. useDipole == .true. ) then
              loc_x = xcent(i)
              loc_y = ycent(j)
              frac = (loc_x**2 + loc_y**2)**2.5
              if (xcent(i)<=-sim_targAmplitude) then
                 sim_magx = -1.77147/loc_x**2*3*loc_x*loc_y/frac
                 sim_magy = -1.77147/loc_x**2*(2*loc_y**2-loc_x**2)/frac
!               else if (xcent(i)<=sim_targAmplitude ) then
!                  if (frac==0) then
!                     sim_magx = 0.0
!                     sim_magy = 150000.0
!                  else
!                     if (loc_x<=0) then
!                        sim_magx = -205761317*loc_x**6*3*loc_x*loc_y/frac
!                        sim_magy = -205761317*loc_x**6*(2*loc_y**2-loc_x**2)/frac + 150000.0
!                     else
!                        sim_magx = 205761317*loc_x**6*3*loc_x*loc_y/frac
!                        sim_magy = 205761317*loc_x**6*(2*loc_y**2-loc_x**2)/frac + 150000.0
!                     end if 
!                  end if
              else
                 loc_x = -sim_targAmplitude
                 frac = (loc_x**2 + loc_y**2)**2.5
                 sim_magx = -1.77147/loc_x**2*3*loc_x*loc_y/frac
                 sim_magy = -1.77147/loc_x**2*(2*loc_y**2-loc_x**2)/frac
              end if
          end if

!            if (NDIM == 2 .and. useDipole == .true. ) then
!                !loc_x = xcent(i) - sim_targDist
!                !loc_y = ycent(j)
!                !frac = (loc_x**2 + loc_y**2)**2.5
!                if (xcent(i)<=sim_targDist-sim_targAmplitude ) then
!                   sim_magy = exp(sim_magfac*xcent(i))
!                   sim_magx = 0.0
!                else if (xcent(i)>sim_targDist-sim_targAmplitude ) then
!                   sim_magx = 0.0
!                   sim_magy = exp(sim_magfac*(sim_targDist-sim_targAmplitude))
!                else
!                   sim_magx = 0.0
!                   sim_magy = 0.0
!                end if
!            end if

           if (NDIM>=2) then
               magx = sim_magx
               magy = sim_magy
               magp = 0.5*(magx*magx + magy*magy)
              if (NDIM==3) then
                  magz = sim_magz
                  magp = 0.5*(magx*magx + magy*magy + magz*magz)
              end if
           end if

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tele)
           
           call Grid_putPointData(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, magx)
           call Grid_putPointData(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, magy)
           if (NDIM>2) then
              call Grid_putPointData(blockId, CENTER, MAGZ_VAR, EXTERIOR, axis, magz)
           end if
           call Grid_putPointData(blockId, CENTER, MAGP_VAR, EXTERIOR, axis, magp)
           call Grid_putPointData(blockId, CENTER, DIVB_VAR, EXTERIOR, axis, 0.0)

#ifdef FLASH_3T
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)

           ! Set up radiation energy density:
           call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)
#endif
           if (NSPECIES > 0) then
              ! Fill mass fractions in solution array if we have any SPECIES defined.
              ! We put nearly all the mass into either the Xe material if XE_SPEC is defined,
              ! or else into the first species.
              do n = SPECIES_BEGIN,SPECIES_END
                 if (n==species) then
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, 1.0e0-(NSPECIES-1)*sim_smallX)
                 else
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, sim_smallX)
                 end if
              enddo
           end if

#ifdef BDRY_VAR
           call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, -1.0)
#endif


#if NFACE_VARS > 0
           !! In this case we initialized Az using the cell-cornered coordinates.
           if (sim_killdivb) then
              if (NDIM >= 2) then
                 facexData(MAG_FACE_VAR,i,j,k)= magx
                 faceyData(MAG_FACE_VAR,i,j,k)= magy
                 if (NDIM>2) facezData(MAG_FACE_VAR,i,j,k)= magz
              endif
           endif
#endif
        enddo
     enddo
  enddo
#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM>2) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif
  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

  return

end subroutine Simulation_initBlock
