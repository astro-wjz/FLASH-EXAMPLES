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
  real :: rho, tele, trad, tion, zbar, abar,bx,by
  integer :: species
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData
#if NDIM > 0
  real, pointer, dimension(:,:,:,:) :: facezData
#endif

#ifndef CHAM_SPEC
  integer :: CHAM_SPEC = 1, TARG_SPEC = 2, TAR2_SPEC=3, TAR3_SPEC=4,  WALL_SPEC=5
#endif
  real, dimension(128,256) :: Bxinit,Byinit
  integer :: row, col
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

  

   open(unit=1,file="BX.txt")
   do row=1,128
        read(1,*) (Bxinit(row,col),col=1,256)
   end do
   close(unit=1)

   open(unit=1,file="BZ.txt")
   do row=1,128
        read(1,*) (Byinit(row,col),col=1,256)
   end do
   close(unit=1)



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
              if(NDIM == 1) then
                 if ( xcent(i) <= sim_targetHeight + sim_vacuumHeight .and. &
                      xcent(i) >= sim_vacuumHeight ) then
                    species = TARG_SPEC
                 end if

              elseif(NDIM == 2) then
#ifdef WALL_SPEC

                if ( ABS(xcent(i)) <= sim_walloutRadius .and. &  !!added for wall
                     ABS(xcent(i)) >= sim_wallinRadius .and. &
                    ycent(j) <= sim_wallYmax .and. &  !!here I just use the Ymin and Ymax to discribe the height of target.
                    ycent(j) >= sim_wallYmin ) then
                    species = WALL_SPEC
                end if
#endif

               ! if ( (ABS(xcent(i)) <= sim_targetRadius) .and. &
               !      ycent(j) <= sim_targetYmax .and. &
               !      ycent(j) >= sim_targetYmin ) then
               !    species = TARG_SPEC
               ! end if
                if ((xcent(i)**2+(ycent(j)+0.03)**2)<0.03**2) then
                        species = TARG_SPEC
                end if
                if ( ABS(xcent(i)) <= 50e-4 .and. & !! coil ring ! ABS(xcent(i)) >= 450e-4.and. &
                     ycent(j) >= 1580e-4 .and. &
                     ycent(j) <= 1630e-4) then
                   species = TARG_SPEC
                end if
                !blocker ! if (ycent(j) >= 950e-4 .and. &
               !     ycent(j) <= 1000e-4 .and. &
               !     ABS(xcent(i))<=150e-4)then
               !     species=TARG_SPEC
               ! end if
                bx=bxinit(FLOOR((ycent(j)+275e-4)/10e-4)+1,FLOOR((xcent(i)-5e-4)/10e-4)+1)*1e4/2/1.7724
                by=byinit(FLOOR((ycent(j)+275e-4)/10e-4)+1,FLOOR((xcent(i)-5e-4)/10e-4)+1)*1e4/2/1.7724


#ifdef TAR2_SPEC
                if (  ABS(xcent(i)) <= sim_target2Radius .and. &
                          ycent(j)>=sim_target2Ymin .and. &
                            ycent(j)<sim_wallYmin) then
                        species = TAR2_SPEC

                        !write(*,*) "cone bottom"
                end if
#endif
#ifdef TAR3_SPEC
                if (  ABS(xcent(i)) <= sim_target3Radius .and. &
                              (ycent(j)<sim_target3Ymax .and. &
                              ycent(j)>sim_wallYmax)) then
                                species = TAR3_SPEC
                end if
#endif

                    !    write(*,*) "cone top"
#ifdef CONE_SPEC
                if( ABS(ycent(j)) > sim_wallYmax .and. &
                             ABS(ycent(j))>(sim_target2a1*ABS(xcent(i))+sim_target2b1) .and. &
                             ABS(ycent(j))<(sim_target2a1*ABS(xcent(i))+sim_target2b2))then
                      !       write(*,*) "cone"

                          !! the sentence below is for cylinder cover
                          !!ABS(xcent(i)) > sim_target2inRadius) .and. &
                          !!   (  (ycent(j)<=sim_target2outYmax .and. &
                          !!     ycent(j)>sim_wallYmax ) .or. &
                          !!    (ycent(j)>=sim_target2outYmin .and. &
                          !!     ycent(j)<sim_wallYmin) )) then
                          !! over
                          species = CONE_SPEC
                end if
#endif
              elseif(NDIM==3) then
                if ( (ABS(xcent(i)) <= sim_targetRadius) .and. &
                     (ABS(ycent(j)) <= sim_targetRadius) .and. &
                     zcent(k) <= sim_targetYmax .and. &
                     zcent(k) >= sim_targetYmin ) then
                   species = TARG_SPEC
                end if

                if ( (ABS(xcent(i)) <= sim_targetRadius) .and. & !!up plate
                    (ABS(ycent(j)) <= sim_targetRadius) .and. &
                    xcent(i)**2+ycent(j)**2 >= 0.02**2 .and. & !! hole
                    zcent(k) <= sim_targetYmax+sim_vacuumHeight .and. &
                    zcent(k) >= sim_targetYmin+sim_vacuumHeight ) then
                  species = TARG_SPEC
                end if

                if ( xcent(i)>= sim_targetRadius .and. &
                     xcent(i)<= sim_targetRadius+0.05 .and. &
                       (ABS(ycent(j)-0.03)<=(sim_targetYmax-sim_targetYmin)/2 .or. &
                       ABS(ycent(j)+0.03)<=(sim_targetYmax-sim_targetYmin)/2) .and. &
                       (ABS(zcent(k)-(sim_targetYmin+sim_targetYmax)/2)<=(sim_targetYmax-sim_targetYmin)/2 .or. &
                       ABS(zcent(k)-(sim_targetYmin+sim_targetYmax)/2-sim_vacuumHeight)<=(sim_targetYmax-sim_targetYmin)/2))then
                    species = TARG_SPEC
                end if
                if (  (ABS(ycent(j)-0.03)<=(sim_targetYmax-sim_targetYmin)/2 .or. &
                       ABS(ycent(j)+0.03)<=(sim_targetYmax-sim_targetYmin)/2)&
                        .and. &
                       (xcent(i)-0.05-sim_targetRadius)**2+(zcent(k)-(sim_targetYmin+sim_vacuumHeight+sim_targetYmax)/2)**2 <= &
                       (sim_targetYmax+sim_vacuumHeight-sim_targetYmin)**2/4 .and. &
                       (xcent(i)-0.05-sim_targetRadius)**2+(zcent(k)-(sim_targetYmin+sim_vacuumHeight+sim_targetYmax)/2)**2 >= &
                       (-sim_targetYmax+sim_vacuumHeight+sim_targetYmin)**2/4 .and. &
                       xcent(i)>= (0.05+sim_targetRadius)                )then
                       species = TARG_SPEC
                end if

              end if
           else !!sphere
               if (sqrt(xcent(i)**2+ycent(j)**2+zcent(k)**2)<= sim_targetRadius) then
                  species = TARG_SPEC
              end if
           end if
           if(species == TARG_SPEC) then
              rho = sim_rhoTarg
              tele = sim_teleTarg
              tion = sim_tionTarg
              trad = sim_tradTarg
#ifdef TAR2_SPEC
           else if(species == TAR2_SPEC) then
                 rho = sim_rhoTar2
                 tele = sim_teleTar2
                 tion = sim_tionTar2
                 trad = sim_tradTar2
                ! write(*,*) xcent(i),ycent(i),sim_target2a1*ABS(xcent(i))+sim_target2b1,i,j,k
              !   write(*,*)sim_rhoTar2
#endif
#ifdef TAR3_SPEC
           else if(species == TAR3_SPEC) then
                  rho = sim_rhoTar3
                  tele = sim_teleTar3
                  tion = sim_tionTar3
                  trad = sim_tradTar3
#endif
#ifdef CONE_SPEC
           else if(species == CONE_SPEC) then
                      rho = sim_rhoCone
                      tele = sim_teleCone
                      tion = sim_tionCone
                      trad = sim_tradCone
#endif
#ifdef WALL_SPEC
           else if(species==WALL_SPEC) then  !!added
              rho = sim_rhoWall
              tele = sim_teleWall
              tion = sim_tionWall
              trad = sim_tradWall
#endif
           else
              rho = sim_rhoCham
              tele = sim_teleCham
              tion = sim_tionCham
              trad = sim_tradCham

           end if
          ! write(*,*) "check point 3",xcent(i),ycent(i),species, rho

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tele)

#ifdef FLASH_3T
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)
           !call Grid_putPointData(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, bx)
           !call Grid_putPointData(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, by)

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
              if (NDIM == 2) then
                 facexData(MAG_FACE_VAR,i,j,k)=0 ! bx
                 faceyData(MAG_FACE_VAR,i,j,k)=0 ! by
                 if (NDIM>2) facezData(MAG_FACE_VAR,i,j,k)= 0.0
              endif
           endif
#endif
!end if
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
