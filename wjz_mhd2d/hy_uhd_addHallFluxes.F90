!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_addHallFluxes
!!
!! NAME
!!
!!  hy_uhd_addHallFluxes
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_addHallFluxes(blockID,blkLimitsGC,ix,iy,iz,sweepDir)
!!
!!  hy_uhd_addHallFluxes(integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM),
!!                            integer(IN) :: ix,
!!                            integer(IN) :: iy,
!!                            integer(IN) :: iz,
!!                            real(INOUT)    :: Flux,
!!                            integer(IN) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds hall flux contributions to total MHD fluxes
!!
!!
!! ARGUMENTS
!!
!!  blockID     - a local blockID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Flux        - array containing MHD fluxes
!!  sweepDir    - direction of sweep
!!
!!***

!!REORDER(4): U

Subroutine hy_uhd_addHallFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,sweepDir)

  use Grid_interface, ONLY : Grid_getBlkPtr, &
                             Grid_releaseBlkPtr, &
                             Grid_getDeltas, &
                             Grid_getCellCoords

  use hy_uhd_slopeLimiters, ONLY : mc, &
                                   minmod, &
                                   vanLeer, &
                                   signum, &
                                   get_upwind

  use Hydro_data, ONLY: hy_geometry, &
                        hy_conserveAngField, &
                        hy_avogadro, &
                        hy_qele, &
                        hy_useHall

  use MagneticResistivity_interface, &
                  ONLY : MagneticResistivity
  use Eos_interface,  ONLY : Eos_getAbarZbar

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Argument List ----------------------------------------------------------
  integer, INTENT(IN) :: blockID,ix,iy,iz
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
  real, dimension(HY_VARINUM), intent(INOUT) :: Flux
  integer, INTENT(IN) :: sweepDir

    !! ----------------------------------------------------------------------
  real :: dxBy, dxBz, dyBx, dyBz, dzBx, dzBy
  real :: drBz, drBphi
  real :: Ex, Ey, Ez, Bx, By, Bz,jx,jy,jz,idx,idy,idz
  real :: jr, jphi, Er, Ephi, Br, Bphi
  real :: jz_dl, jz_dr, jz_ul, jz_ur
  real :: ye, nele_left,nele_right,nele, inv_r, sumY
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: U !,bx,by
  real :: inv_dVr !2/(rp*abs(rp) - rm*abs(rm))
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: xCenter, yCenter
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS)) :: xCenter, yCenter

#endif

  if(.not. hy_useHall) return

  if (hy_useHall .and. (hy_geometry == CYLINDRICAL) .and. (.not.hy_conserveAngField)) then
  print *, "[hy_uhd_addHallFluxes]: Hall works only with hy_conserveAngField=TRUE. Hall effect turned off."
    return
  end if


  !! Get deltas
  call Grid_getDeltas(blockID,del)

  idx=1./del(DIR_X)
  if (NDIM >= 2) then
     idy=1./del(DIR_Y)
     if (NDIM == 3) then
        idz=1./del(DIR_Z)
     endif
  endif

  !! Get pointer
  call Grid_getBlkPtr(blockID,U,CENTER)

  !! Compute Electron density at (i,j,k)
  !if (.not. oneTemp) then
     call Eos_getAbarZbar(U(:,ix,iy,iz), Ye=ye, sumy=sumY)
     nele_right = ye * hy_avogadro * U(DENS_VAR,ix,iy,iz)
  !else
    ! nele_right = U(DENS_VAR,ix,iy,iz)
  !endif

  !! Get cell x-coords for this block
 ! if (hy_geometry /= CARTESIAN) then
   call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
   call Grid_getCellCoords(JAXIS,blockID, CENTER,.true.,yCenter, blkLimitsGC(HIGH,JAXIS))

  !endif

  if (hy_geometry == CARTESIAN) then

  !! Compute Hall parts and add them to flux components
     select case(sweepDir)

     case(DIR_X)
     !! Electron density on the face
      	!if (.not. oneTemp) then
         call Eos_getAbarZbar(U(:,ix-1,iy,iz), Ye=ye)
         nele_left = ye * hy_avogadro * U(DENS_VAR,ix-1,iy,iz)
      	!else
         !nele_left = U(DENS_VAR,ix-1,iy,iz)
      	!endif
        nele = (nele_left + nele_right)*0.5

     !! face-centered B-field
        Bx = (U(MAGX_VAR,ix,iy,iz) + U(MAGX_VAR,ix-1,iy,iz))*0.5
        By = (U(MAGY_VAR,ix,iy,iz) + U(MAGY_VAR,ix-1,iy,iz))*0.5
        Bz = (U(MAGZ_VAR,ix,iy,iz) + U(MAGZ_VAR,ix-1,iy,iz))*0.5

     !! 1D case : d/dy=d/dz=0
        jx = 0.0
     !! jy = -dBz/dx
        jy = -(U(MAGZ_VAR,ix,iy,iz)-U(MAGZ_VAR,ix-1,iy,iz))*idx
     !! jz = dBy/dx
        jz =  (U(MAGY_VAR,ix,iy,iz)-U(MAGY_VAR,ix-1,iy,iz))*idx
#if NDIM >= 2
     !! 2D case : d/dy .ne. 0 but d/dz=0
     !! jx = dBz/dy
        jx = (U(MAGZ_VAR,ix,  iy+1,iz) - U(MAGZ_VAR,ix,  iy-1,iz) &
                +  U(MAGZ_VAR,ix-1,iy+1,iz) - U(MAGZ_VAR,ix-1,iy-1,iz))*0.25*idy
     !! jz = dBy/dx - dBx/dy
        jz = jz - (U(MAGX_VAR,ix,  iy+1,iz) - U(MAGX_VAR,ix,  iy-1,iz) &
                +  U(MAGX_VAR,ix-1,iy+1,iz) - U(MAGX_VAR,ix-1,iy-1,iz))*0.25*idy
#if NDIM == 3
     !! jx = dBz/dy - dBy/dz
        jx = jx - (U(MAGY_VAR, ix,  iy,iz+1) - U(MAGY_VAR, ix,  iy, iz-1)&
                +  U(MAGY_VAR, ix-1,iy,iz+1) - U(MAGY_VAR, ix-1,iy, iz-1))*0.25*idz
     !! jy = dBx/dz - dBz/dx
        jy = jy + (U(MAGX_VAR, ix,  iy,iz+1) - U(MAGX_VAR, ix,  iy, iz-1)&
                +  U(MAGX_VAR, ix-1,iy,iz+1) - U(MAGX_VAR, ix-1,iy, iz-1))*0.25*idz


#endif
#endif

        !! Hall electric field
       ! Ey = 0!i(jz*Bx - jx*Bz)/(hy_qele*nele)
        if (ycenter(iy)>0.14 .and.((U(DENS_VAR,ix,iy,iz)+U(DENS_VAR,ix-1,iy,iz)) &
        +(U(DENS_VAR,ix-1,iy+1,iz)+U(DENS_VAR,ix-1,iy-1,iz)+U(DENS_VAR,ix,iy+1,iz)+U(DENS_VAR,ix,iy-1,iz))/2)>0.2 )  then
             Ez =  5e12 !(jx*By - jy*Bx)/(hy_qele*nele)
             Ey = 0
        else 
             Ez = (jx*By - jy*Bx)/(hy_qele*nele)
             Ey = (jz*Bx - jx*Bz)/(hy_qele*nele)
        endif 

        !! Hall contributions to the magnetic x-flux
        Flux(F07MAGY_FLUX) = Flux(F07MAGY_FLUX) - Ez
        Flux(F08MAGZ_FLUX) = Flux(F08MAGZ_FLUX) + Ey

        !! Hall x-component of the Poynting vector
        Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) + Ey*Bz - Ez*By


#if NDIM >= 2
     case(DIR_Y)
     !! Electron density on the face
      	!if (.not. oneTemp) then
         call Eos_getAbarZbar(U(:,ix,iy-1,iz), Ye=ye)
         nele_left = ye * hy_avogadro * U(DENS_VAR,ix,iy-1,iz)
      	!else
         !nele_left = U(DENS_VAR,ix,iy-1,iz)
      	!endif
        nele = (nele_left + nele_right)*0.5

     !! face-centered B-field
        Bx = (U(MAGX_VAR,ix,iy,iz) + U(MAGX_VAR,ix,iy-1,iz))*0.5
        By = (U(MAGY_VAR,ix,iy,iz) + U(MAGY_VAR,ix,iy-1,iz))*0.5
        Bz = (U(MAGZ_VAR,ix,iy,iz) + U(MAGZ_VAR,ix,iy-1,iz))*0.5

     !! jx = dBz/dy
        jx = (U(MAGZ_VAR,ix,iy,iz) - U(MAGZ_VAR,ix,iy-1,iz))*idy

     !! jy = -dBz/dx
        jy = -(U(MAGZ_VAR,ix+1,iy,iz) - U(MAGZ_VAR,ix-1,iy,iz) &
                +  U(MAGZ_VAR,ix+1,iy-1,iz) - U(MAGZ_VAR,ix-1,iy-1,iz))*0.25*idx

     !! jz = dBy/dx - dBx/dy
        jz = - (U(MAGX_VAR,ix  ,iy  ,iz ) - U(MAGX_VAR,ix  ,iy-1,iz))*idy &
             + (U(MAGY_VAR,ix+1,iy  ,iz ) - U(MAGY_VAR,ix-1,iy  ,iz)      &
             +  U(MAGY_VAR,ix+1,iy-1,iz ) - U(MAGY_VAR,ix-1,iy-1,iz))*0.25*idx



#if NDIM == 3
     !! jx = dBz/dy - dBy/dz
        jx = jx - (U(MAGY_VAR,ix,iy,iz+1 ) - U(MAGY_VAR,ix,iy,iz-1)        &
              +    U(MAGY_VAR,ix,iy-1,iz+1 ) - U(MAGY_VAR,ix,iy-1,iz-1))*0.25*idz   !!! ATTENTION: dans hy_uhd_addResistiveFluxes cetait idx au lieu de idz

     !! jy = dBx/dz - dBz/dx
        jy = jy + (U(MAGX_VAR,ix,iy,iz+1 ) - U(MAGX_VAR,ix,iy,iz-1)        &
              +    U(MAGX_VAR,ix,iy-1,iz+1 ) - U(MAGX_VAR,ix,iy-1,iz-1))*0.25*idz

#endif

        !! Hall electric field
        !Ex =0! (jy*Bz - jz*By)/(hy_qele*nele)
        if (ycenter(iy)>0.14 .and.((U(DENS_VAR,ix,iy,iz)+U(DENS_VAR,ix,iy-1,iz)) &
        +(U(DENS_VAR,ix+1,iy,iz)+U(DENS_VAR,ix-1,iy,iz)+U(DENS_VAR,ix+1,iy-1,iz)+U(DENS_VAR,ix-1,iy-1,iz))/2)>0.2 ) then
             Ez = 5e12!2e12!(jx*By - jy*Bx)/(hy_qele*nele)
             Ex = 0
        else
             Ex = (jy*Bz - jz*By)/(hy_qele*nele)
             Ez = (jx*By - jy*Bx)/(hy_qele*nele)
        endif

        !! Hall contributions to the magnetic y-flux
        Flux(F06MAGX_FLUX) = Flux(F06MAGX_FLUX) + Ez
        Flux(F08MAGZ_FLUX) = Flux(F08MAGZ_FLUX) - Ex

        !! Hall y-component of the Poynting vector
        Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) + Ez*Bx - Ex*Bz


#if NDIM == 3
     case(DIR_Z)
    !! Electron density on the face
      	!if (.not. oneTemp) then
         call Eos_getAbarZbar(U(:,ix,iy,iz-1), Ye=ye)
         nele_left = ye * hy_avogadro * U(DENS_VAR,ix,iy,iz-1)
      	!else
         !nele_left = U(DENS_VAR,ix,iy,iz-1)
      	!endif
        nele = (nele_left + nele_right)*0.5

     !! face-centered B-field
        Bx = (U(MAGX_VAR,ix,iy,iz) + U(MAGX_VAR,ix,iy,iz-1))*0.5
        By = (U(MAGY_VAR,ix,iy,iz) + U(MAGY_VAR,ix,iy,iz-1))*0.5
        Bz = (U(MAGZ_VAR,ix,iy,iz) + U(MAGZ_VAR,ix,iy,iz-1))*0.5

     !! jx = dBz/dy - dBy/dz
        jx = (U(MAGZ_VAR,ix,iy+1,iz  ) -U(MAGZ_VAR,ix,iy-1,iz  )           &
            + U(MAGZ_VAR,ix,iy+1,iz-1) -U(MAGZ_VAR,ix,iy-1,iz-1))*0.25*idy &
            -(U(MAGY_VAR,ix,iy,iz) -U(MAGY_VAR,ix,iy,iz-1))*idz

     !! jy = dBx/dz - dBz/dx
        jy = (U(MAGX_VAR,ix,iy,iz)-U(MAGX_VAR,ix,iy,iz-1))*idz &
            -(U(MAGZ_VAR,ix+1,iy,iz)-U(MAGZ_VAR,ix-1,iy,iz)      &
            + U(MAGZ_VAR,ix+1,iy,iz-1)-U(MAGZ_VAR,ix-1,iy,iz-1))*0.25*idx


     !! jz = dBy/dx - dBx/dy
        jz = (U(MAGY_VAR,ix+1,iy,iz)-U(MAGY_VAR,ix-1,iy,iz)      &
            + U(MAGY_VAR,ix+1,iy,iz-1)-U(MAGY_VAR,ix-1,iy,iz-1))*0.25*idx   &
            - (U(MAGX_VAR,ix,iy+1,iz  ) -U(MAGX_VAR,ix,iy-1,iz  )           &
            + U(MAGX_VAR,ix,iy+1,iz-1) -U(MAGX_VAR,ix,iy-1,iz-1))*0.25*idy

        !! Hall electric field
        Ex = (jy*Bz - jz*By)/(hy_qele*nele)
        Ey = (jz*Bx - jx*Bz)/(hy_qele*nele)

        !! Hall contributions to the magnetic z-flux
        Flux(F06MAGX_FLUX) = Flux(F06MAGX_FLUX) - Ey
        Flux(F07MAGY_FLUX) = Flux(F07MAGY_FLUX) + Ex

        !! Hall z-component of the Poynting vector
        Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) + Ex*By - Ey*Bx

#endif
#endif
     end select
  endif



   if (hy_geometry == CYLINDRICAL) then
!   !! Compute Hall parts and add them to flux components
!   !! Notice that X == R, Y == Z, Z == PHI. Be aware of signs
!   !! when calculating curls

    !! WARNING: "real" cylindrical coordinates below...
    !! For rflux, I need Ephi = jz*Br - jr*Bz and Ez = jr*Bphi - jphi*Br
    !! with:
    !! jr = 0 (1D) or -dBphi/dz (2D)
    !! jphi = - dBz/dr (1D) or dBr/dz - dBz/dr (2D)
    !! jz = 1/r*d(r*Bphi)/dr (1D)


  !! Compute Hall parts and add them to flux components
     select case(sweepDir)

     case(DIR_X)
     !! Electron density on the face
        !if (.not. oneTemp) then
         call Eos_getAbarZbar(U(:,ix-1,iy,iz), Ye=ye)
         nele_left = ye * hy_avogadro * U(DENS_VAR,ix-1,iy,iz)
        !else
         !nele_left = U(DENS_VAR,ix-1,iy,iz)
        !endif
        nele = (nele_left + nele_right)*0.5

    !! face-centered B-field
    Br = (U(MAGX_VAR,ix,iy,iz) + U(MAGX_VAR,ix-1,iy,iz))*0.5
    Bphi = (U(MAGZ_VAR,ix,iy,iz) + U(MAGZ_VAR,ix-1,iy,iz))*0.5
    Bz = (U(MAGY_VAR,ix,iy,iz) + U(MAGY_VAR,ix-1,iy,iz))*0.5

    !!! IF 1D...
    jr = 0.0

    jphi = -(U(MAGY_VAR,ix,iy,iz) - U(MAGY_VAR,ix-1,iy,iz))*idx

    inv_dVr = xCenter(ix)*xCenter(ix) - xCenter(ix-1)*abs(xCenter(ix-1))
    inv_dVr = 2.0/inv_dVr

    jz = (U(MAGZ_VAR,ix  ,iy,iz)*xCenter(ix) &
             -  U(MAGZ_VAR,ix-1,iy,iz)*abs(xCenter(ix-1)))*inv_dVr

#if NDIM >= 2
    jr = -(U(MAGZ_VAR,ix,  iy+1,iz) - U(MAGZ_VAR,ix,  iy-1,iz) &
                +  U(MAGZ_VAR,ix-1,iy+1,iz) - U(MAGZ_VAR,ix-1,iy-1,iz))*0.25*idy   !! TODO: check this, do I need to take care of the fact that Bphi changes sign

    jphi = jphi + (U(MAGX_VAR,ix,  iy+1,iz) - U(MAGX_VAR,ix,  iy-1,iz) &
                +  U(MAGX_VAR,ix-1,iy+1,iz) - U(MAGX_VAR,ix-1,iy-1,iz))*0.25*idy
#endif
    !!! endif NDIM >= 2

    if ((xCenter(ix)>0.10) .and. (U(DENS_VAR,ix,iy,iz)+U(DENS_VAR,ix-1,iy,iz))>1.0) then
      Ephi =  1.0E13!! add a potential as the driver
    else
      Ephi = 0
    endif
    Ez   = 0!(jr*Bphi - jphi*Br)/(hy_qele*nele)

    Flux(F07MAGY_FLUX) = Flux(F07MAGY_FLUX) + Ephi
    if (hy_conserveAngField) then
      Flux(F08MAGZ_FLUX) = Flux(F08MAGZ_FLUX) - Ez  !!this is still NOT added as a source if not conserveAngField (TODO ?)
    endif

    Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) + Ephi*Bz - Ez*Bphi


    case(DIR_Y)
     !! Electron density on the face
        !if (.not. oneTemp) then
         call Eos_getAbarZbar(U(:,ix,iy-1,iz), Ye=ye)
         nele_left = ye * hy_avogadro * U(DENS_VAR,ix,iy-1,iz)
        !else
         !nele_left = U(DENS_VAR,ix-1,iy,iz)
        !endif
        nele = (nele_left + nele_right)*0.5

    !! face-centered B-field
    Br = (U(MAGX_VAR,ix,iy,iz) + U(MAGX_VAR,ix,iy-1,iz))*0.5
    Bphi = (U(MAGZ_VAR,ix,iy,iz) + U(MAGZ_VAR,ix,iy-1,iz))*0.5
    Bz = (U(MAGY_VAR,ix,iy,iz) + U(MAGY_VAR,ix,iy-1,iz))*0.5

    !! jr = 0 (1D) or -dBphi/dz (2D)
    !! jphi = - dBz/dr (1D) or dBr/dz - dBz/dr (2D)
    !! jz = 1/r*d(r*Bphi)/dr (1D)

    jr = -(U(MAGZ_VAR,ix,iy,iz) - U(MAGZ_VAR,ix,iy-1,iz))*idy

    jphi = (U(MAGX_VAR,ix,iy,iz) - U(MAGX_VAR,ix,iy-1,iz))*idy
    jphi = jphi - (U(MAGY_VAR,ix+1,iy,iz) - U(MAGY_VAR,ix-1,iy,iz) &
                +  U(MAGY_VAR,ix+1,iy-1,iz) - U(MAGY_VAR,ix-1,iy-1,iz))*0.25*idx

    !!!! To compute the z-current at the z-interface with compute the four intermediate currents:
    !!!! Z-current in upper left
    inv_dVr = xCenter(ix)*xCenter(ix) - xCenter(ix-1)*abs(xCenter(ix-1))
    inv_dVr = 2.0/inv_dVr

    jz_ul = (U(MAGZ_VAR,ix  ,iy,iz)*xCenter(ix) &
             -  U(MAGZ_VAR,ix-1,iy,iz)*abs(xCenter(ix-1)))*inv_dVr

    !!!! Z-current in lower left (same inv_dVr than above)
    jz_dl = (U(MAGZ_VAR,ix  ,iy-1,iz)*xCenter(ix) &
             -  U(MAGZ_VAR,ix-1,iy-1,iz)*abs(xCenter(ix-1)))*inv_dVr

    !!!! Z-current in upper right
    inv_dVr = xCenter(ix+1)*xCenter(ix+1) - xCenter(ix)*abs(xCenter(ix))
    inv_dVr = 2.0/inv_dVr

    jz_ur = (U(MAGZ_VAR,ix+1  ,iy,iz)*xCenter(ix+1) &
             -  U(MAGZ_VAR,ix,iy,iz)*abs(xCenter(ix)))*inv_dVr

    !!!! Z-current in lower right (same inv_dVr than above)
    jz_dr = (U(MAGZ_VAR,ix+1  ,iy-1,iz)*xCenter(ix+1) &
             -  U(MAGZ_VAR,ix,iy-1,iz)*abs(xCenter(ix)))*inv_dVr

    !!!! Computation of the z-current on the correct interface
    jz = 0.25*(jz_dl + jz_dr + jz_ul + jz_ur)


    Er = 0!(jphi*Bz - jz*Bphi)/(hy_qele*nele)
    if ((xCenter(ix)>0.10) .and. (U(DENS_VAR,ix,iy,iz)+U(DENS_VAR,ix-1,iy,iz))>1.0) then
      Ephi = 1.0E13!! add a potential as the driver
    else
      Ephi =0
    endif


    Flux(F06MAGX_FLUX) = Flux(F06MAGX_FLUX) - Ephi
    Flux(F08MAGZ_FLUX) = Flux(F08MAGZ_FLUX) + Er

    Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) + Er*Bphi - Ephi*Br


     end select

   endif
      !!! endif hy_geometry == CYLINDRICAL


  if ((hy_geometry .ne. CYLINDRICAL) .and. (hy_geometry .ne. CARTESIAN)) then
      print *, &
        "[hy_uhd_addHallFluxes]: Hall effect only works for either Cartesian or Cylindrical coordinates!"
  endif

  !! Release pointer
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  End Subroutine hy_uhd_addHallFluxes
