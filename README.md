# FLASH-EXAMPLES
Contains several examples for FLASH4.6.2 code

- folder `/cn4files` contains several .cn4 files may be used
- folder `/wjz_mhd3d` is an example of 3-D MHD or 3-D HD simulation with proton imaging module. 

```
./setup -auto wjz_mhd3d -3d +hdf5typeio -nxb=32 -nyb=32 -nzb=32 species=cham,targ,tar2,tar3 +mtmmmt +laser +usm3t +mgd mgd_meshgroups=6 -parfile=flash.par -objdir=astrowjz/wjz_mhd3d_data +protonImaging pi_maxBeams=1 pi_maxDetectors=1 threadProtonTrace=True`
```
- folder `/wjz_mhd2d` is an example of 2-D MHD or 2-D HD simulation with proton imaging module. 
```
./setup -auto wjz_mhd2d -2d +hdf5typeio -nxb=32 -nyb=32 -nzb=32 species=cham,targ +mtmmmt +laser +uhd3t +mgd mgd_meshgroups=6 -parfile=flash.par -objdir=astrowjz/wjz_mhd3d_data
```

- folder `/wjz_protonimaging` is an example of proton imaging module. 
```
./setup wjz_protonimaging -auto -3d -geometry=cartesian -parfile=flash.par -maxblocks=4000 -nxb=16 -nyb=16 -nzb=16 +protonImaging --without-unit=physics/Hydro --without-unit=physics/Eos -objdir=astrowjz/proton16 pi_maxBeams=1 pi_maxDetectors=1 threadProtonTrace=True
 ```
