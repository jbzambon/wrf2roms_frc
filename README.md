# wrf2roms_frc

Programs to convert WRF output to ROMS forcing files.  Need to run wrfout_to_cf.ncl first on WRF output and then the wrf2roms_frc.m file.  Default selections in both should work.

Need up fix bug with NCL running out of memory.  Difficult to use on less than 4GB of RAM.

Example usage...
```
conda create -n ncl_stable -c conda-forge ncl
conda activate ncl
ncl 'file_in="wrfout_day1.nc"' 'file_out="wrfout_day1_post.nc"' wrfout_to_cf.ncl
matlab wrf2roms_frc.m
```

Joseph B. Zambon  
jbzambon@ncsu.edu  
1 June 2021
