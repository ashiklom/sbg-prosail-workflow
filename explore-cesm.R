library(ncdf4)

bfile <- "data/cesm/b.e11.B1850C5CN.f09_g16.005.clm2.h0.LAISHA.200001-209912.nc"
ffile <- "data/cesm/f.e11.F1850C5CN.f09_f09.001.clm2.h0.LAISHA.200001-209912.nc"

bnc <- nc_open(bfile)
