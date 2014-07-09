pro get_bkgd, dat,pks,rbin,  bkgd,bk_err

; get dimensions of roi in sub-pixels
xdim=(size(dat, /dim))[0]
ydim=(size(dat, /dim))[1]

; get x and y locations of psf
xc=pks[0]/float(rbin)
yc=pks[1]/float(rbin)

; check xc and yc are far from the image edges - if not then reselect the background region
roi_xc=float(xdim)/2.
roi_yc=float(ydim)/2.

; get distance between pks and roi center
xdif=pks[0]-roi_xc
ydif=pks[1]-roi_yc

if abs(xdif) lt 8.*rbin AND abs(ydif) lt 8.*rbin then begin ; target is roughly in center
  ; define a region around the perimeter
  
  temp_dat=dat
  temp_dat[2*rbin:xdim-(2*rbin)-1,2*rbin:ydim-(2*rbin)-1]=-9999
  
  xx=where(temp_dat ne -9999, nxx)
  
  x_loc=(array_indices(dat, xx))[0,*]
  y_loc=(array_indices(dat, xx))[1,*]
  
;  plot_image, dat
;  oplot, x_loc, y_loc, color=cgcolor('orange'), psym=2
  
;  stop
  
endif else begin  ; target is within a quadrent
  ; define a region based in the other 3 quads
  
  temp_dat=dat
  
  ; define 4 quads - lower-left, right and upper-left, right
  if xdif lt 0. and ydif lt 0. then temp_dat[0:xdim-(4*rbin)-1,0:ydim-(4*rbin)-1]=-9999 ; lower-left
  if xdif gt 0. and ydif gt 0. then temp_dat[4*rbin:xdim-1, 0:ydim-(4*rbin)-1]=-9999  ; lower-right
  if xdif lt 0. and ydif gt 0. then temp_dat[0:xdim-(4*rbin)-1,4*rbin:ydim-1]=-9999  ; upper-left
  if xdif gt 0. and ydif gt 0. then temp_dat[4*rbin:xdim-1,4*rbin:ydim-1]=-9999  ; upper-right
  
  xx=where(temp_dat ne -9999, nxx)
  
  x_loc=(array_indices(dat, xx))[0,*]
  y_loc=(array_indices(dat, xx))[1,*]
  
;  plot_image, dat
;  oplot, x_loc, y_loc, color=cgcolor('orange'), psym=2
  
;  stop
  
endelse


bkpix=[dat[xx]]

bkgd=robust_mean(bkpix,3)

bk1=bkpix/bkgd

bad=where(finite(bk1) eq 0, nbad)

if nbad gt 0 OR n_elements(bkpix) lt 50 then bk_err=9.9999 else bk_err=robust_sigma(bk1) 

end