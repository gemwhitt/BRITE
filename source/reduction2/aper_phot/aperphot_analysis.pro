pro aperphot_analysis

; Use to plot results of aper_phot and investigate the variables with the original images and apertures
; 
Compile_opt idl2

nbin=10

targets=['HD127973','HD129056','HD135379','HD136415']

txtdir='~/BRITE/TOR/CENTAURUS/data/aper_lc/'

savdir='~/BRITE/TOR/CENTAURUS/data/p5/'

filein=file_search(txtdir+'*.txt', count=nf)

for f=0, nf-1 do begin
  
readcol, filein[f], time, im, flux1, flux2, bkgd, bkgd_err, xcom, ycom, temperature, vmag, bmag, resid1, resid2, medimg, medcol, $
  format='d,i,d,d,f,f,f,f,f,f,f,d,d,f,f'
  
  print, file_basename(filein[f], '.txt')
  print, vmag[0]
  fname=file_basename(filein[f], '.txt')
  
; restore data1 and psf_aper
obj=obj_new('IDL_Savefile', savdir+fname+'*.sav')
obj->restore, 'data1'
obj->restore, 'psf_aper'  
  
time1=time-time[0]

pcflux1=flux1/(flux1+resid1)*100.
pcflux2=flux2/(flux2+resid2)*100.

pcdiff=pcflux2-pcflux1

npt=n_elements(time)

plotsym, 0, /fill, 0.8
wset, 0
plot, flux2, color=cgcolor('black'), psym=8, /nodata, yrange=[min(flux1),max(flux2)], xrange=[0,150]
oplot, flux1, color=cgcolor('orange'), psym=8
oplot, flux2, color=cgcolor('purple'), psym=8

wset, 1
plot, pcflux1, color=cgcolor('black'), psym=8, xrange=[0,150], yrange=[70,100], /nodata
oplot, pcflux1, color=cgcolor('orange'), psym=8
oplot, pcflux2, color=cgcolor('purple'), psym=8

; count number of pixels in each aperture from psf_aper
aper1=where(psf_aper[*,*,0] ne 0, nap1)
aper2=where(psf_aper[*,*,1] ne 0, nap2)

wset, 2
plot,  medimg, color=cgcolor('black'), psym=8, xrange=[0,150], /ynozero
oplot, medcol, color=cgcolor('blue'), psym=8


stop

endfor


stop
for i=0, 100 do begin
  

  im0=data1[*,*,im[i]]
  
  s=size(im0, /dim)
  
  im1=congrid(im0, nbin*s[0], nbin*s[1], /interp, /center)
  
  ; now modify xy_com to obtain the center of PSF in this binned image
  xcen=fix(round(xcom[im[i]]*nbin))
  ycen=fix(round(ycom[im[i]]*nbin))
  
  ; get the center of the binned image
  xcen2=s[0]*nbin/2.
  ycen2=s[1]*nbin/2.
  
  ; calculate the shifts
  xdiff=xcen2-xcen
  ydiff=ycen2-ycen
  
  ap=0
  
  modelapp=reform(psf_aper[*,*,ap])
    
  aper=modelapp*0
    
  ; find the locations of ALL PSF pixels gt 0
  psfpix=where(modelapp gt 0, npsf)
  xpix=(array_indices(modelapp,psfpix))[0,*]
  ypix=(array_indices(modelapp,psfpix))[1,*]
    
  aper[xpix-xdiff,ypix-ydiff]=modelapp[xpix,ypix]
  
  modpix=where(aper gt 0, npsf)
  xpix=(array_indices(aper,modpix))[0,*]
  ypix=(array_indices(aper,modpix))[1,*]
  
  wset, 2  
  plot_image, bytscl(im1, 0, 500), title=im[i], charsize=0.8, color=cgcolor('black')
  
  oplot, xpix, ypix, color=cgcolor('orange'), psym=2
  
  stop  
  
endfor

stop



print, 'End of program'
end