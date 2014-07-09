pro psf_fwhm, psf, iloc, sdat, nfrm, fwhm_psf

  ; rebin psf back to original bin size
  psf=rebin(psf1, xsize1, ysize1)
  
  ;use the PSF to get an average flux value and pixel count
  xsize=(size(psf, /dim))[0]
  ysize=(size(psf, /dim))[1]
  psf1d=reform(psf, xsize1*ysize1)
  ord1=reverse(sort(psf1d))
  signal1=psf1d[ord1]
  
  fwhm_avg=fullwid_halfmax(psf) ; average FWHM for the average PSF
  
  for i=0, nfrm-1 do fwhm_psf[*,iloc[i]]=fullwid_halfmax(reform(sdat[x1:x2,y1:y2,i]), peak_index=pi) ; FWHM for each image
  
  if fwhm_avg[0] eq 0 or fwhm_avg[1] eq 0 then stop
  
  for i=0, nfrm-1 do begin
  
    plot_image, bytscl(sdat[*,*,i], 20, 1500)
    
    pkimg=array_indices(sdat[*,*,i], where(sdat[*,*,i] eq max(sdat[*,*,i])))
    
    oplot, [xy_psf[0,i]-fwhm_psf[0,iloc[i]]/2.,xy_psf[0,i]+fwhm_psf[0,iloc[i]]/2., $
      xy_psf[0,i]+fwhm_psf[0,iloc[i]]/2.,xy_psf[0,i]-fwhm_psf[0,iloc[i]]/2.,xy_psf[0,i]-fwhm_psf[0,iloc[i]]/2.], $
      [xy_psf[1,i]-fwhm_psf[1,iloc[i]]/2.,xy_psf[1,i]-fwhm_psf[1,iloc[i]]/2.,$
      xy_psf[1,i]+fwhm_psf[1,iloc[i]]/2.,xy_psf[1,i]+fwhm_psf[1,iloc[i]]/2.,xy_psf[1,i]-fwhm_psf[1,iloc[i]]/2.], $
      color=cgcolor('blue'), thick=3
    oplot, [xy_psf[0,i]], [xy_psf[1,i]], color=cgcolor('red'), psym=2
    
    stop
    
    
  endfor
  
  minpix1=round(fwhm_psf[0,*]*fwhm_psf[1,*])
  maxpix1=round(2.5*fwhm_xy[0]*2.5*fwhm_xy[1])
  stop
  minflux1=total(signal1[0:minpix1])
  stop
  if maxpix1 gt (xsize*ysize)-1 then maxpix1=(xsize*ysize)-1
  maxflux1=total(signal1[0:maxpix1])
  stop
  sigcum=lonarr(maxpix1)
  
  for i=0, maxpix1-1 do sigcum[i]=total(signal1[0:i])
  
  xx=(where(sigcum ge 0.9*maxflux1))[0]
  
  mostflux=0.9*maxflux1
  mostpix=xx
  stop
  npix_psf[*,iloc]=[minpix1,mostpix,maxpix1]
  fl_psf[*,iloc]=[minflux1,mostflux,maxflux1]


end