pro sos_aper, model, nbin, sos_aper1

plotsym, 0, /fill, 0.5

; program to take model PSF and use SNR to determine best aperture, or perhaps 3 apertures?
; 
; Called by psf_aper
; 
; Output: psf_aper[*,*,3] fo r3 apertures where *,* are dimensions of ccd raster
; 
; get total number of binned pixels in the frame
npixel=n_elements(model)

;reorder pixels by data number
sb=size(model, /dim)            ; size of binned array

mod1d=reform(model, sb[0]*sb[1])
sort1=reverse(sort(mod1d))
mod1d=mod1d[sort1]

count=2 ; start on count > 0
sig1=total(mod1d[0:1])*3.25
noise1=sqrt(mod1d[1] + 1*10. + (15.^2) + 30.)
snr=sig1/noise1

; plot, [count], [snr], color=cgcolor('black'), psym=8, xrange=[0, npixel], yrange=[0,20000]

repeat begin

  count=count+1
  snr1=snr
  
  sig1=total(mod1d[0:count])*3.25
  
  noise1=sqrt(mod1d[count] + count*10. + (15.^2) + 30.)
  
  snr=sig1/noise1
  
;   oplot, [count], [snr], color=cgcolor('black'), psym=8
  
endrep until snr lt snr1
; stop
xx=count-1          ; initial aperture
psfp=sort1[0:xx]
xi=(array_indices(model, psfp))[0,*]
yi=(array_indices(model, psfp))[1,*]

;wset, 0
;plot_image, bytscl(model, 0, 500)
;wset, 1
;plot_image, bytscl(model, 0, 500)
;plotsym, 0, 0.1
;oplot, xi, yi, color=cgcolor('purple'), psym=8

sos_aper1=model*0
sos_aper1[xi,yi]=model[xi,yi]

;plot_image, bytscl(sos_aper1, 0, 500)
;xx=where(sos_aper1 ne 0, nxx)
;xi=(array_indices(model, xx))[0,*]
;yi=(array_indices(model, xx))[1,*]
;oplot, xi, yi, color=cgcolor('purple'), psym=8
; 
; 
end