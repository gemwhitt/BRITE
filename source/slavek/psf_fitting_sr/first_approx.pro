function first_approx, cen,sig,nbins,pks,xsize,ysize
  ; SMR Oct.2013
  ;
  ; prepares centred Gaussian PSF over nsmall_pixels x nsmall_pixels
  ; input:
  ;   cen = central strength, typically 1.0 or some value
  ;   sig = Gaussian sigma, FWHM = 2.354 sigma
  ; output:
  ;   psf0 = Gusssian PSF
  ; uses:
  ;   gs2_smooth.pro
  ; run:
  ;   psf0 = B2_im3(8000.,1.5, nbins)
  
  x = fltarr(xsize,ysize)
  x = shift(gs2_smooth(xsize,sig),pks[0],pks[1])
  x = x/max(x)
  x = cen*x
  
  return,x
end

