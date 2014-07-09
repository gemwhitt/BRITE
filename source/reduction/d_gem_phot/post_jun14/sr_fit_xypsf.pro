pro sr_fit_xypsf, sdat,iloc,g, xy_psf

  ; Slavek fitting programs - use to get PSF centers
  ;
  ; sdat=subset of data - one per orbit
  
  nfrm=n_elements(iloc)
  xsize0=(size(sdat, /dim))[0]
  ysize0=(size(sdat, /dim))[1]
  
  ; get an average PSF to determine approx. center
  x = fltarr(xsize0,ysize0)                     ; average image
  m = nfrm                                      ; how many images
  for i=0,m-1 do x = x + sdat[*,*,i]            ; adding images
  x = x/m                                       ; average image
  
  max_pix=where(x eq max(x), nmax)
  if nmax gt 1 then max_pix=max_pix[0]
  
  pks=[array_indices(x, max_pix), max(x)]
  
  rbin=10.
  
  ; rebin the images
  xb=rebin_data(x,rbin) ; binned 10x10  ; binned average
  sdat2=rebin_data(sdat, rbin)    ; binned 10x10  ; binned images
  
  pks2=pks                      ; pks in binned data images and average
  pks2[0]=(pks[0])*rbin
  pks2[1]=(pks[1])*rbin
  
  xsize=xsize0*rbin ; binned size
  ysize=ysize0*rbin ; binned size
  
  ; Use a Gaussian as the 0-th approximation for the PSF...
  cen=max(xb)
  sig=2.*rbin
  psf0 = first_approx(cen,sig,rbin,pks2,xsize,ysize)
  
  ; Determine the shifts in small pixels from the image centre for individual images assuming a given PSF.
  ; note that these shifts are given in "small pixels", so shx/8. is the shift in big pixels
  get_shifts, sdat2,psf0,rbin,xsize,ysize, shx1,shy1
  
  count=-1
  repeat begin
  
    count=count+1
    shx=shx1
    shy=shy1
    
    ; Using the shx, shy, adjust individual images to the centre. The determine the improved 2-D PSF.
    ; Note the background is subtracted (for use as an PSF), but is given by "bck".
    psf1 = get_psf(sdat2,shx,shy,rbin,xsize,ysize, dud,bck,w)  ; w is 2d-indices of background pixels
    
    get_shifts, sdat2,psf1,rbin,xsize,ysize, shx1,shy1
    
    ; compute the residuals between shifts of old fit and shifts of new fit
    shx_res=abs(shx1-shx)
    shy_res=abs(shy1-shy)
    
    ; determine if there are any shifts and if so, repeat
    rs1=where(shx_res gt 0, nrs1)
    rs2=where(shy_res gt 0, nrs2)
    
    if nrs1 gt 0 OR nrs2 gt 0 then cycle=1 else cycle=0
    
  endrep until cycle eq 0  ; step 7 is to iterate until convergence
  
  psf=psf1  ; binned 10x10
  
  ; save PSF center locations
  xy_psf[0,iloc]=pks2[0]+shx
  xy_psf[1,iloc]=pks2[1]+shy
  
  endoffitting:
end
