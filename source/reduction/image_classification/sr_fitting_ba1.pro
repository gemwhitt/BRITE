pro sr_fitting_ba1, data2,iloc,orbit,outdir,hdname,rbin, xy_psf,mod1,npix1

; Slavek fitting programs - use to get PSF centers
;
; call by modelpsf_ba1
  
nfrm=n_elements(iloc)

s0=size(data2, /dim)
  
; get an average PSF AND determine approx. center
x = fltarr(s0[0],s0[1])                     ; average image
m = nfrm                                      ; how many images
for i=0,m-1 do x = x + data2[*,*,i]            ; adding images
x = x/m                                       ; average image

xtot=total(x)

xcen=total(total(x, 2)*indgen(s0[0]))/xtot
ycen=total(total(x, 1)*indgen(s0[1]))/xtot
  
pks=[xcen, ycen]

; rebin the images
xb=rebin_data(x,rbin)           ; binned average
dat2=rebin_data(data2, rbin)    ; binned images
  
pks2=pks                      ; pks in binned data images and average
pks2[0]=(pks[0])*rbin
pks2[1]=(pks[1])*rbin
  
s1=size(xb, /dim)
  
; Use a Gaussian as the 0-th approximation for the PSF...
cen=max(xb)
sig=2.*rbin
psf0 = first_approx(cen,sig,rbin,pks2,s1[0],s1[1])

; Determine the shifts in small pixels from the image centre for individual images assuming a given PSF.
; note that these shifts are given in "small pixels", so shx/8. is the shift in big pixels
get_shifts, dat2,psf0,rbin,s1[0],s1[1], shx1,shy1
stop  
  count=-1
  repeat begin
  
    count=count+1
    shx=shx1
    shy=shy1
    
    ; Using the shx, shy, adjust individual images to the centre. The determine the improved 2-D PSF.
    ; Note the background is subtracted (for use as an PSF), but is given by "bck".
    psf1 = get_psf(sdat2,shx,shy,rbin,s1[0],s1[1], dud,bck,w)  ; w is 2d-indices of background pixels
    
    get_shifts, dat2,psf1,rbin,s1[0],s1[1], shx1,shy1
    
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
  
  
  ; get characteristics of PSF in each image
  data1d=reform(psf, s1[0]*s1[1])
  sort1=reverse(sort(data1d))
  data1d=data1d[sort1]/float(rbin^2)
  
  count=2
  sig1=(data1d[0:1])*3.25
  snr=0
  
  repeat begin
  
    count=count+1
    snr1=snr
    
    sig1=total(data1d[0:count])*3.25
    
    noise1=sqrt(sig1 + count*(1/float(rbin^2) + (15^2)/float(rbin^2) + 21/float(rbin^2)))
    
    snr=sig1/noise1
    
  endrep until snr lt snr1
  
  xx=count
  psfp=sort1[0:xx]
  xflx=(array_indices(psf, psfp))[0,*]
  yflx=(array_indices(psf, psfp))[1,*]
  
  modelpsf=lonarr(xsize0*rbin,ysize0*rbin)
  xc=xsize0*rbin/2.
  yc=ysize0*rbin/2.
  modelpsf[xflx+(xc-pks2[0]), yflx+(yc-pks2[1])]=psf[xflx, yflx]
  
  mod1[*,*,orbit]=modelpsf
  
  ;plot_image, bytscl(modelpsf, 0, 1)
  
  npix1[orbit]=count
  
  
  endoffitting:
end