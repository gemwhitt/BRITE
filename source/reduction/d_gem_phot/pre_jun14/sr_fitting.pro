pro sr_fitting, sdat,iloc,orb,outdir,hdname, xy_psf,npixel,avgfl,duds,fwhm_orb

; Slavek fitting programs - use to get PSF centers
; 
; call by get_psf_info, which is called by p3_p4_gem6
; 
; sdat=subset of data - one per orbit

nfrm=(size(sdat, /dim))[2]
xsize0=(size(sdat, /dim))[0]
ysize0=(size(sdat, /dim))[1]

; find approx center of PSF
x = fltarr(xsize0,ysize0)                             ; average image
m = nfrm                                      ; how many images
for i=0,m-1 do x = x + sdat[*,*,i]            ; adding images
x = x/m                                       ; average image

max_pix=where(x eq max(x), nmax)
if nmax gt 1 then max_pix=max_pix[0]

pks=[array_indices(x, max_pix), max(x)]

; make a cutout around maxpix
if pks[0]-10. lt 0. then x1=0. else x1=pks[0]-10. 
if pks[0]+10. gt xsize0-1 then x2=xsize0-1 else x2=pks[0]+10. 
if pks[1]-10. lt 0. then y1=0. else y1=pks[1]-10. 
if pks[1]+10. gt ysize0-1 then y2=ysize0-1 else y2=pks[1]+10.

cutout=x[x1:x2,y1:y2] 

; use this to rebin back later
xsize1=(size(cutout, /dim))[0]
ysize1=(size(cutout, /dim))[1]

rbin=10.
pks2=pks
pks2[0]=(pks[0]-x1)*rbin
pks2[1]=(pks[1]-y1)*rbin

sdat2=sdat[x1:x2,y1:y2,*]

sdat2=rebin_data(sdat2, rbin)

; rebin the images 
cutoutb=rebin_data(cutout,rbin)             ; ORIGINAL!

xsize=float((size(cutoutb, /dim))[0])
ysize=float((size(cutoutb, /dim))[1])

; Use a Gaussian as the 0-th approximation for the PSF...
cen=max(cutoutb)
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

psf=psf1

; save PSF center locations
xy_psf[0,iloc]=((pks2[0]+shx)/rbin) + x1
xy_psf[1,iloc]=((pks2[1]+shy)/rbin) + y1

; rebin psf back to original bin size
psf=rebin(psf1, xsize1, ysize1)

;use the PSF to get an average flux value and pixel count
xsize=(size(psf, /dim))[0]
ysize=(size(psf, /dim))[1]
psf1d=reform(psf, xsize1*ysize1)
ord1=reverse(sort(psf1d))
signal1=psf1d[ord1]
bkgd=median(signal1)

fwhm_xy=fullwid_halfmax(psf, centroid=cxy, /inner_max)

fwhm_orb[*,orb]=fwhm_xy

if fwhm_xy[0] eq 0 or fwhm_xy[1] eq 0 then goto, endoffitting

minpix1=round(fwhm_xy[0]*fwhm_xy[1])
maxpix1=round(3.*fwhm_xy[0]*3*fwhm_xy[1])

minflux1=total(signal1[0:minpix1])

if maxpix1 gt (xsize*ysize)-1 then maxpix1=(xsize*ysize)-1
maxflux1=total(signal1[0:maxpix1])

sigcum=lonarr(maxpix1)

for i=0, maxpix1-1 do sigcum[i]=total(signal1[0:i])

xx=(where(sigcum ge 0.9*maxflux1))[0]

mostflux=0.9*maxflux1
mostpix=xx

npixel[*,orb]=[minpix1,mostpix,maxpix1]
avgfl[*,orb]=[minflux1,mostflux,maxflux1]

;goto, endoffitting

plotout=outdir+hdname+'_psf.pdf'
plotsym, 0, /fill, 1.5

;check file exists
chk1=file_search(plotout, count=nf)

if nf gt 0 AND orb eq 0 then begin
  spawn, 'rm '+plotout
  nf=0
endif

if nf eq 0 then begin
  
  temp_plot=outdir+'temp.ps'
  
  ps_on, temp_plot, xsize=15, ysize=15
  plot_image, bytscl(psf, 20, 500), title=hdname, charsize=0.7, color=cgcolor('black')
  psfpix1=array_indices(psf, ord1[0:minpix1-1])
  psfpix2=array_indices(psf, ord1[minpix1:mostpix-1])
  psfpix3=array_indices(psf, ord1[mostpix:maxpix1-1])
  oplot, psfpix1[0,*], psfpix1[1,*], psym=8, color=cgcolor('green')
  oplot, psfpix2[0,*], psfpix2[1,*], psym=8, color=cgcolor('purple')
  oplot, psfpix3[0,*], psfpix3[1,*], psym=8, color=cgcolor('orange')
  ps_off 
  
   spawn, 'convert '+temp_plot+' '+plotout
  
endif else begin
  
  temp1=outdir+'temp1.pdf'
  spawn, 'mv '+plotout+' '+temp1
  
  temp_plot=outdir+'temp.ps'
  
  ps_on, temp_plot, xsize=15, ysize=15
  plot_image, bytscl(psf, 20, 500), title=hdname, charsize=0.7, color=cgcolor('black')
  psfpix1=array_indices(psf, ord1[0:minpix1-1])
  psfpix2=array_indices(psf, ord1[minpix1:mostpix-1])
  psfpix3=array_indices(psf, ord1[mostpix:maxpix1-1])
  oplot, psfpix1[0,*], psfpix1[1,*], psym=8, color=cgcolor('green')
  oplot, psfpix2[0,*], psfpix2[1,*], psym=8, color=cgcolor('purple')
  oplot, psfpix3[0,*], psfpix3[1,*], psym=8, color=cgcolor('orange')
  ps_off
  
  temp2=outdir+'temp2.pdf'
  spawn, 'convert '+temp_plot+' '+temp2
  
  ; append files and destroy temp.ps
  spawn, 'PDFconcat -o '+plotout+' '+temp1+' '+temp2
  
  spawn, 'rm '+temp1
  spawn, 'rm '+temp2
  spawn, 'rm '+temp_plot
  
endelse

endoffitting:
end