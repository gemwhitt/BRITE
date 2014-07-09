pro get_psf_info, dat,gap,ngap,nfrm,outdir,hdname, xy_psf,npix,e_npix,fl,e_fl,fwhm

; 20 March 2014
; called by p3_p4_gem6
; use slavek programs to get PSF centers
; get average flux per orbit and averahe number of pixels per orbit
; 
gap=[0,gap,nfrm-1]

npixel=lonarr(3,ngap+1)
avgfl=lonarr(3,ngap+1)
fwhm_orb=fltarr(2,ngap+1)

for g=0, ngap do begin
  
  if g eq 0 then iloc=indgen(gap[g+1]+1) else iloc=indgen(gap[g+1]-gap[g])+gap[g]+1
  
  sdat=dat[*,*,iloc]  ; subset of dat - for one orbit
  
  ; get the PSF centers
  sr_fitting, sdat,iloc,g,outdir,hdname, xy_psf,npixel,avgfl,duds,fwhm_orb
  
endfor

; remove bad orbits
bad=where(avgfl[0,*] lt 500, nbad, complement=good)

if nbad gt 0 then begin
  npixel=npixel[*,good]
  avgfl=avgfl[*,good]
  fwhm_orb=fwhm_orb[*,good]
endif

; get average values for the whole observation sequence
for k=0, 2 do npix[k]=median(npixel[k,*])
for k=0, 2 do e_npix[k]=robust_sigma(npixel[k,*]/npix[k])

for k=0, 2 do fl[k]=median(avgfl[k,*])
for k=0, 2 do e_fl[k]=robust_sigma(avgfl[k,*]/fl[k])

fwhm[0]=median(fwhm_orb[0,*])
fwhm[1]=median(fwhm_orb[1,*])

end
  
