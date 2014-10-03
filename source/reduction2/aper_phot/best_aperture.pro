pro best_aperture, sat, field, target

; determine the best aperture from the flux and residuals - use sigma clipping first
; 
Compile_opt idl2

;sat='BA'

;field='CENTAURUS'

indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_sav/'

outdir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_best/'
chkout=file_search(outdir, count=nchk)
if nchk eq 0 then spawn, 'mkdir -p '+outdir

statsdir='~/BRITE/'+sat+'/'+field+'/reduction/best_aper/'
chk=file_search(statsdir, count=nchk)
if nchk eq 0 then spawn, 'mkdir -p '+statsdir

statsfile=statsdir+'best_ap.txt'
chk=file_search(statsfile, count=nchk)
if nchk ne 0 then spawn, 'rm '+statsfile

filesin=file_search(indir+target+'*.sav', count=nf)

if nf eq 0 then stop

plotsym, 0, /fill, 0.8

for f=0, nf-1 do begin
  
  restore, filesin[f] ;time, frame, flux, bkgd, bkgd_err, xcom, ycom, temperature, $
  ;vmag, bmag, resid, medimg, psf_aper, data
  
  fname=file_basename(filesin[f],'.sav')
  
  time1=time-time[0]
  npts=n_elements(time1)
  
  nap=(size(flux, /dim))[1] ; number of light curves to analyse
   
  time2=time1[1:npts-1]
  diff=time2-time1
  gap=where(diff gt 0.015, ngap)
  gap=[-1,gap,npts-1]
  pporb=gap[1]-gap[0]
  
  for ap=0, nap-1 do begin
    
    for orb=0, ngap do begin
      
      iloc=indgen(gap[orb+1]-gap[orb])+(gap[orb]+1) ; index of points in this orbit
      
      ni=n_elements(iloc)
      
      if ni lt 5 then begin
        ; make points Nan values
        flux[iloc,ap]=!values.F_nan
        resid[iloc,ap]=!values.F_nan
        goto, next_orbit
      endif
      
      ; smooth over the orbit
      flux1=flux[iloc,ap]
      sm1=smooth(flux1, ni/1.5, /edge_mirror)
      
;      plot, flux1, color=cgcolor('black'), psym=8, /ynozero
;      oplot, sm1, color=cgcolor('orange')
;      stop
      flux2=flux1-sm1+median(flux1)
      
      sigmaf=robust_sigma(flux2)
      
      clip=where(flux2 ge median(flux2)+(3*sigmaf) OR flux2 le median(flux2)-(3*sigmaf), nclip, complement=keep)
      
;      if nclip gt 0 then oplot, clip, flux1[clip], color=cgcolor('red'), psym=8 else print, 'none'
      
      ; make clipped points Nan values
      flux[iloc[clip],ap]=!values.F_nan
      resid[iloc[clip],ap]=!values.F_nan
      
    next_orbit:
    endfor  ; end loop over orbit
  
  endfor  ; end loop over all apertures
    
  ; now compute light curve scatter (RMS) and pcflux using residuals
  sig_lc=fltarr(nap)
  pcflux=fltarr(nap)

  for ap=0, nap-1 do begin
    
    sig2=[]
  
    for orb=0, ngap do begin
    
      iloc=indgen(gap[orb+1]-gap[orb])+(gap[orb]+1) ; index of points in this orbit
      
      flux1=flux[iloc,ap]
      flux1=flux1[where(finite(flux1) eq 1, ngood)]
           
      if ngood lt 3 then continue
      
      flux2=flux1/median(flux1)
      ni=n_elements(flux2)
      
      sig2=[sig2,(stddev(flux2)/sqrt(ni))^2.]
      
    endfor  ; end loop over orbit
    
    sig_lc[ap]=(sqrt(total(sig2)/float(ngap+1)))
    flux0=flux[*,ap]
    flux0=flux0[where(finite(flux0) eq 1)]
    res0=resid[*,ap]
    res0=res0[where(finite(flux0) eq 1)]
    pcflux[ap]=median(flux0/(flux0+res0)*100)
    
  endfor  ; end loop over all apertures
  
  ; print results to file
  openw, lun, statsfile, /get_lun, /append
  printf, lun, fname, sig_lc, format='(a,x,f,x,f,x,f,x,f)'
  printf, lun, fname, pcflux, format='(a,x,f,x,f,x,f,x,f)'
  printf, lun, '', format='(a)'
  free_lun, lun
  
  ; choose the best light curve by going with the lowest sigma
  best=where(sig_lc eq min(sig_lc), nbest)  ; best indicates best ap
  
  if nbest gt 1 then stop ; no choose highest pcflux
  
  ; resave new flux
  newflux=reform(flux[*,best])
  
  keep=where(finite(newflux) eq 1, nkeep)
  
  ; modify other arrays accordingly
  time=time[keep]
  flux=newflux[keep]
  frame=frame[keep]
  bkgd=bkgd[keep]
  bkgd_err=bkgd_err[keep]
  xcom=xcom[keep]
  ycom=ycom[keep]
  temperature=temperature[keep]
  vmag=vmag[keep]
  bmag=bmag[keep]
  resid=reform(resid[*,best])
  resid=resid[keep]
  medimg=medimg[keep]
  
  mags=-2.5*alog10(flux)
  
  ; write out text file with new LC  
  fileout=outdir+fname+'.txt'
  openw, lun, fileout, /get_lun
  for ii=0, nkeep-1 do printf, lun, time[ii], frame[ii], flux[ii], mags[ii], bkgd[ii], bkgd_err[ii], xcom[ii], $
    ycom[ii], temperature[ii], vmag[ii], bmag[ii], resid[ii], medimg[ii], $
    format='(d14.6,x,i,x,f,x,f,x,f,x,f7.4,x,f,x,f,x,f7.3,x,f7.3,x,f7.3,x,f,x,i)'
  free_lun, lun
   
endfor  ;end loop over all files
 
print, 'End of program'
print, 'Do int pix correction on best light curve'
end