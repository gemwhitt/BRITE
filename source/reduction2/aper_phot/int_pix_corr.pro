pro int_pix_corr, sat, field, target

  ; PURPOSE: Apply intra-pixel corrections using flux values and x and y PSF centers
  ;
  Compile_opt idl2
  
;  sat='BA'
if n_elements(sat) eq 0 then sat='TOR'
if n_elements(field) eq 0 then field='CENTAURUS'
if n_elements(target) eq 0 then target='HD127973'
  
;  field='CENTAURUS'
  
  res=10  ; resolution of intrapixel stuff
  
  indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_best/'
  ;indir='~/BRITE/'+sat+'/'+field+'/data/sos_lc_sav/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_intpix/'
  chk=file_search(outdir, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+outdir
  
  filein=file_search(indir+target+'*', count=nf)
  
  if nf eq 0 then stop
;  stop
  for ff=0, nf-1 do begin
  
    readcol, filein[ff], time, frame, flux, mags, bkgd, bkgd_err, xcom, $
      ycom, temperature, vmag, bmag, resid, medimg, $
      format='(d,i,f,f,f,f,f,f,f,f,f,f,i)'
 ; restore, filein[ff]
;      
      caldat, time, mon, day, yr
    xx=where(yr lt 2013, nxx, complement=keep)
    
    time=time[keep]
    frame=frame[keep]
    flux=flux[keep]
    mags=mags[keep]
    bkgd=bkgd[keep]
    bkgd_err=bkgd_err[keep]
    xcom=xcom[keep]
    ycom=ycom[keep]
    temperature=temperature[keep]
    resid=resid[keep]
    medimg=medimg[keep]
    
    nimg=n_elements(time)
    
    ; put pixel coords into actual pixel size
    xpsf=xcom
    ypsf=ycom
    
    npts=n_elements(flux)
    
    eflux1=dblarr(npts)
    eflux2=dblarr(npts)
    
    ; get fraction of pixel coord
    xpsf1=(xpsf mod 1)
    ypsf1=(ypsf mod 1)
    
    xpsf2=xpsf1*res
    ypsf2=ypsf1*res
    
    mask0=fltarr(res,res)
    
    for x=0, res-1 do begin
      for y=0, res-1 do begin
      
        mask0[x,y]=total(flux[where(fix(xpsf2) eq x AND fix(ypsf2) eq y, npix)])/float(npix)
        
        if npix eq 0 then mask1[x,y]=mean(flux)
        
      endfor
    endfor
    
    mask1=mask0 - median(mask0)    ; original
    mask2=mask0/median(mask0)    ; ????
    
;    wset, 0
    plot, xcom, flux, color=cgcolor('black'), psym=2, /ynozero
;    
;    plot, time, flux, color=cgcolor('black'), psym=2
;    stop
;    wset, 1
;    plot_image, mask1
;    
;    stop
    
    for n=0, npts-1 do begin
    
      eflux1[n]=flux[n]-mask1[xpsf2[n],ypsf2[n]]  ;- original
      eflux2[n]=flux[n]*mask2[xpsf2[n],ypsf2[n]]
      
    endfor  ; end loop over points
    
    plot, xcom, flux, color=cgcolor('black'), psym=2, /ynozero, xrange=[10,12]
    stop
;    
;;    wset, 0
    plot, xcom, eflux1, color=cgcolor('black'), psym=2, /ynozero, xrange=[10,12]
    stop
;  ;  wset, 1
    plot, xcom, eflux2, color=cgcolor('black'), psym=2, /ynozero, xrange=[10,14]
    stop
;;stop
flux=eflux1
mags=-2.5*alog10(flux)

; write out new file with correction applied
fileout=outdir+file_basename(filein[ff])
openw, lun, fileout, /get_lun
for ii=0, npts-1 do  printf, lun, time[ii], frame[ii], flux[ii], mags[ii], bkgd[ii], bkgd_err[ii], xcom[ii], $
    ycom[ii], temperature[ii], vmag[ii], bmag[ii], resid[ii], medimg[ii], $
    format='(d14.6,x,i,x,f,x,f,x,f,x,f7.4,x,f,x,f,x,f7.3,x,f7.3,x,f7.3,x,f,x,i)'
free_lun, lun

    
endfor  ; end loop over file
  
print, 'End of program'
print, 'Plot Light curves and run analysis'
end