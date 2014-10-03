pro sigclip_sosphott, sat, field, target

  ; Do a final sigma clip on the corrected light curves - per orbit rather than by phase folding
  ;
  Compile_opt idl2
  
    sat='BA'
  ;
  field='CENTAURUS'
  
  target='HD127973'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/sos_lc_sav/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/sos_lc_final/'
  chk=file_search(outdir, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+outdir
  
  filein=file_search(indir+target+'*.sav', count=nf)
  
  plotsym, 0, /fill, 0.6
  
  for ff=0, nf-1 do begin
  
    restore, filein[ff] ;time, frame, flux, xcom, ycom, temperature, $
    ;vmag, bmag, resid, psf_aper, data
    
    fname=file_basename(filein[ff], '.sav')
      
    caldat, time, mon, day, yr
    xx=where(yr lt 2013, nxx, complement=keep)
    
    time=time[keep]
    frame=frame[keep]
    flux=flux[keep]
    xcom=xcom[keep]
    ycom=ycom[keep]
    temperature=temperature[keep]
    resid=resid[keep]
    
    nimg=n_elements(time)
    
    time1=time-time[0]
    time2=time1[1:nimg-1]
    diff=time2-time1
    gap=where(diff gt 0.015, ngap)
    xx=indgen(gap[0])
    gap=[-1,gap,nimg-1]
    cadence=median(diff[xx])  ; in days
    totdur=time1[nimg-1]-time1[0]
    
    ;  wset, 0
    plot, time1, flux, color=cgcolor('black'), psym=2, /YNOZERO;, xrange=[40,60]
    ;  stop
    ;
    flux2=flux*0
    ; first get the median of the points in each orbit
    med_orb=fltarr(ngap+1)
    med_time=dblarr(ngap+1)
    for orb=0, ngap do begin
      iloc=indgen(gap[orb+1]-gap[orb])+(gap[orb]+1)
      
      med_orb[orb]=median(flux[iloc])
      
      med_time[orb]=median(time1[iloc])
      
      flux2[iloc]=flux[iloc]-med_orb[orb]
      
    endfor
    
;    oplot, med_time, med_orb, color=cgcolor('green'), psym=8
;    stop
;    plot, time1, flux2, color=cgcolor('black'), psym=2, /YNOZERO;, xrange=[40,60]
;    stop
    med_res=median(flux2)
    sig_res=robust_sigma(flux2)
    
    ;  plotsym, 0, /fill, 0.8
    out=where(flux2 le med_res-(5*sig_res) OR flux2 ge med_res+(5*sig_res), nout, complement=keep)
;    if nout gt 0 then oplot, time1[out], flux2[out], psym=8, color=cgcolor('red')
;    stop
    
    
    ;;
    time=time[keep]
    flux=flux[keep]
    frame=frame[keep]
    xcom=xcom[keep]
    ycom=ycom[keep]
    temperature=temperature[keep]
    vmag=vmag[keep]
    bmag=bmag[keep]
    resid=resid[keep]
    
    
    ; plot, time1, flux, color=cgcolor('black'), psym=8, /ynozero
    
    ; write out new file with correction applied
    fileout=outdir+fname+'.txt'
    openw, lun, fileout, /get_lun
    for ii=0, n_elements(keep)-1 do  printf, lun, time[ii], frame[ii], flux[ii], xcom[ii], $
      ycom[ii], temperature[ii], vmag[ii], bmag[ii], resid[ii], $
      format='(d14.6,x,i,x,f,x,f,x,f,x,f7.3,x,f7.3,x,f7.3,x,d9.2)'
    free_lun, lun
    
  endfor
  
  print, 'End of Progam'
  print, 'Now plot LCs and do stats'
end
