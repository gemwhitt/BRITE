pro aperphot_edit

  ; Use resid to edit the light curves - by doing a "sigma clip"
  ;
  Compile_opt idl2
  
  sat='TOR'
  ;sat='BA'
  
  field='CENTAURUS'
  ;field='ORION'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc/'
  
  filein=file_search(indir+'*.txt', count=nf)
  
  if nf eq 0 then stop
  
  for f=0, nf-1 do begin  
    readcol, filein[f], time, flux, bkgd, bkgd_err, xcom, ycom, temperature, vmag, bmag, resid, $
      format='d,d,f,f,f,f,f,f,f,d'
      
    time1=time-time[0]
    
    pcflux=flux/(flux+resid)*100.
    
    ;calculate median, min and max
    med=median(pcflux)
    minv=med-(3*robust_sigma(pcflux))
    maxv=med+(3*robust_sigma(pcflux))
    
    wset, 0  
    plotsym, 0, /fill, 0.9
    plot, time1, pcflux, color=cgcolor('black'), psym=8, /ynozero;, xrange=[3,3.04] 
    oplot, [0,40], [minv,minv], color=cgcolor('purple'), thick=2
    oplot, [0,40], [maxv,maxv], color=cgcolor('purple'), thick=2
    
    wset, 1
    plot, time1, temperature, color=cgcolor('black'), psym=8, /ynozero  ;xrange=[3,3.04], 
    
    wset, 2
    plot, time1, flux, color=cgcolor('black'), psym=8, /ynozero ;xrange=[3,3.04], /ynozero
    
    wset, 3
    plot, time1, bkgd, color=cgcolor('black'), psym=8, /ynozero ;xrange=[3,3.04], /ynozero
    
      
    stop
      
  endfor 
  
  
  
  print, 'End of program'
  print, 'Now do interpixel corrections'
end