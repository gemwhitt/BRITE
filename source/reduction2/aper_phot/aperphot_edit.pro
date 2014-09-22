pro aperphot_edit

  ; Use resid to edit the light curves - by doing a "sigma clip"
  ;
  Compile_opt idl2
  
  sat='BA'
  
  field='ORION'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc/'
  
  filein=file_search(indir+'*.txt', count=nf)
  
  if nf eq 0 then stop
  
  for f=0, nf-1 do begin  
    readcol, filein[f], time, flux, xcom, ycom, temperature, vmag, bmag, resid, $
      format='d,d,f,f,f,f,f,f'
      
    time1=time-time[0]
      
    plotsym, 0, /fill, 0.6
    plot, time1, resid, color=cgcolor('black'), psym=8
    oplot, time1, temperature*133., color=cgcolor('purple'), psym=8
    oplot, time1, flux/5.5, color=cgcolor('orange'), psym=8
      
    stop
      
  endfor 
  
  
  
  print, 'End of program'
  print, 'Now do interpixel corrections'
end