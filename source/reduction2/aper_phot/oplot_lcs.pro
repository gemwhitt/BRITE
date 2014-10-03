pro oplot_lcs

  Compile_opt idl2
  
  sat=['BA','UB','TOR']
  
  field='CENTAURUS'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_final/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/plots/aper_lc_final/'
  
  target='HD127973' ;HD127973
  
  filein=file_search(indir+target+'*.txt', count=nf)
  
  plotsym, 0, /fill, 0.5
  cols1=[cgcolor('red'), cgcolor('blue'), cgcolor('green'),cgcolor('purple')]
  
  for ff=0, nf-1 do begin
    
    readcol, filein[ff], time, frame, flux, mags, bkgd, bkgd_err, xcom, $
      ycom, temperature, vmag, bmag, resid, medimg, $
      format='(d,i,f,f,f,f,f,f,f,f,f,f,i)'
      
      time1=time-2454000D
      
    if ff eq 0 then plot, time1, flux, /ynozero, xrange=[2740,2870], /nodata, color=cgcolor('black'), xstyle=1

    oplot, time1, flux, color=cols1[ff], psym=8
    print, filein[ff]
      
      stop
      
  endfor
  
  stop
  
  
  
  end