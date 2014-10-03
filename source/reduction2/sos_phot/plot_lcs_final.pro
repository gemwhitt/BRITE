pro plot_lcs_final

  Compile_opt idl2
  
  sat='BA'
  
  field='CENTAURUS'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/sos_lc_final/'
    
  filein=file_search(indir+'HD127973.txt', count=nf)
  
  for ff=0, nf-1 do begin
  
    readcol, filein[ff], time, frame, flux, xcom, $
      ycom, temperature, vmag, bmag, resid, $
      format='(d,i,f,f,f,f,f,f,f)'
      
    time1=time-time[0]
    
    nimg=n_elements(time)
    
    wset, 0
    plot, time1, flux, color=cgcolor('black'), psym=2, /YNOZERO;, xrange=[80,100]
    
    ; now calculate the error
    time2=time1[1:nimg-1]
    tdiff=time2-time1
    gap=where(tdiff gt 0.015, ngap)
    gap=[-1,gap,nimg-1]
    
    err1=[]
    
    for orbit=0, ngap do begin
    
      iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
      ni=n_elements(iloc)
      
      flux=flux/median(flux)
      mags=-2.5*alog10(flux)
      ;mags=mags/median(mags)
      
      err1=[err1,(stddev(mags[iloc])/sqrt(float(ni)))^2.]
      
    endfor
    
    xx=where(finite(err1) eq 0, nxx, complement=keep)
    err1=err1[keep]
    nn=ngap+1-nxx
    
    err2=sqrt(total(err1, /nan)/float(nn))
    
    print, vmag[0], err2*1000., ' mmag'
    
    stop
    
  endfor
  
  print, 'End of Program'
end