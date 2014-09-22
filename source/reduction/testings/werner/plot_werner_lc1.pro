pro plot_werner_lc1

indir='~/BRITE/TESTSETS/werner4lc/lc1/'

filein=file_search(indir+'*.txt', count=nf)

for f=0, nf-1 do begin
  
  readcol, filein[f], jd, flux, xcen, ycen, avgtemp, format='d,d,f,f,f)'
  
  mags=-2.5*alog10(flux)
  
  jd1=jd-jd[0]
  
  nimg=n_elements(jd)
  
  wset, 0
  plot, jd1, mags, color=cgcolor('black'), psym=2, xrange=[0,1.], /YNOZERO
  stop
  
  ; now calculate the error 
  jd2=jd1[1:nimg-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  gap=[-1,gap,nimg-1]
  
  err1=0
  
  for orbit=0, ngap do begin
        
    iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
    ni=n_elements(iloc)
    
    flux=flux/median(flux)
    mags=-2.5*alog10(flux)
    ;mags=mags/median(mags)
    
    err1=err1+(stddev(mags[iloc])/sqrt(float(ni)))^2.
    
    ; do x-correlation with the temp
;    print, correlate(avgtemp[iloc], flux[iloc])
;    
;    lag=indgen(ni)
;    
;    result=c_correlate(avgtemp[iloc], flux[iloc], lag)
    
;    wset, 1
;    plot, result, color=cgcolor('black'), psym=8
;    stop
  endfor
  
  err2=sqrt(err1/float(ngap+1))
    
  print, err2*1000.
  
  
  
  
  stop
  
endfor

print, 'End of Program'
end