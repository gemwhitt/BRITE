pro warsaw_snr

  Compile_opt idl2
  
  sat='TOR'
  field='CENTAURUS'
  
  plotsym, 0, /fill, 1.
  
  indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_final/'
  outdir='~/BRITE/reports/presentations/warsaw/'
  
  filein=file_search(indir+'*.txt', count=nf)
  fileout=outdir+sat+'_'+field+'_snr.ps'
  
  if nf eq 0 then stop
  
  sigma=fltarr(nf)
  vmags=fltarr(nf)
  
  for ff=0, nf-1 do begin
    
    readcol, filein[ff], time, frame, flux, xcom, $
      ycom, temperature, vmag, bmag, resid, $
      format='(d,i,f,f,x,x,x,f,f,f,f,f,x)'
      
    vmags[ff]=vmag[0]
      
    time1=time-time[0]
    
    nimg=n_elements(time)
    
    mags=-2.5*alog10(flux/median(flux))
    nmags=mags
    
    time2=time1[1:nimg-1]
    diff=time2-time1
    gap=where(time2 gt 0.015, ngap)
    gap=[-1,gap,nimg-1]
    
    err1=[]
    
    for gp=0, ngap-1 do begin
      iloc=indgen(gap[gp+1]-gap[gp])+(gap[gp]+1)
      
      ni=n_elements(iloc)
      
      if ni le 5 then continue
      
      err1=[err1,(stddev(mags[iloc])/sqrt(float(ni)))^2.]      
      
    endfor
    
    nn=n_elements(err1)
    
    err2=sqrt(total(err1, /nan)/float(nn))
    
    sigma[ff]=err2
    
  endfor
  
  plotsym, 0, /fill, 1.2
  plot, vmags, sigma, color=cgcolor('black'), psym=8, xrange=[2,5]
  
  stop
  
  

end

