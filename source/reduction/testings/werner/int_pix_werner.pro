pro int_pix_werner

  ; PURPOSE: Apply intra-pixel corrections using flux values and x and y PSF centers
  ;
  Compile_opt idl2
  
  res=10  ; resolution of intrapixel stuff
  
  indir='~/BRITE/TESTSETS/werner4lc/lc1/'
  
  outdir='~/BRITE/TESTSETS/werner4lc/lc1b/'
  
  filein=file_search(indir+'*.txt', count=nf)
  
  for ff=0, nf-1 do begin
  
    readcol, filein[ff], jd, flux, x_com, y_com, avgtemp, format='d,d,f,f,f)'
        
    ; put pixel coords into actual pixel size
    xpsf=x_com
    ypsf=y_com
    
    npts=n_elements(flux)
    
    eflux=dblarr(npts)
    
    ; get fraction of pixel coord
    xpsf=(xpsf mod 1)
    ypsf=(ypsf mod 1)
    
    xpsf=xpsf*res
    ypsf=ypsf*res
    
    mask1=fltarr(res,res)
    
    for x=0, res-1 do begin
      for y=0, res-1 do begin
      
        mask1[x,y]=total(flux[where(fix(xpsf) eq x AND fix(ypsf) eq y, npix)])/float(npix)
        
        if npix eq 0 then mask1[x,y]=mean(flux)
        
      endfor
    endfor
    
    mask1=mask1 - median(mask1)    ; original
    ;mask1=mask1/median(mask1)    ; ????
    
    for n=0, npts-1 do begin
    
      eflux[n]=flux[n]-mask1[xpsf[n],ypsf[n]]  ;- original
      ;eflux[n]=flux[n]*mask1[xpsf[n],ypsf[n]]
      
    endfor
    
    ; convert flux to magnitudes and calculate the error per orbit
    flux=flux/median(flux)
    eflux=eflux/median(eflux)
    
    mags1=(-2.5)*alog10(flux)
    emags1=(-2.5)*alog10(eflux)
    
    wset, 0
    plot, jd, (-1)*emags1, color=cgcolor('black'), psym=2, /YNOZERO
    stop
    mags=mags1
    emags=emags1
    
    jd1=jd-jd[0]
    jd2=jd1[1:npts-1]
    jdiff=jd2-jd1
    gap=where(jdiff gt 0.015, ngap)
    gap=[-1,gap,npts-1]
    
    fsig=fltarr(ngap+1)
    
    efsig=fltarr(ngap+1)
    
    for orbit=0, ngap do begin
    
      iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
      
      ni=n_elements(iloc)
      
      fsig[orbit]=(stddev(mags[iloc])/sqrt(ni))^2.
      
      efsig[orbit]=(stddev(emags[iloc])/sqrt(ni))^2.
      
    endfor
    
    flxrms=(sqrt(total(fsig)/float(ngap+1)))*1000.
    
    eflxrms=(sqrt(total(efsig)/float(ngap+1)))*1000.
    
    print, flxrms, eflxrms
    stop
    ; print out results
    openw, lun, outdir+'int_pix_stats.txt', /get_lun, /append
    printf, lun, name, flxrms, eflxrms, format='(a,x,f,x,f)'
    free_lun, lun
    
  endfor
  
  
  
  
  print, 'End of program'
end