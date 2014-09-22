pro aper_phot_intpix

  ; PURPOSE: Apply intra-pixel corrections using flux values and x and y PSF centers
  ;
  Compile_opt idl2
  
  res=10  ; resolution of intrapixel stuff
  
  sat='BA'
  
  field='ORION'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_intpix/'
  
  chkout=file_search(outdir, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+outdir
  
  filein=file_search(indir+'*.txt', count=nf)
  
  for ff=0, nf-1 do begin
  
    readcol, filein[ff], jd, flux, x_com, y_com, avgtemp, vmag, bmag, npix_psf, format='d,d,f,f,f,f,f,i)'
    
    fname=file_basename(filein[ff], '.txt')
    
    npts=n_elements(jd)
    
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
    
    ; save new flux
    fileout=outdir+fname+'.txt'
    openw, lun, fileout, /append, /get_lun
    for ii=0, npts-1 do printf, lun, jd[ii], flux[ii], eflux[ii], x_com[ii], y_com[ii], avgtemp[ii], vmag[ii], bmag[ii], npix_psf[ii], $
      format='(d14.6,x,d14.3,x,d14.3,x,f7.3,x,f7.3,x,f7.3,x,f7.3,x,f7.3,x,i)'
    free_lun, lun
    
    ; convert flux to magnitudes and calculate the error per orbit
    flux=flux/median(flux)
    eflux=eflux/median(eflux)
    
    mags1=(-2.5)*alog10(flux)
    emags1=(-2.5)*alog10(eflux)
    
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
    
    
  endfor
  
  
  
  
  print, 'End of program'
end