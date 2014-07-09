pro int_pix1

; PURPOSE: Apply intra-pixel corrections using flux values and x and y PSF centers
; 
Compile_opt idl2

sat='BA'
field='ORION'

nb='1'

res=10  ; resolution of intrapixel stuff

indir='~/BRITE/data/'+sat+'/sos_phott'+nb+'/'+field+'/lc_txt/'

outdir='~/BRITE/data/'+sat+'/sos_phott'+nb+'/'+field+'/stats/'

filesin=file_search(indir+'*.txt', count=nf)

for ff=0, nf-1 do begin
  
  readcol, filesin[ff], jd, flux, x_com, y_com, format='d,d,x,x,x,f,f,x,x,x,x', comment='#'
  
  readcol, filesin[ff], name, format='x,a', numline=1
  
  ; put pixel coords into actual pixel size
  xpsf=x_com/float(nb)
  ypsf=y_com/float(nb)
  
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
 mags=(-2.5)*alog10(flux)
 emags=(-2.5)*alog10(eflux)
 
 mags=mags/median(mags)
 emags=emags/median(emags)
 
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
      
      fsig[orbit]=(stddev(mags[iloc]))^2.
      
      efsig[orbit]=(stddev(emags[iloc]))^2.
            
    endfor
    
    flxrms=(sqrt(total(fsig)/float(ngap+1)))*1000.
    
    eflxrms=(sqrt(total(efsig)/float(ngap+1)))*1000.
    
    ; print out results
    openw, lun, outdir+'int_pix_stats.txt', /get_lun, /append
    printf, lun, name, flxrms, eflxrms, format='(a,x,f,x,f)'
    free_lun, lun

endfor




print, 'End of program'
end