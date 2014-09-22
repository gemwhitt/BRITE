pro aper_phot

; Extract photometry from images using aper_psf and xy_com with data1

; Save to a text file with all variables

Compile_opt idl2

nbin=10

sat='BA'

field='ORION'

indir='~/BRITE/'+sat+'/'+field+'/data/p5/'

outdir='~/BRITE/'+sat+'/'+field+'/data/aper_lc/'

chkout=file_search(outdir, count=nchk)
if nchk eq 0 then spawn, 'mkdir -p '+outdir

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

for f=0, 2 do begin ;nf-1 do begin
  
  fname=file_basename(filesin[f], '_p3_p4.sav')
  
  restore, filesin[f] ;roi_name, jd, data1, roi_loc, ccd_temp, $
  ;vmag, bmag, flag, psf_loc, npix_psf, modelpsf, xy_com, psf_aper
  
  nimg=n_elements(jd)
  s=size(data1, /dim)
  
  jd1=jd-jd[0]
  
  good=where(flag eq 2, ngood)
  
  for im=0, nimg-1 do begin
    
    if flag[im] ne 2 then continue
      
    ; define new arrays to save the binned and shifted PSFs for this orbit
    binarray=fltarr(s[0]*nbin,s[1]*nbin)
    shiftarray=fltarr(s[0]*nbin,s[1]*nbin)
  
    im0=data1[*,*,im]
    
    ; calulcate the background using psf_loc (to avoid the PSF pixels) 
    backimg=float(im0)
    backimg[psf_loc[0,im]:psf_loc[1,im],psf_loc[2,im]:psf_loc[3,im]]=!Values.F_NAN
    
    backgd=backimg[where(finite(backimg) eq 1, nbkgd)]
        
    ; do histogram of background pixels and pick the value which is most represented but not zero
    result=histogram(backgd, min=1, binsize=1, locations=loc)
    
    sort1=reverse(sort(result))
    
    ; choose the most populated value
    bkgd=loc[sort1[0]]
    bkgd_err=sqrt((total((backgd-mean(backgd))^2))/float(nbkgd))
    
    ; remove the background value from im0
    im0=im0-bkgd
    
    ; bin the image
    im1=congrid(im0, nbin*s[0], nbin*s[1], /interp, /center)
    
    ; now modify xy_com to obtain the center of PSF in this binned image
    xcen=fix(round(xy_com[0,im]*nbin))
    ycen=fix(round(xy_com[1,im]*nbin))
    
    ; get the center of the binned image
    xcen2=s[0]*nbin/2.
    ycen2=s[1]*nbin/2.
    
    ; calculate the shifts
    xdiff=xcen2-xcen
    ydiff=ycen2-ycen
    
    ; now shift aperture to the center of the PSF in the image
    aper=psf_aper*0
    
    ; find the locations of ALL PSF pixels gt 0
    psfpix=where(psf_aper gt 0, npsf)
    xpix=(array_indices(psf_aper,psfpix))[0,*]
    ypix=(array_indices(psf_aper,psfpix))[1,*]
    
    aper[xpix-xdiff,ypix-ydiff]=psf_aper[xpix,ypix]
   
    modpix=where(aper gt 0, npsf)
    xpix=(array_indices(aper,modpix))[0,*]
    ypix=(array_indices(aper,modpix))[1,*]
    
    
;        plot_image, bytscl(im1, 0, 500)
;        
;        plotsym, 0, /fill, 0.2
;        oplot, xpix, ypix, color=cgcolor('purple'), psym=8
        
        ; do a goodness of fit by calculating the residual of model to image
;        res1=im1
;        res1[xpix,ypix]=0
;        
;        res2=(im1-aper)/nbin^2
;        q1=max(res2)-min(res2)
;        qual1=[qual1,q1]
;        
;        if q1 gt 10 then begin
;          plot_image, bytscl(im1, 0, 500), title=q1, charsize=0.8, color=cgcolor('black')
;          plotsym, 0, /fill, 0.2
;          
;          oplot, xpix, ypix, color=cgcolor('purple'), psym=8
;          stop
;        endif

; calculate the residual in the image after subtracting the flux within the model aperture
    im2=im1
    im2[xpix,ypix]=0
    resid=total(im2)/nbin^2.
    
   ; if resid ge 7000 then stop
    
    if jd1[im] ge 16. AND resid ge 7000 then stop
          
    flux=total(im1[xpix,ypix])/nbin^2.
    time=jd[im]
    xcom=xy_com[0,im]
    ycom=xy_com[1,im]
    temperature=average(ccd_temp[*,im])
    
    ; save results
    fileout=outdir+fname+'.txt'
    
    openw, lun, fileout, /get_lun, /append
    printf, lun, time, flux, xcom, ycom, temperature, vmag, bmag, resid, $
      format='(d14.6,x,d14.3,x,f7.4,x,f7.4,x,f7.3,x,f7.3,x,f7.3,x,f7.2)'
    free_lun, lun
    
    npts=n_elements(time)
    
    endfor  ; end loop over image
    

;    plotsym, 0, /fill, 0.2
;    plot, time, flux, color=cgcolor('black'), psym=8
;    
    

endfor  ; end loop over file

print, 'End of Program'
print, 'Now do sigma clipping with resid, then interpix, then sigma clip flux using smooth'

end

