pro aper_phot

; Extract photometry from images using aper_psf and xy_com with data1

; Save to a text file with all variables

; 2 apertures in psf_aper

Compile_opt idl2

nbin=10

sat='TOR'

field='CENTAURUS'

indir='~/BRITE/'+sat+'/'+field+'/data/p5/'

outdir='~/BRITE/'+sat+'/'+field+'/data/aper_lc/'

chkout=file_search(outdir, count=nchk)
if nchk eq 0 then spawn, 'mkdir -p '+outdir

targets=['HD127973','HD129056','HD135379','HD136415']

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

for f=0, nf-1 do begin
  
  fname=file_basename(filesin[f], '.sav')
  
  restore, filesin[f] ;roi_name, jd, data1, roi_loc, ccd_temp, $
  ;vmag, bmag, flag, psf_loc, npix_psf, modelpsf, xy_com, psf_aper
  
  nimg=n_elements(jd)
  s=size(data1, /dim)
  
  jd1=jd-jd[0]
  
  good=where(flag eq 2, ngood)
  
  ; if there is already an output file for this target, then delete it
  chk=file_search(outdir+fname+'.txt', count=nchk)
  if nchk gt 0 then spawn, 'rm '+outdir+fname+'.txt'
  
  
  for im=0, nimg-1 do begin
    
    flux=[]
    resid=[]
    
    if flag[im] ne 2 then continue
      
    ; define new arrays to save the binned and shifted PSFs for this orbit
    binarray=fltarr(s[0]*nbin,s[1]*nbin)
    shiftarray=fltarr(s[0]*nbin,s[1]*nbin)
  
    im0=data1[*,*,im]
    
    ; calulcate the background using psf_loc (to avoid the PSF pixels) 
    backimg=float(im0)
    backimg[psf_loc[0,im]:psf_loc[1,im],psf_loc[2,im]:psf_loc[3,im]]=!Values.F_NAN
    
    backgd=backimg[where(finite(backimg) eq 1, nbkgd)]
    
    fr = fractile(backgd, [.25,.75])
    
    xx=where(backgd ge fr[0] AND backgd le fr[1], nxx)
    
 ;   print, total(backgd[xx])/float(nxx)
    
    bkgd=total(backgd[xx])/float(nxx)
    
    bkgd_err=sqrt((total((backgd[xx]-bkgd)^2))/float(nxx))
    
    
    
    ; HISTOGRAM METHOD.....    
    ; do histogram of background pixels and pick the value which is most represented but not zero
;    result=histogram(backgd, min=1, binsize=1, locations=loc)
;    
;    sort1=reverse(sort(result))
;    
;    ; choose the most populated value
;    bkgd=loc[sort1[0]]
;    bkgd_err=sqrt((total((backgd-mean(backgd))^2))/float(nbkgd))

; AVERAGE method

    
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
    sap=size(psf_aper, /dim)
    nap=sap[2]
    
    for ap=0, nap-1 do begin
      modelapp=reform(psf_aper[*,*,ap])
      
      aper=modelapp*0
      
      ; find the locations of ALL PSF pixels gt 0
      psfpix=where(modelapp gt 0, npsf)
      xpix=(array_indices(modelapp,psfpix))[0,*]
      ypix=(array_indices(modelapp,psfpix))[1,*]
      
      aper[xpix-xdiff,ypix-ydiff]=modelapp[xpix,ypix]
      
      modpix=where(aper gt 0, npsf)
      xpix=(array_indices(aper,modpix))[0,*]
      ypix=(array_indices(aper,modpix))[1,*]
      
      flux=[flux,total(im1[xpix,ypix])/nbin^2.]
      
      ; calculate the residual in the image after subtracting the flux within the model aperture
      im2=im1
      im2[xpix,ypix]=0
      resid=[resid,total(im2)/nbin^2.]
      
    endfor
    
    time=jd[im]
    xcom=xy_com[0,im]
    ycom=xy_com[1,im]
    temperature=average(ccd_temp[*,im])
    
    ; save results
    fileout=outdir+fname+'.txt'
    
    openw, lun, fileout, /get_lun, /append
    printf, lun, time, im, flux[0], flux[1], bkgd, bkgd_err, xcom, ycom, temperature, vmag, bmag, resid[0], resid[1], $
      medimg0[im], average(medcol[*,im]), $
      format='(d14.6,x,i6,x,d14.3,x,d14.3,x,f7.2,x,f7.2,x,f7.4,x,f7.4,x,f7.3,x,f7.3,x,f7.3,x,d14.3,x,d14.3,x,f7.1,x,f7.2)'
    free_lun, lun
    
    undefine, flux, resid
    
    npts=n_elements(time)
    
    endfor  ; end loop over image
    

;    plotsym, 0, /fill, 0.2
;    plot, time, flux, color=cgcolor('black'), psym=8
;    
    

endfor  ; end loop over file

print, 'End of Program'
print, 'Now do sigma clipping with resid, then interpix, then sigma clip flux using smooth'

end

