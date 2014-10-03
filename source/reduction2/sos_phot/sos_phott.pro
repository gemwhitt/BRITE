pro sos_phott, sat, field, target

  ; Get light curves, plots and figures using sum-of-squares method - based on a SNR in EACH image
  
  Compile_opt idl2

  rbin=10.
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p5/'
  
  outdir='~/BRITE/'+sat+'/'+field+'/data/sos_lc_sav/'
  
  chkout=file_search(outdir, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+outdir
  
  filesin=file_search(indir+target+'*.sav', count=nf)
  
  if nf eq 0 then stop
  
  for ff=0, nf-1 do begin
  
    fname=file_basename(filesin[ff], '.sav')
      
    restore, filesin[ff]  ;roi_name, jd, data1, roi_loc, ccd_temp, $
   ; vmag, bmag, flag, psf_loc, npix_psf, modelpsf, xy_com, psf_aper, medimg0, medcol
    
    jd0=jd-jd[0]
    hdname=roi_name[0]
    frame=[]
    
    nimg=n_elements(jd)
    
    ; set up arrays for collecting results
    npix=lonarr(nimg)         ; number of pixels (binned or otherwise) in PSF
    flux=dblarr(nimg)         ; flux array
    resid=fltarr(nimg)
    
    for im=0, nimg-1 do begin
    
      if flag[im] ne 2 then continue
      
      im0=data1[*,*,im]
      
      ; get dimensions of image
      xdim=(size(im0, /dim))[0]
      ydim=(size(im0, /dim))[1]
            
      ; Use psf_loc [*,im] to make a cutout around the PSF
      x1=(psf_loc[0,im]-2 > 0)
      x2=(psf_loc[1,im]+2 < xdim-1)
      y1=(psf_loc[2,im]-2 > 0)
      y2=(psf_loc[3,im]+2 < ydim-1)
      
      im1=im0[x1:x2,y1:y2]  ; smaller cutout around psf
      s=size(im1, /dim)
      
      ; bin the data according to rbin
      ;im1=rebin_data(im1,rbin)
      
      im1=congrid(im1, rbin*s[0], rbin*s[1], /interp, /center)
      
      ; get a value for the total number of "pixels" in imb
      npixel=n_elements(im1)
      
      ;reorder pixels by data number
      sb=size(im1, /dim)            ; size of binned array
      
      im1d=reform(im1, sb[0]*sb[1])
      sort1=reverse(sort(im1d))
      im1d=im1d[sort1]
      
      count=2 ; start on count > 0
      sig1=total(im1d[0:1])*3.25
      snr=0
      
      ; plot, [count], [snr], color=cgcolor('black'), psym=8, xrange=[0,80], yrange=[0,2000]
      
      repeat begin
      
        count=count+1
        snr1=snr
        
        sig1=total(im1d[0:count])*3.25
        
        noise1=sqrt(im1d[count] + count + (15.^2))
        
        snr=sig1/noise1
        
        ; oplot, [count], [snr], color=cgcolor('black'), psym=8
        
      endrep until snr lt snr1
      ; stop
      xx=count-1          ; initial aperture
      psfp=sort1[0:xx]
      xflx=(array_indices(im1, psfp))[0,*]
      yflx=(array_indices(im1, psfp))[1,*]
      
      npix[im]=xx
      flux[im]=total(im1d[0:xx])/float(rbin^2)
      resid[im]=total(im0)-flux[im]
      
;      plot_image, bytscl(im1, 20, 500)
;      plotsym, 0, /fill, 0.9
;      oplot, xflx, yflx, color=cgcolor('purple'), psym=8
      
      frame=[frame,im]      
      
    endfor  ; end loop over images
    
    good=where(npix ne 0, ngood)
        
    ; keep only good points
    jd=jd[good]
    data1=data1[*,*,good]
    flag=flag[good]
    flux=flux[good]
    npix=npix[good]
    ccd_temp=ccd_temp[*,good]
    xy_com=xy_com[*,good]
    frame=frame[good]
    resid=resid[good]
    n=ngood

    
    ; do a sigma clip per orbit
    jd0=jd-jd[0]
    jd1=jd0[1:n-1]
    jdiff=jd1-jd0
    gap=where(jdiff gt 0.015, ngap)
    gap=[-1,gap,n-1]
    
    good_bad=intarr(ngood)
    
    for orbit=0, ngap do begin
    
      iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
      
      if n_elements(iloc) lt 5 then good_bad[iloc]=1
      if n_elements(iloc) lt 5 then CONTINUE
      
      ; nsub pix per orbit
      orbpx=npix[iloc]
      oflx=flux[iloc]
      med1=median(orbpx)
      sig1=robust_sigma(orbpx)
      med2=median(oflx)
      sig2=robust_sigma(oflx)
      
      xx=where(orbpx lt med1-3*sig1 OR orbpx gt med1+3*sig1 OR oflx lt med2-3*sig2 OR oflx gt med2+3*sig2, nbad)
      
      if nbad gt 0 then good_bad[iloc[xx]]=1
      
    endfor
    
    bad=where(good_bad eq 1, nbad, complement=good)
    
    n=n_elements(good)
    
    ; keep only good points
    jd=jd[good]
    data1=data1[*,*,good]
    flag=flag[good]
    flux=flux[good]
    npix=npix[good]
    temperature=fltarr(ngood)
    for ii=0, n-1 do temperature=average(ccd_temp[*,good[ii]])
    frame=frame[good]
    xcom=xy_com[0,good]
    ycom=xy_com[1,good]
    data=data1
    time=jd
    resid=resid[good]
    
    fname=file_basename(filesin[ff])
    fileout=outdir+fname
    save, filename=fileout, time, frame, flux, xcom, ycom, temperature, $
      vmag, bmag, resid, psf_aper, data
 
  endfor  ; end loop over file
  
  print, 'end of program'
end


