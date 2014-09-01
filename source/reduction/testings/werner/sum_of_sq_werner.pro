pro sum_of_sq_werner

  ; Get light curves, plots and figures using sum-of-squares method - based on a SNR in EACH image
  
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 0.9
    
  rbin=10.
  
  indir='~/BRITE/TESTSETS/werner4lc/p2/subsets/'
  
  filesin=file_search(indir+'*.sav', count=nfiles)
  
  outdir='~/BRITE/TESTSETS/werner4lc/lcs/'
  
  plotdir='~/BRITE/TESTSETS/werner4lc/p2/subsets/plots/'
  
  for ff=2, nfiles-1 do begin
  
    fname=file_basename(filesin[ff], '.sav')
    
    print, fname
    print, systime()
    
    restore, filesin[ff]  ;jd, data1, ccd_temp, medcol1, medcol2, medimg0, roi_loc, vmag, bmag, flag, xy_com
    
    jd0=jd-jd[0]
    
    dashpos=strsplit(fname, '_')
    hdname=strmid(fname, dashpos[2], dashpos[3]-1-dashpos[2])
     
    nimg=n_elements(jd)
    
    ccd_temp2=fltarr(nimg)
    for im=0, nimg-1 do ccd_temp2[im]=average(ccd_temp[*,im])
    
    ; set up arrays for collecting results
    npix=lonarr(nimg)         ; number of pixels (binned or otherwise) in PSF
    flux=dblarr(nimg)         ; flux array
    max_dn=lonarr(nimg)       ; max pixel number on image
 stop   
    for im=0, nimg-1 do begin
    
      if flag[im] ne 2 then continue
      
      im0=data1[*,*,im]
      
      ; get dimensions of image
      xdim=(size(im0, /dim))[0]
      ydim=(size(im0, /dim))[1]
      
      max_dn[im]=(max(im0))[0]
      
      ; Use psf_loc [*,im] to make a cutout around the PSF
      x1=(psf_loc[0,im]-2 > 0)
      x2=(psf_loc[1,im]+2 < xdim-1)
      y1=(psf_loc[2,im]-2 > 0)
      y2=(psf_loc[3,im]+2 < ydim-1)
      
      im1=im0[x1:x2,y1:y2]  ; smaller cutout around psf
      
      ; check that size of the region is not too small ???
      sim1=size(im1, /dim)
      
      if sim1[0]*sim1[1] le 50 then begin
        flag[im]=0
        goto, next_im
      endif
      
      ; bin the data according to rbin
      imb=rebin_data(im1,rbin)
      
      ; get a value for the total number of "pixels" in imb
      npixel=n_elements(imb)
      
      ;reorder pixels by data number
      sb=size(imb, /dim)            ; size of binned array
      
      im1d=reform(imb, sb[0]*sb[1])
      sort1=reverse(sort(im1d))
      im1d=im1d[sort1]
      
      count=2 ; start on count > 0
      sig1=total(im1d[0:1])*3.25
      snr=0
      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      
;      dat=imb
;    
;      s=size(dat, /dim)
;      
;      ; add a border to the image
;      dat2=lonarr(s[0]+2,s[1]+2)
;      dat2[1,1]=dat
;      
;      thr=50
;      
;      ; use label region to determine number of illuminated regions and number of PSF pixels
;      r2=label_region(dat2 ge thr, /all_neighbors)
;      r2=r2[1:s[0],1:s[1]]                         ; trim off border
;      
;      hist2=histogram(r2, binsize=1, locations=loc2, reverse_indices=ri2)
;      
;      if n_elements(hist2) gt 2 then stop
;
;      ; now overplot 2nd biggest group of pixels
;      ; sort hist1 in decreasing order
;      sort2=reverse(sort(hist2))
;      ind2=ri2[ri2[sort2[1]]:ri2[sort2[1]+1]-1]
;      i2=array_indices(dat, ind2)
;      
;      xi=i2[0,*]
;      yi=i2[1,*]
;      
;      plot_image, bytscl(dat, 20, 200)
;      
;      oplot, xi, yi, color=cgcolor('orange'), psym=2
;      
;      stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      
     ;  plot, [count], [snr], color=cgcolor('black'), psym=8, xrange=[0,80], yrange=[0,20000]
      
      repeat begin
      
        count=count+1
        snr1=snr
        
        sig1=total(im1d[0:count])*3.25
        
        ;noise1=sqrt(im1d[count] + count*10. + (15.^2) + ccd_temp2[im]) ;- original
        noise1=sqrt(im1d[count]*3.25 + ccd_temp2[im]/rbin^2*count + (15.^2)/rbin^2*count) ;- mod
        ;noise1=sqrt(im1d[count] + ((ccd_temp2[im]*2)*count) + (15.^2)*count)
                  ; signal in the pixel + DC + RN + TEMP variation
        
        snr=sig1/noise1
        
    ;     oplot, [count], [snr], color=cgcolor('black'), psym=8
        
      endrep until snr lt snr1
       
      xx=count-1          ; initial aperture
      psfp=sort1[0:xx]
      xflx=(array_indices(imb, psfp))[0,*]
      yflx=(array_indices(imb, psfp))[1,*]
      
      npix[im]=xx
      flux[im]=total(im1d[0:xx])/float(rbin^2)
      
      ;plot_image, bytscl(imb, 20, 200)
      ;oplot, xflx, yflx, color=cgcolor('purple'), psym=8
      ;print, systime()
      ;stop
      next_im:
    endfor  ; end loop over images
    
    good=where(npix ne 0, ngood)
    
    ; keep only good points
    jd=jd[good]
    jd0=jd0[good]
    flux=flux[good]
    npix=npix[good]
    max_dn=max_dn[good]
    ccd_temp2=ccd_temp2[good]
    psf_loc=psf_loc[*,good]
    xy_com=xy_com[*,good]
    
    n=ngood
    
    ; do a sigma clip per orbit
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
    
    ; do clips
    jd=jd[good]
    jd0=jd0[good]
    flux=flux[good]
    npix=npix[good]
    max_dn=max_dn[good]
    ccd_temp2=ccd_temp2[good]
    psf_loc=psf_loc[*,good]
    xy_com=xy_com[*,good]
    
    time1=jd-2451544.5
    
    ; Make Light Curve - save to plots
    lcout=outdir+'plots/'+fname+'_sos_b'+strtrim(fix(rbin),2)+'.ps'
    plotsym, 0, /fill, 0.4
    ps_on, lcout, xsize=16, ysize=10
    plot, time1, flux, color=cgcolor('black'), psym=8, xtitle='JD-JD2000', ytitle='SOS Flux', title=fname, $
      charsize=0.7, /ynozero
    ps_off
   
    bminv=bmag-vmag
    
    ; save result in a .txt file with a header
    txtout=outdir+'lc_txt/'+fname+'_sos_b'+strtrim(fix(rbin),2)+'.txt'
    openw, lun, txtout, /get_lun
    printf, lun, '#name', fname, format='(a10,x,a10)'
    printf, lun, '#Vmag', vmag, format='(a10,x,a10)'
    printf, lun, '#B-V ', bminv, format='(a10,x,a10)'
    printf, lun, '#JD', 'FLUX', 'NPIX', 'MAX_DN', 'CCD_TEMP', 'X_COM', 'Y_COM', $
      'X_MIN_PSF', 'X_MAX_PSF', 'Y_MIN_PSF', 'Y_MAX_PSF', $
      format='(a14,x,a10,x,a10,x,a10,x,a10,x,a10,x,a10,x,a10,x,a10,x,a10,x,a10)'
    for out=0, n-1 do printf, lun, $
      jd[out], flux[out], npix[out], fix(max_dn[out]), ccd_temp2[out], xy_com[0,out], xy_com[1,out], $
      psf_loc[0,out], psf_loc[1,out], psf_loc[2,out], psf_loc[3,out], $
      format='(d14.6,x,d10.2,x,i10,x,i10,x,d10.2,x,d10.4,x,d10.4,x,d10.4,x,d10.4,x,d10.4,x,d10.4)'
    free_lun, lun
    
    
    ; measure RMS values
    jd1=jd0[1:n-1]
    jdiff=jd1-jd0
    gap=where(jdiff gt 0.015, ngap)
    gap=[-1,gap,n-1]
    
    mags=(-2.5)*alog10(flux)
    mags=mags/median(mags)
    
    fl_std=fltarr(ngap+1)
    fl_sig=fltarr(ngap+1)
    duty_cyc=float(n)/float(nimg)*100.
    
    for orbit=0, ngap do begin
    
      iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
      
      ni=n_elements(iloc)
      
      fl_std[orbit]=(stddev(mags[iloc]))^2.
      
      fl_sig[orbit]=robust_sigma(mags[iloc])
      
    endfor
    
    rms_error=(sqrt(total(fl_std)/float(ngap+1)))*1000. ; mmag
    
    sigma=median(fl_sig)*1000.  ; mmag
    
    ; write out
    openw, lun, outdir+'stats/lc_error_b'+strtrim(fix(rbin),2)+'.txt', /get_lun, /append
    printf, lun, fname, vmag, rms_error, sigma, format='(a10,x,f,x,f,x,f)'
    free_lun, lun
    
  endfor  ; end loop over file
  
  print, 'end of program'
end


