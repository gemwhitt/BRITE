pro p4_gem9

  ; modified from p3_p4_gem8 - produce output figs and stats
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 0.4
  
  indir='/Users/gemmawhittaker/BRITE/data/UB/p2/ORION/class/sav/'
  
  filesin=file_search(indir+'*.sav', count=nfiles)
  
  outdir='/Users/gemmawhittaker/BRITE/data/UB/p4/ORION/'
  
  plotdir='/Users/gemmawhittaker/BRITE/data/UB/p2/ORION/class/plots/lcs/'
  
  ; choose number of apertures
  nap=1
  
  rbin=10.
  
  dates=[0]
  
  for ff=0, nfiles-1 do begin
  
    fname=file_basename(filesin[ff], '_p2.sav')
    
    print, fname
    
    restore, filesin[ff] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
    ;medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, $
    ;parlax, otype, sptype, iflag, ccd_temp2, xy_psf
    
    jd0=jd-jd[0]
  
    flxrms=fltarr(3)
    deltadn=fltarr(3)
    deltatemp=fltarr(3)
    duty=fltarr(3)
    error2=fltarr(3)
    
    hdname=roi_name[0]
    
    for ll=0, n_elements(dates)-1 do begin

    sg=where(jd0 ge dates[ll] AND jd0 lt dates[ll]+1000., nsg)
    
    ;get total number of frames in this file
    nfrm=nsg
    
    jd=jd[sg]
    jd1=jd0[sg]
    dat=data1[*,*,sg]
    temps=ccd_temp2[sg]  
    xy_psfs=xy_psf[*,sg]
    iflag0=iflag[sg]
       
    jd2=jd1[1:nsg-1]
    jdiff=jd2-jd1
    ; get number of orbits = ngap+1
    gap=where(jdiff gt 0.015, ngap)
    ; calculate cadence of observations
    cadence=robust_mean(jdiff,2)  ; in days
    cadence_m=cadence*24.*60.       ; in minutes
    cadence_s=cadence*24.*3600.       ; in seconds
    n_img_per_orbit=fix(15./cadence)
    
    ; get x and y pixel dimensions
    xdim=(size(dat, /dim))[0]
    ydim=(size(dat, /dim))[1]
    
    ; set up arrays for collecting results
    nsubpix=lonarr(nfrm)      ; number of subpixels in PSF
    flux=dblarr(nfrm)         ; flux array - using signal2
    nsat=intarr(nfrm)           ; number of saturated pixels - left out of flux determination
    max_dn=lonarr(nfrm)
    
    for img=0, nfrm-1 do begin
    
      if xy_psfs[0,img] le 0 OR xy_psfs[1,img] le 0 OR xy_psfs[0,img] ge 320. OR xy_psfs[1,img] ge 320. then CONTINUE
      
      dat1=dat[*,*,img]
      
      max_dn[img]=(max(dat1))[0]
      
      ; determine the background - away from the PSF center
      ;get_bkgd, dat1,pks1,rbin, bkgd,bk_err
      ;dat1=dat1-bkgd
      
      ; bin the image 10x10 and use xy_psf and the model to make a cutout around the target
      ; make it bigger
      bindat=rebin_data(dat1, rbin)
      
      ; use xy_psf to make image up with model and make cutout
      xc=(size(bindat, /dim))[0]/2.
      yc=(size(bindat, /dim))[1]/2.
      
      ; get shifts
      xsh=xc-xy_psfs[0,img]
      ysh=yc-xy_psfs[1,img]
      
      ; get where model has pixels
      mod_ipix=where(model0 ne 0.0, nmodp)
      xm_ipix=(array_indices(model0, mod_ipix))[0,*]
      ym_ipix=(array_indices(model0, mod_ipix))[1,*]
      
      ; get where image has PSF pixels
      xi_ipix=xm_ipix-xsh
      yi_ipix=ym_ipix-ysh
      
      ; get min and max in x and y
      xmin=min(xi_ipix)
      xmax=max(xi_ipix)
      ymin=min(yi_ipix)
      ymax=max(yi_ipix)
      
      if xmin-(3*rbin) lt 0.0 then x1=0 else x1=xmin-(3*rbin)
      if xmax+(3*rbin) gt (xc*2.)-rbin then x2=(xc*2.)-rbin else x2=xmax+(3*rbin)
      if ymin-(3*rbin) lt 0.0 then y1=0 else y1=ymin-(3*rbin)
      if ymax+(3*rbin) gt (yc*2.)-rbin then y2=(yc*2.)-rbin else y2=ymax+(3*rbin)
      
      cutout1=bindat[x1:x2,y1:y2]
      
     ; ps_on, '~/BRITE/reports/mini_meeting2/eta_ori_psf.ps', xsize=15, ysize=15
     ; plot_image, bytscl(dat1[x1/rbin:x2/rbin,y1/rbin:y2/rbin], 20, 500), charsize=0.7
     ; ps_off
      
     ; stop
      
      ; deal with saturated/nonlinear pixels
      sat=where(dat1 ge 9000, numsat, complement=lin)
      
      ; if numsat gt 0 then stop
      nsat[img]=numsat
      
      ncut=n_elements(cutout1)
      
      ; TAKE THE NEXT BIT FROM THE MODEL PSF PROGRAMS
      
      ;reorder pixels by data number
      xsize=(size(cutout1, /dim))[0]
      ysize=(size(cutout1, /dim))[1]
      data1d=reform(cutout1, xsize*ysize)
      sort1=reverse(sort(data1d))
      data1d=data1d[sort1]/float(rbin^2)
      
      count=2
      sig1=(data1d[0:1])*3.25
      snr=0
      
      repeat begin
      
        count=count+1
        snr1=snr
        
        if count ge nmodp*2.5 then goto, skipover
        
        sig1=total(data1d[0:count])*3.25
        
        noise1=sqrt(data1d[count] + count*(50./float(rbin^2) + ((15.^2)/float(rbin^2))^2. + (50.)/float(rbin^2)))
        
        snr=sig1/noise1
        
      endrep until snr lt snr1
      
      xx=count-1          ; initial aperture
      psfp=sort1[0:xx]
      xflx=(array_indices(cutout1, psfp))[0,*]
      yflx=(array_indices(cutout1, psfp))[1,*]
      
   
         
      nsubpix[img]=xx
      flux[img]=total(data1d[0:xx])
     ;stop
      skipover:
    endfor  ; end loop over images
    
    good=where(nsubpix ne 0, ngood)
    
    ; DO A SIGMA CLIP PER ORBIT!!!!!!!!!
    jd2=jd1[good]
    flux=flux[good]
    nsubpix=nsubpix[good]
    max_dn=max_dn[good]
    temps=temps[good]
    jd=jd[good]
    xy_psfs=xy_psfs[*,good]
    
    jd3=jd2[1:ngood-1]
    jdiff=jd3-jd2
    gap2=where(jdiff gt 0.015, ngap)
    gap2=[-1,gap2,ngood-1]
    
    good_bad=intarr(ngood)
    
    for orbit=0, ngap do begin
    
      iloc=indgen(gap2[orbit+1]-gap2[orbit])+gap2[orbit]+1
      
      if n_elements(iloc) lt 5 then good_bad[iloc]=1
      if n_elements(iloc) lt 5 then CONTINUE
      
      ; nsub pix per orbit
      orbpx=nsubpix[iloc]
      oflx=flux[iloc]
      med1=median(orbpx)
      sig1=robust_sigma(orbpx)
      med2=median(oflx)
      sig2=robust_sigma(oflx)
      
      xx=where(orbpx lt med1-3*sig1 OR orbpx gt med1+3*sig1 OR oflx lt med2-3*sig2 OR oflx gt med2+3*sig2, nbad)
      
      if nbad gt 0 then good_bad[iloc[xx]]=1
      
    endfor
    
    bad=where(good_bad eq 1, nbad, complement=good)
    
    ; do clips
    flux=flux[good]
    mags=-2.5*alog10(flux)
    mags=mags/robust_mean(mags,2)
    jd2=jd2[good]
    jd=jd[good]
    nsubpix=nsubpix[good]
    max_dn=max_dn[good]
    max_dn=max_dn/robust_mean(max_dn,2)
    temps=temps[good]
     xy_psfs=xy_psfs[*,good]
;    stop
    ngood=n_elements(good)
    xpsf=reform(xy_psfs[0,*])
    ypsf=reform(xy_psf[1,*])
    
    ; save result in a .txt file
    txtout=outdir+fname+'_p4b.txt'
    openw, lun, txtout, /get_lun
    for out=0, ngood-1 do printf, lun, jd[out], vmag[0], flux[out], xpsf[out], ypsf[out], format='(d14.5,x,f7.2,x,d10.3,x,f7.1,x,f7.1)'
    free_lun, lun
    
    goto, next_roi
    
    ; measure RMS values
    jd3=jd2[1:ngood-1]
    jdiff=jd3-jd2
    gap2=where(jdiff gt 0.015, ngap)
    gap2=[-1,gap2,ngood-1]
    
    fsig=fltarr(ngap+1)
    tsig=fltarr(ngap+1)
    dnsig=fltarr(ngap+1)
    er2=fltarr(ngap+1)
    
    duty[ll]=ngood/float(nsg)*100.
    
    for orbit=0, ngap do begin
    
      iloc=indgen(gap2[orbit+1]-gap2[orbit])+gap2[orbit]+1
      
      ni=n_elements(iloc)
      
      fsig[orbit]=(stddev(mags[iloc]))^2.
      tsig[orbit]=max(temps[iloc])-min(temps[iloc])
      dnsig[orbit]=max(max_dn[iloc])-min(max_dn[iloc])
      er2[orbit]=robust_sigma(mags[iloc])
      
    endfor
    
    flxrms[ll]=sqrt(total(fsig)/float(ngap+1))
    deltadn[ll]=max(dnsig)
    deltatemp[ll]=max(tsig)
    
    error2[ll]=median(er2)
    stop
    ps_on, '~/BRITE/reports/mini_meeting2/eta_ori_4d.ps', xsize=16, ysize=11
    plot, jd, mags, color=cgcolor('black'), psym=8, /ynozero, xtitle='Julian Date', ytitle='Normalized magnitudes',  $
      charsize=0.7
      ps_off
    
   stop
  endfor
  
  flxrms=flxrms*1000.
  
  deltadn=deltadn*100.

  ; write out 
  openw, lun, '~/Desktop/phot_gem9.txt', /get_lun, /append
  printf, lun, fname, vmag, flxrms, deltatemp, format='(a10,x,f,x,f,x,f,x,f,x,f,x,f,x,f)'
  free_lun, lun
  
  next_roi:
  endfor
  
  
  print, 'end of program'
  stop
end


