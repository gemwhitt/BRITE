pro p4_gem8

  ; modified from p3_p4_gem7 - for Orion
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 0.8
    
  indir='/Users/gemmawhittaker/BRITE/data/UB/p2/ORION/class/sav/'
  
  filesin=file_search(indir+'*.sav', count=nfiles)
  
  outdir='/Users/gemmawhittaker/BRITE/data/UB/p2/ORION/class/lcs/'
  
  plotdir='/Users/gemmawhittaker/BRITE/data/UB/p2/ORION/class/plots/lcs/'
  
  ; choose number of apertures
  nap=1
  
  rbin=10.
   
  for ff=0, nfiles-1 do begin
    
    fname=file_basename(filesin[ff], '_p2.sav')
    
    print, fname
    
    restore, filesin[ff] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
                         ;medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, $
                         ;parlax, otype, sptype, iflag, ccd_temp2, xy_psf
    
    hdname=roi_name[0]
    
    ;get total number of frames in this file
    nfrm=n_elements(jd)
    
    jd1=jd-jd[0]  
    jd2=jd1[1:n_elements(jd1)-1]
    jdiff=jd2-jd1
    ; get number of orbits = ngap+1
    gap=where(jdiff gt 0.015, ngap)
    ; calculate cadence of observations
    cadence=robust_mean(jdiff,2)  ; in days
    cadence_m=cadence*24.*60.       ; in minutes
    cadence_s=cadence*24.*3600.       ; in seconds
    n_img_per_orbit=fix(15./cadence)
    
    ; get x and y pixel dimensions
    xdim=(size(data1, /dim))[0]
    ydim=(size(data1, /dim))[1]
    
    ; set up arrays for collecting results
    nsubpix=lonarr(nfrm)      ; number of subpixels in PSF
    flux=dblarr(nfrm)         ; flux array - using signal2
    nsat=intarr(nfrm)           ; number of saturated pixels - left out of flux determination
    max_dn=lonarr(nfrm)
    ratio1=fltarr(nfrm)
    
    for img=0, 2000 do begin ;nfrm-1 do begin
      
      if xy_psf[0,img] le 0 OR xy_psf[1,img] le 0 OR xy_psf[0,img] ge 320. OR xy_psf[1,img] ge 320. then CONTINUE
    
      dat1=data1[*,*,img]
      
      maxdn=(max(dat1))[0]
            
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
      xsh=xc-xy_psf[0,img]
      ysh=yc-xy_psf[1,img]
      
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
      
 ;     model2=model0[x1+xsh:x2+xsh,y1+ysh:y2+ysh]
      
;      ratio1[img]=total(cutout1)/total(model2)*100.

      ;cutout2=cutout1 
      ;for i=0, (size(cutout0, /dim))[0]-1 do begin      ; column
      ;  for j=0, (size(cutout0, /dim))[1]-1 do begin    ; row
        
      ;    cutout2[i*rbin:((i+1)*rbin)-1,j*rbin:((j+1)*rbin)-1]=cutout0[i,j] ; this array is equivalent in dimension to cutout
          
      ;  endfor
      ;endfor
      
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
      
      if count ge nmodp*2 then goto, skipover
    
      sig1=total(data1d[0:count])*3.25
      
      noise1=sqrt(data1d[count] + count*(50./float(rbin^2) + ((15.^2)/float(rbin^2))^2. + (50.)/float(rbin^2))) 
      
      snr=sig1/noise1
       
  endrep until snr lt snr1
      
  xx=count-1          ; initial aperture
  psfp=sort1[0:xx]
  xflx=(array_indices(cutout1, psfp))[0,*]
  yflx=(array_indices(cutout1, psfp))[1,*]
  
  ;plot_image, bytscl(cutout1, 20, 500)
  ;stop
  ;oplot, xflx, yflx, color=cgcolor('purple'), psym=2
  ;stop
  
  nsubpix[img]=xx
  flux[img]=total(data1d[0:xx])
  
  ; plot images    
  ; if nap gt 1 then stop
  ;         
  ;    if plotpsf eq 0 then begin
  ;    
   ;     plotout=plotdir+fname+'_aper.ps'
   ;     if nap gt 1 then !p.multi=[0,2,2,0,0] else !p.multi=0
   ;     ps_on, plotout, xsize=18, ysize=18
   ;     
   ;     loc2d=array_indices(cutout1, order1[0:xx])
   ;     xpix=loc2d[0,*]
   ;     ypix=loc2d[1,*]
        
   ;     plot_image, bytscl(cutout1, 0, 500), color=cgcolor('black'), $
   ;       title=hdname+',V='+strtrim(vmag,2)+', NPIX='+strtrim(float(nsubpix[0,img])/(rbin^2),2), charsize=0.5
   ;     oplot, xpix, ypix, psym=8, color=cgcolor('purple')
        
   ;     ps_off
        
   ;     plotpsf=1
   ;     !p.background=cgcolor('white')
   ;   endif
   ;   
   skipover:
  endfor  ; end loop over images
  
  good=where(nsubpix ne 0, ngood)
  
  ; DO A SIGMA CLIP PER ORBIT!!!!!!!!!
  jd2=jd1[good]
  flux2=flux[good]
  nsubpix2=nsubpix[good]
  
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
    orbpx=nsubpix2[iloc]
    oflx=flux2[iloc]
    med1=median(orbpx)
    sig1=robust_sigma(orbpx)
    med2=median(oflx)
    sig2=robust_sigma(oflx)
    
    xx=where(orbpx lt med1-3*sig1 OR orbpx gt med1+3*sig1 OR oflx lt med2-3*sig2 OR oflx gt med2+3*sig2, nbad)
    
    if nbad gt 0 then good_bad[iloc[xx]]=1
  
  endfor
  
  bad=where(good_bad eq 1, nbad)
  
  wset, 0
  plot, jd2, flux2, color=cgcolor('black'), psym=8, /ynozero
  oplot, jd2[bad], flux2[bad], color=cgcolor('red'), psym=8
  stop
  
continue    
   stop   
    
    nsat=robust_mean(nsat,2)
    
    exp_time=float(1)
    
    bminv=bmag-vmag
    fwhm=[fwhmx,fwhmy]
    
    ; save .txt file
   ; openw, lun, txtout, /get_lun
    printf, lun, hdname, format='(a)'
    printf, lun, 'UB', format='(a)'
    printf, lun, fname, format='(a)'
    printf, lun, 'Processed by GNW', format='(a)'
    printf, lun, 'Version '+version, format='(a)'
    printf, lun, 'Notes:', format='(a)'
    printf, lun, 'Vmag='+strtrim(vmag,2), format='(a)'
    printf, lun, 'BminV='+strtrim(bminv,2), format='(a)'
    printf, lun, 'Num saturated pixels'+strtrim(nsat,2), format='(a)'
    printf, lun, 'Best aperture: '+strtrim(xx,2), format='(a)'
    printf, lun, 'Number of pixels in aperture: '+strtrim(npix_in_ap,2)+' +/- '+strtrim(err_npix,2)
    printf, lun, 'Average FWHM: '+strtrim(fwhm[0],2)+', '+strtrim(fwhm[1],2)+' pixels'
    printf, lun, 'Average sigma across each orbit is: '+strtrim(min(sigflux),2)+'%', format='(a)'
    printf, lun, 'Data below: JD, Magnitudes, DN, ExpTime, CCDtemp, x_psf, y_psf, max_dn, medimg', format='(a)'
    for i=0, n-1 do printf, lun, jd[i], mags[i], flux[i], exp_time, ccd_temp[i], xy_psf[0,i], xy_psf[1,i], max_dn[i], medimg[i], $
      format='(d20.8, x, f, x, f, x, f, x, f7.3, x, f, x, f, x, f, x, f)'
    free_lun, lun
    
    
  endfor
  
  
  print, 'end of program'
  stop
end


