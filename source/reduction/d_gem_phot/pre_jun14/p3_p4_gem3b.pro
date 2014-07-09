pro p3_p4_gem3b

  ; modified from p3_p4_gem3b - apply a filter
  ; produce p4 light curves
  ; this version treats saturated stars and gives 4 versions of light curve using a ratio of snr/xx
  
  ; 31 Jan 2014
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, 0.8, /fill
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  indir='~/BRITE/data/UB/p3_2/'
  filesin=file_search(indir+'Orion*.sav', count=nfiles)
  fname=file_basename(filesin, '_p3.sav')
  
  outdir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot3b/'
  
  plotdir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot3b/plots/'
  
  for bb=0, nfiles-1 do begin
  
    plotpsf=0
    
    print, fname[bb]
    
    hdname=fname[bb]
    
    restore, filesin[bb] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend, sig_trend
    
    jd1=jd-jd[0]
    keep=where(jd1 ge 0.0)
    
    jd1=jd1[keep]
    data1=data1[*,*,keep]
    roi_name=roi_name[keep]
    exp_num=exp_num[keep]
    ra_dec1=ra_dec1[*,keep]
    ccd_temp=ccd_temp[keep]
    p2_trend=p2_trend[keep]
    
    dat=data1 ; no effect to data1 OR times 10. for 0.1s exp
    
    jd2=jd1[1:n_elements(jd1)-1]
    jdiff=jd2-jd1
    cadence=robust_mean(jdiff,2)  ; in days
    cadence=cadence*24.*60.       ; in minutes
    n_img_per_orbit=fix(15./cadence)
    
    ; print the magnitude
    vmag=strmid(simbad_mags, strpos(simbad_mags, 'V=')+2)
    print, 'vmag='+vmag
    
    ;get total number of frames in this file
    nfrm=(size(dat, /dim))[2]
    
    ; get x and y pixel dimensions
    xdim=(size(dat, /dim))[0]
    ydim=(size(dat, /dim))[1]
    
    ; set up arrays for collecting results
    nsubpix=intarr(4,nfrm)      ; number of subpixels in PSF
    max_snr=fltarr(nfrm)        ; SNR at nsubpix
    max_dn=fltarr(nfrm)         ; maximum data number on image
    xy_psf=fltarr(2,nfrm)       ; x and y coords of PSF center
    flux=fltarr(4,nfrm)         ; flux array - using signal2
    nsat=intarr(nfrm)           ; number of saturated pixels - left out of flux determination
    
    bkgd_cum=fltarr(43)
    
    for img=0, 40 DO BEGIN  ;nfrm-1 do begin
    
     ; print, img
      
      dat1=dat[*,*,img]
      
      ; filter the image - gaussian
      ;dat1=filter_image(dat1, fwhm=1.5, /all_pixels)
      
      ; filter image in 1d
      dat_1d=reform(dat1, xdim*ydim)
      sort1=reverse(sort(dat_1d))
      dat_1d=dat_1d[sort1]
      
      ;wset, 2
      if img eq 0 then plot, dat_1d, xrange=[-1,20], yrange=[0,2200], psym=2, color=cgcolor('black') else $
        oplot, dat_1d, psym=2, color=cgcolor('black')
      
      cgminmax, dat1
      
      ; apply a filter...
      sm1=smooth(dat_1d, 3., /edge_truncate)
      
      cgminmax, sm1
      
      ;oplot, sm1, psym=4, color=cgcolor('red')
      
      ; put pixels back in original order
      datsm=dat_1d
      datsm[sort1]=sm1
      ; put pixels back in 2d-array
      sm_2d=reform(datsm, xdim, ydim)

      ;wset, 0
      ;plot_image, bytscl(dat1, 50 , 500)
      
      ;wset, 1
      ;plot_image, bytscl(sm_2d, 50 , 500)
      

      continue
      
     
      
      ; find approx center of PSF
      mthr=2.5*median(dat1)
      spr=2
      pks=brite_findpeaks(dat1, minthresh=mthr, spread=spr)
      
      if pks[0] eq -999.9 then CONTINUE
      if pks[0] lt 0 OR pks[1] lt 0 then CONTINUE
      if finite(pks[0]) eq 0 or finite(pks[1]) eq 0 then continue
      
      npks=float(n_elements(pks))/3.
      
      if npks gt 1 then begin
        brightest=(where(pks[2,*] eq max(pks[2,*])))[0]
        pks=pks[*,brightest]
      endif
      ;stop
      ; use centroid function to get new x and y location?            ??????
      ; centroid, dat1, XCEN, YCEN, XY_PEAK=[pks[0],pks[1]], FWHM=1.5, /PEAK_LOC
      
      xy_psf[0,img]=pks[0]  ;xcen
      xy_psf[1,img]=pks[1]  ;ycen
      
      pks1=round(pks) ; round of pks to nearest whole pixel
      
      ; determine the background - away from the PSF center - subtract this from the image
      rbin=1
      get_bkgd, dat1,pks1,rbin, bkgd,bk_err
      ; subtract the background measurement from the image
      
      ;dat1=dat1-bkgd
      
      ; make cutout - centered on PSF center
      if pks1[0]-(10) lt 0 then x1=0 else x1=pks1[0]-(10)
      if pks1[0]+(10) gt xdim-1 then x2=xdim-1 else x2=pks1[0]+(10)
      if pks1[1]-(10) lt 0 then y1=0 else y1=pks1[1]-(10)
      if pks1[1]+(10) gt ydim-1 then y2=ydim-1 else y2=pks1[1]+(10)
      cutout0=dat1[x1:x2,y1:y2]
;      if img eq 377 then stop
      ; rebin the cutout - 4x4, 8x8, 16x16...
      rbin=8. ;8.
      cutout1=rebin_data(cutout0,rbin)
      
      cutout2=cutout1
      ;stop
      for i=0, (size(cutout0, /dim))[0]-1 do begin      ; column
        for j=0, (size(cutout0, /dim))[1]-1 do begin    ; row
        
          cutout2[i*rbin:((i+1)*rbin)-1,j*rbin:((j+1)*rbin)-1]=cutout0[i,j] ; this array is equivalent in dimension to cutout
          
        endfor
      endfor
      
      ; deal with saturated pixels
      sat=where(cutout2 ge 9000, numsat, complement=lin)
      
      if numsat gt 0 then begin
        ;stop
        nsat[img]=numsat/(rbin^2)
        ;  print, nsat
        cutout1[sat]=-9999
      endif else nsat[img]=0
      
      ; divide pixel counts by e.g. 64 to add together later for flux values
      cutout2=cutout2/(rbin^2)
      
      ncut=n_elements(cutout1)
      
      ;reorder pixels by data number
      order1=reverse(sort(cutout1))
      temp_arr=order1
      signal1=cutout1[order1]
      signal2=cutout2[order1]
      
      count=1 ; start on count=2, not count=0
      snr1=0
      snr_cum=[signal1[0],total(signal1[0:1])]
      sig_cum=[snr_cum[0]/sqrt(signal1[0]), snr_cum[1]/sqrt(signal1[1] + (p2_trend[img]/(rbin^2)))]
      
      repeat begin
      
        count=count+1
        
        if count gt 6000 then goto, skipover
        
        snr=snr1
        
        sig=total(signal1[0:count])
        sig_cum=[sig_cum,sig]
        
        ;snr1=sig/sqrt(signal1[count] + (bkgd/(rbin^2)*count))  ; original
        snr1=sig/sqrt(signal1[count] + (p2_trend[img]/(rbin^2)*count))  ; using p2_trend - different depending on temp
        snr_cum=[snr_cum,snr1]
        
      endrep until snr1 lt snr
      
      if count lt 400 then goto, skipover
      
      xx=count-1 ; location of maximum SNR - i.e. number of subpixels
      
      max_snr[img]=snr_cum[xx]
      
      ; for plotting subpixels
      loc_2d=array_indices(cutout1, order1[0:xx])
      resid=array_indices(cutout1, order1[xx+1:ncut-1])
      
      xpix=loc_2d[0,*]
      ypix=loc_2d[1,*]
      
      ;window, 0, xsize=500, ysize=500, xpos=1500, ypos=200
      
      ; save other variables
      nsubpix[0,img]=xx
      
      flux[0,img]=total(signal2[0:xx])
      
      max_dn[img]=signal1[0]
      ;      stop
      ; GET RATIO - max_snr/nsubpix
      ratio=max_snr[img]/float(xx)
      scale1=1.1  ; 10%  
      scale2=1.2  ; 20%
      scale3=1.3  ; 30%
      
      ; get 3 more locations for cutoffs of nsubpix
      yy1=round(scale1*xx)
      yy2=round(scale2*xx)
      yy3=round(scale3*xx)
      
      nsubpix[1,img]=yy1
      nsubpix[2,img]=yy2
      nsubpix[3,img]=yy3
      ; stop
      flux[1,img]=total(signal2[0:yy1])
      flux[2,img]=total(signal2[0:yy2])
      flux[3,img]=total(signal2[0:yy3])
      
      ; for plotting.....
      if plotpsf eq 0 then begin
      
        plotout=plotdir+hdname+'_psf.ps'
        !p.multi=[0,2,2,0,0]
        ps_on, plotout, xsize=18, ysize=18
        
        for gem=0, 3 do begin
        
          loc_2d=array_indices(cutout1, order1[0:nsubpix[gem,img]])
          
          xpix=loc_2d[0,*]
          ypix=loc_2d[1,*]
          
          plot_image, bytscl(cutout1, 0, 500), color=cgcolor('black'), title=hdname+',V='+strtrim(vmag,2)+', scale= 0-500', charsize=0.5
          oplot, xpix, ypix, psym=8, color=cgcolor('purple')
          
        endfor
        
        ps_off
        
        ;spawn, 'open '+plotout+' &'
        ;stop
        plotpsf=1
        !p.background=cgcolor('white')
        ; stop
      endif
   STOP   
      skipover:
      STOP
    endfor  ; end loop over images
    stop
    ; reject points with really low or really high: nsubpix, max_dn, max_snr....
    mean_count=robust_mean(nsubpix[0,*],2)
    mean_dn=robust_mean(max_dn,2)
    mean_snr=robust_mean(max_snr,2)
    
    rej=where(reform(nsubpix[0,*]) le mean_count*0.5 OR reform(nsubpix[0,*]) ge mean_count*1.5 OR $
      max_dn le mean_dn*0.5 OR max_dn ge mean_dn*1.5 OR $
      max_snr lt mean_snr*0.5 OR max_snr ge mean_snr*1.5, nrej, complement=keep)
    ;stop
    !p.multi=[0,1,3,0,0]
    !p.background=cgcolor('white')
    plotsym, 0, /fill, 0.3
    plotout=plotdir+hdname+'_clips.ps'
    ps_on, plotout, xsize=16, ysize=25
    
    ;window, 0, xsize=700, ysize=400
    plot, jd1, flux[0,*], color=cgcolor('black'), psym=8, title=hdname+', Vmag='+strtrim(vmag,2), xtitle='Time (days)', $
      ytitle='Flux (DN)', charsize=0.8, xrange=[min(jd1)-0.1, max(jd1)+0.1], xstyle=1
    if nrej gt 0 then  oplot, jd1[rej], flux[0,rej], color=cgcolor('red'), psym=8
    
    ;window, 1, xsize=700, ysize=400
    plot, jd1[keep], flux[0,keep], color=cgcolor('black'), psym=8, title='After clips', xtitle='Time (days)', $
      ytitle='Flux (DN)', charsize=0.8, /ynozero, xrange=[min(jd1)-0.1, max(jd1)+0.1], xstyle=1
      
    mags=-2.5*alog10(reform(flux[0,keep]))
    plot, jd1[keep], mags, color=cgcolor('black'), psym=8, title='After clips', xtitle='Time (days)', $
      ytitle='Flux (Mags)', charsize=0.8, /ynozero  , xrange=[min(jd1)-0.1, max(jd1)+0.1], xstyle=1
      
    ps_off
    
    ;     spawn, 'open '+plotout+' &'
    ;stop
    ; modify results - rejecting bad points...
    nsubpix=nsubpix[*,keep]
    max_snr=max_snr[keep]
    max_dn=max_dn[keep]
    xy_psf=xy_psf[*,keep]
    flux=flux[*,keep]
    jd1=jd1[keep]
    ccd_temp=ccd_temp[keep]
    roi_name=roi_name[keep]
    exp_num=exp_num[keep]
    ra_dec1=ra_dec1[*,keep]
    p2_trend=p2_trend[keep]
    nsat=nsat[keep]
    
    jd=jd1
    ;stop
    ; save results
    fileout=outdir+fname[bb]+'_p4_gem3b.sav'
    ;    stop
    save, filename=fileout, vmag, nsubpix, max_snr, max_dn, xy_psf, flux, jd, nsat, $
      roi_name, exp_num, ra_dec1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
      simbad_mags, parlax, otype, sptype, p2_trend, med_trend, exp_ttl, exp_time
    ; stop
  endfor
  
  
  print, 'end of program'
  stop
end


