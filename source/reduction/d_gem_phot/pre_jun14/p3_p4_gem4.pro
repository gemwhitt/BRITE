pro p3_p4_gem4

  ; modified from p3_p4_gem3 - change from "interpolating" the data (with rebin) to just dividing up
  ; pixels (8x8)/64 and smoothing the image (filter_image)
  ; 
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
  
  filesin=file_search(indir+'HD*.sav', count=nfiles)
  fname=file_basename(filesin, '_p3.sav')
  
  outdir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot4/'
  
  plotdir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot4/plots/'
  
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

    for img=0, nfrm-1 do begin
    
      dat1=dat[*,*,img]
      
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
      
      ; use centroid function to get new x and y location?            ??????
      ; centroid, dat1, XCEN, YCEN, XY_PEAK=[pks[0],pks[1]], FWHM=1.5, /PEAK_LOC
      
      xy_psf[0,img]=pks[0]  ;xcen
      xy_psf[1,img]=pks[1]  ;ycen
      
      pks1=round(pks) ; round of pks to nearest whole pixel
      
      ; determine the background - away from the PSF center - subtract this from the image
      rbin=1
      ;get_bkgd, dat1,pks1,rbin, bkgd,bk_err
      ; subtract the background measurement from the image
      ;dat1=dat1-bkgd
      
      ; make cutout - centered on PSF center
      if pks1[0]-(9) lt 0 then x1=0 else x1=pks1[0]-(9)
      if pks1[0]+(9) gt xdim-1 then x2=xdim-1 else x2=pks1[0]+(9)
      if pks1[1]-(9) lt 0 then y1=0 else y1=pks1[1]-(9)
      if pks1[1]+(9) gt ydim-1 then y2=ydim-1 else y2=pks1[1]+(9)
      cutout0=dat1[x1:x2,y1:y2]
      
      ; filter the image here - unless there are saturated pixels
      ; deal with saturated pixels
      sat=where(cutout0 ge 8000, numsat, complement=lin)
      
      if numsat gt 0 then begin
        nsat[img]=numsat
        cutout0[sat]=7999
      endif else nsat[img]=0
      
      ; filter image
      cutoutf=filter_image(cutout0, /smooth, /all_pixels)
      
      ; rebin the cutout - 4x4, 8x8, 16x16...
      rbin=8. ;8.
      cutout1=rebin_data(cutoutf,rbin)
      
      cutout2=cutout1
      ;stop
      for i=0, (size(cutout0, /dim))[0]-1 do begin      ; column
        for j=0, (size(cutout0, /dim))[1]-1 do begin    ; row
        
          cutout2[i*rbin:((i+1)*rbin)-1,j*rbin:((j+1)*rbin)-1]=cutout0[i,j] ; this array is equivalent in dimension to cutout
          
        endfor
      endfor
    
      ; divide pixel counts by e.g. 64 to add together later for flux values
      cutout2=cutout2/(rbin^2)
      cutout2_orig=cutout2  ; save a copy of cutout2 that hasn't been filtered
      
      ncut=n_elements(cutout2)
      
      ;reorder pixels by data number
      order1=reverse(sort(cutout2))
      temp_arr=order1
      signal1=cutout1[order1]
      signal2=cutout2[order1]
       
      count=1 ; start on count=2, not count=0
      snr1=0
      sig_cum=[signal2[0],total(signal2[0:1])]
      snr_cum=[sig_cum[0]/sqrt(signal2[0]), sig_cum[1]/sqrt(signal2[1]+(p2_trend[img]/(rbin^2)))]
      
      repeat begin
      
        count=count+1
        
        ;if count gt 5000 then goto, skipover  ;=78 pixels
        
        snr=snr1
        
        sig=total(signal2[0:count])
        sig_cum=[sig_cum,sig]
        
        snr1=sig/sqrt(signal2[count] + (p2_trend[img]/(rbin^2)*count))  ; original
        snr_cum=[snr_cum,snr1]
        
      endrep until snr1 lt snr  ;ORIGINAL
     
      ;if count lt 400 then goto, skipover  ; < 6.25 pixels
      
      xx=count-1 ; location of maximum SNR - i.e. number of subpixels  ;ORIGINAL
     
      x2=xx*1.1
      x3=xx*1.2
      x4=xx*1.4
   
      max_snr[img]=snr_cum[xx]
      max_dn[img]=max(signal2)
      
      ; save other variables
      nsubpix[*,img]=[xx,x2,x3,x4]
      
     if plotpsf eq 0 then begin
      
        plotout=plotdir+hdname+'_psf.ps'
        !p.multi=[0,2,2,0,0]
        ps_on, plotout, xsize=18, ysize=18
        
        for gem=0, 3 do begin
          
          loc_2d=array_indices(cutout2, order1[0:nsubpix[gem,img]])
          
          xpix=loc_2d[0,*]
          ypix=loc_2d[1,*]
          
          plot_image, bytscl(cutout2, 0, 500), color=cgcolor('black'), title=hdname+',V='+strtrim(vmag,2)+', scale= 0-500', charsize=0.5
          oplot, xpix, ypix, psym=8, color=cgcolor('purple')
          
        endfor
       
        ps_off
        
       ; spawn, 'open '+plotout+' &'
        
        plotpsf=1
        !p.background=cgcolor('white')
       ; stop
      endif
        
      flux[*,img]=[total(signal2[0:xx]),total(signal2[0:x2]),total(signal2[0:x3]),total(signal2[0:x4])]
           
      skipover:
    endfor  ; end loop over images
    
    ; remove 0's
    keep=where(nsubpix[0,*] ne 0.0, nkeep)
    nsubpix=nsubpix[*,keep]
    flux=flux[*,keep]
    max_snr=max_snr[keep]
    max_dn=max_dn[keep]
    xy_psf=xy_psf[*,keep]
    jd1=jd1[keep]
    ccd_temp=ccd_temp[keep]
    roi_name=roi_name[keep]
    exp_num=exp_num[keep]
    ra_dec1=ra_dec1[*,keep]
    p2_trend=p2_trend[keep]
    nsat=nsat[keep]
    
    ; reject points with really low or really high: nsubpix, max_dn, max_snr....
    mean_count=robust_mean(nsubpix[0,*],2)
    mean_dn=robust_mean(max_dn,2)
    mean_snr=robust_mean(max_snr,2)
    
    rej=where(reform(nsubpix[0,*]) le mean_count*0.5 OR reform(nsubpix[0,*]) ge mean_count*1.5 OR $
      max_dn le mean_dn*0.5 OR max_dn ge mean_dn*1.5 OR $
      max_snr lt mean_snr*0.5 OR max_snr ge mean_snr*1.5, nrej, complement=keep)
  
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
    
   ; spawn, 'open '+plotout+' &'
       
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
    fileout=outdir+fname[bb]+'_p4_gem1.sav'
    ;    stop
    save, filename=fileout, vmag, nsubpix, max_snr, max_dn, xy_psf, flux, jd, nsat, $
      roi_name, exp_num, ra_dec1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
      simbad_mags, parlax, otype, sptype, p2_trend, med_trend
    
  endfor
  
  
  print, 'end of program'
  stop
end


