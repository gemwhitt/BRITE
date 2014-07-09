pro p3_p4_snrcurves

  ; modified from p3_p4_gem3 
  ; for each target get SNR curves - overplot
  ; plot xx (number of subpixels) versus max DN (DN in brightest pixel) to see how distribution of flux affects the SNR curve
  ; plot DN for first 30 pixels - with scatter (error bars)
  
  ; 31 Jan 2014
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, 0.8, /fill
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  targets=['37128','37742','35468','38771','36486']
  
  indir='~/BRITE/data/UB/p3_2/'
  filesin=file_search(indir+'HD'+targets+'*.sav', count=nfiles)
  fname=file_basename(filesin, '_p3.sav')
  
  outdir='/Users/gemmawhittaker/BRITE/reports/myreports/'
  
  fileout=outdir+'snr_comp.ps'
  
  cols1=cgcolor(['blue', 'orange red', 'purple', 'forest green', 'spring green'])
  
;  ps_on, fileout, xsize=15, ysize=13
;  plot, indgen(8000), indgen(8000), xrange=[0,8000], yrange=[20000,110000], /nodata, color=cgcolor('black'), $
;    charsize=0.7, xtitle='Number of pixels in PSF', ytitle='Cumulative SNR'
    
  for bb=0, nfiles-1 do begin
      
    print, fname[bb]
    
    hdname=fname[bb]
    
    restore, filesin[bb] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend
 ;   stop
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
    
    bkgd_cum=fltarr(43)
      
    for img=0, 0 do begin  ;nfrm-1 do begin  ;nfrm-1 do begin
    
      ; print, img
      
      dat1=dat[*,*,img]
      
      ; find approx center of PSF
      mthr=2.5*median(dat1)
      spr=2
      pks=brite_findpeaks(dat1, minthresh=mthr, spread=spr)
      ;stop
      if pks[0] eq -999.9 then CONTINUE
      
      npks=float(n_elements(pks))/3.
      
      if npks gt 1 then begin
        brightest=(where(pks[2,*] eq max(pks[2,*])))[0]
        pks=pks[*,brightest]
      endif
      ;stop
      ; use centroid function to get new x and y location?            ??????
      ; centroid, dat1, XCEN, YCEN, XY_PEAK=[pks[0],pks[1]], FWHM=1.5, /PEAK_LOC
      
      pks1=round(pks) ; round of pks to nearest whole pixel
      
      ; determine the background - away from the PSF center - subtract this from the image
      rbin=1
      get_bkgd, dat1,pks1,rbin, bkgd,bk_err
      ; subtract the background measurement from the image
      ;dat1=dat1-bkgd
    
      bkgd_cum[img]=bkgd
      
      ; make cutout - centered on PSF center
      if pks1[0]-(9) lt 0 then x1=0 else x1=pks1[0]-(9)
      if pks1[0]+(9) gt xdim-1 then x2=xdim-1 else x2=pks1[0]+(9)
      if pks1[1]-(9) lt 0 then y1=0 else y1=pks1[1]-(9)
      if pks1[1]+(9) gt ydim-1 then y2=ydim-1 else y2=pks1[1]+(9)
      cutout0=dat1[x1:x2,y1:y2]
      
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
        nsat=numsat/(rbin^2)
        ;  print, nsat
        cutout1[sat]=-9999
      endif else nsat=0
      
      ; divide pixel counts by e.g. 64 to add together later for flux values
      cutout2=cutout2/(rbin^2)
      
      ; smooth cutout2
      cutout2=filter_image(cutout2, /smooth, /all_pixels)
      
      ncut=n_elements(cutout1)
      
      ;reorder pixels by data number
      order1=reverse(sort(cutout2))
      temp_arr=order1
      signal1=cutout1[order1]
      signal2=cutout2[order1]
      
      count=2 ; start on count=2, not count=0
      snr1=0
      snr_cum=[intarr(2)]
      sig_cum=[intarr(2)]
      
      repeat begin
      
        count=count+1
                
        snr=snr1
        
        sig=total(signal2[0:count])
        sig_cum=[sig_cum,sig]
        
        ;snr1=sig/sqrt(signal2[count] + (p2_trend[img]/(rbin^2)*count))  ; ORIGINAL
        snr1=sig/sqrt(signal2[count] + (bkgd*bk_err/(rbin^2)*count)) 
        snr_cum=[snr_cum,snr1]
;        stop
      endrep until count eq 8000
            
      xx=where(snr_cum eq max(snr_cum), nxx)
      if nxx gt 1 then xx=xx[nxx-1] ; location of maximum SNR - i.e. number of subpixels
      
     ;xx=xx*2.
;      stop
       temp_arr[0:xx-1]=-9999
       sp=order1[where(temp_arr eq -9999)]
      
      loc_2d=array_indices(cutout2, sp)
      
      wset, 0
      plot_image, bytscl(cutout2*64.,50, 500)
      oplot, loc_2d[0,*], loc_2d[1,*], color=cgcolor('orange'), psym=2
          
      max_snr=snr_cum[xx]
      
      max_dn=signal1[0]
     stop 
      ; plot SNR curve
      wset, 2
      plot, indgen(8000), snr_cum, color=cgcolor('black'), $
        charsize=0.7, xtitle='Number of pixels in PSF', ytitle='Cumulative SNR'
        oplot, indgen(8000), sig_cum, color=cgcolor('blue'), thick=2
        
      
      stop
      
      ;oplot, indgen(8000), , color=cols1[bb], thick=3
      ;oplot, [xx,xx], [0, 25000], color=cgcolor('purple')
      
     ; stop
   
    endfor  ; end loop over images
   
  endfor
  
  al_legend, ['HD37128, V=1.70', 'HD37742, V=1.77', 'HD35468, V=1.64', 'HD38771, V=2.05', 'HD36486, V=2.41'], $
    colors=cols1, textcolor=cgcolor('black'), /right, linestyle=0, charsize=0.65, thick=3
  
  ps_off
  
  spawn, 'open '+fileout+' &'
  
  
  print, 'end of program'
  stop
end


