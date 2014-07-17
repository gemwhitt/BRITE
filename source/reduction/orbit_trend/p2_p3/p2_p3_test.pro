pro p2_p3_test

  ; modified from p2_p3_new to test is changing the background region changes the results for the calculation of p2_trend
  ; i.e. if the same ring of pixels is chosen

  ; modifed from p2_p3 on Mon 17th Feb - use alternative way of defining background region - USE MASK of 11x11 pixels
  ; around target center - mask this off then select median value of remaining pixels - use this to subtract from image and
  ; store this value with the scatter of this value
  
  ; remove a systematic trend in the data by subtracting a median value of the pixels which do not have the target flux
  ; save the median value removed as a variable in the output file as p2_trend - one for each data point
  ;
  ; Program uses measure_systend.pro
  
  Compile_opt idl2
  !p.background=cgcolor('white')
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  ; input directory
  indir='~/BRITE/data/UB/p2/' ; HP cleaned
  outdir='~/BRITE/data/UB/p3_2/'  ; HP cleaned + p2_trend removed
  
  filesin=file_search(indir+'HD*_p2.sav', count=nsav)
  fname=file_basename(filesin, '_p2.sav')
  
  if nsav eq 0 then stop
  
  for bb=0, nsav-1 do begin
  
    print, bb
    print, 'Start time ', systime()
    
    restore, filesin[bb]  ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    ;simbad_mags, parlax, otype, sptype
    
    ;data1=data1/10.
    
    ; define new variables for saving in output file - normalized time and header-board temps
    nfrm=n_elements(jd)
    
    hdname=fname[bb]
    
    ; define p2_trend array - save this in output file
    p2_trend=fltarr(nfrm)
    sig_trend=fltarr(nfrm)
    
    for rr=0, nfrm-1 do begin
    
      dat=reform(data1[*,*,rr])
      dat2=dat
      
      x1=2
      x2=29
      y1=2
      y2=29
      
      ; mask off target
      temp_dat=dat2
      
      temp_dat[x1:x2,y1:y2]=-9999.
      
      xx=where(temp_dat ne -9999)
      
      x_loc=(array_indices(dat2, xx))[0,*]
      y_loc=(array_indices(dat2, xx))[1,*]
      
      ;      if rr eq 0 then begin
      ;plot_image, bytscl(dat2, 0, 1000)
      ;oplot, x_loc, y_loc, color=cgcolor('orange'), psym=2
      ;stop
      ;      endif
      
      ; calculate the median in pixels away from the target
      p2_trend[rr]=median(dat2[xx])
      sig_trend[rr]=robust_sigma(dat2[xx]/median(dat2[xx]))
      
      next_frame:
    endfor  ; end loop over individual observations
    
    ; subtract the p2_trend from each image/frame and add a median value
    med_trend=median(p2_trend)
    for rr=0, nfrm-1 do data1[*,*,rr]=data1[*,*,rr]-p2_trend[rr];+med_trend
    
    ; save p3 file
    outfile=outdir+fname[bb]+'_p3_2.sav'
    save, filename=outfile, roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
      simbad_mags, parlax, otype, sptype, p2_trend, med_trend, sig_trend
      
    print, 'Start time ', systime()
    
  endfor
  
  print, 'End of Program'
  
end