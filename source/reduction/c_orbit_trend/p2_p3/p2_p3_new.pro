pro p2_p3_new

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
  indir='~/BRITE/data/UB/testing/temp_effects/p2/' ; HP cleaned
  outdir='~/BRITE/data/UB/testing/temp_effects/p3/'  ; HP cleaned + p2_trend removed
  
  filesin=file_search(indir+'Orion-CF1-*_p2.sav', count=nsav)
  fname=file_basename(filesin, '_p2.sav')
    
  if nsav eq 0 then stop
  
  for bb=0, nsav-1 do begin
  
    print, bb
    print, 'Start time ', systime()
    
    restore, filesin[bb]  ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
                          ;simbad_mags, parlax, otype, sptype, exp_time, exp_ttl
                          
                          ;jd, jd1, data1, ccd_temp, roi_dim, simbad_mags  - testsets
                          
;                          stop
    
   ; stacked=exp_ttl/exp_time
    
    ;data1=data1/stacked
    
 ;   if bb mod 15 eq 0 then stop    
    ; define new variables for saving in output file - normalized time and header-board temps
    nfrm=n_elements(jd)
    
    hdname=fname[bb]
    
    ; define p2_trend array - save this in output file
    p2_trend=fltarr(nfrm)
    sig_trend=fltarr(nfrm)
    
    for rr=0, nfrm-1 do begin
    
      dat=reform(data1[*,*,rr])
      dat2=dat
       
      ; find approx center of PSF
      mthr=1.5*robust_mean(dat2,3)
      spr=2
      pks=brite_findpeaks(dat2, minthresh=mthr, spread=spr)
      
      if pks[0] eq -999.9 then goto, next_frame
      
      npks=n_elements(pks)/3.
      
      if npks gt 1 then pks=pks[*,where(pks[2,*] eq max(pks[2,*]))]
      
      pks1=round(pks) ; center of target PSF
      
      ; get image dimensions
      xdim=(size(dat2, /dim))[0] 
      ydim=(size(dat2, /dim))[1]
      
      ; check target proximity to edge of frame
      if pks1[0]-10 lt 0 then x1=0 else x1=pks1[0]-10
      if pks1[0]+10 gt xdim-1 then x2=xdim-1 else x2=pks1[0]+10
      if pks1[1]-10 lt 0 then y1=0 else y1=pks1[1]-10
      if pks1[1]+10 gt ydim-1 then y2=ydim-1 else y2=pks1[1]+10
      
      ; mask off target        
      temp_dat=dat2
               
      temp_dat[x1:x2,y1:y2]=-9999.
      
      xx=where(temp_dat ne -9999)
        
      x_loc=(array_indices(dat2, xx))[0,*]
      y_loc=(array_indices(dat2, xx))[1,*]
        
    ;  if rr eq 0 then begin
    ;    window, 0, xsize=500, ysize=500, xpos=900, ypos=600
     ;   fileout='~/BRITE/reports/myreports/p2_backgd_image.ps'
     ;   ps_on, fileout, xsize=15, ysize=15
    ;    plot_image, bytscl(dat2, 0, 1000)
    ;    oplot, x_loc, y_loc, color=cgcolor('orange'), psym=2
     ;   ps_off
     ;   spawn, 'open '+fileout
     ;   stop
     ; endif
              
      ; calculate the median in pixels away from the target
      p2_trend[rr]=robust_mean(dat2[xx],2)
      sig_trend[rr]=robust_sigma(dat2[xx]/p2_trend[rr])
      
      next_frame:
    endfor  ; end loop over individual observations
    
    ; subtract the p2_trend from each image/frame and add a median value
    med_trend=median(p2_trend)
    for rr=0, nfrm-1 do begin
      
      data1[*,*,rr]=data1[*,*,rr]-p2_trend[rr]          ;+med_trend
      
      ; find pixels with DN < 0.
      data2=data1[*,*,rr]
      xx=where(data2 lt 0.0, nxx)
      
      if nxx gt 0 then begin
        
        xx_2d=array_indices(data2, xx)
        xloc=xx_2d[0,*]
        yloc=xx_2d[1,*]
        
        data2[xloc, yloc]=0.0
        
        data1[*,*,rr]=data2
        
      endif
     
    endfor
    
    ; save p3 file
    outfile=outdir+fname[bb]+'_p3.sav'
  ;  save, filename=outfile, roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
  ;    simbad_mags, parlax, otype, sptype, exp_time, exp_ttl, p2_trend, med_trend, sig_trend
      
     save, filename=outfile, jd, jd1, data1, ccd_temp, roi_dim, simbad_mags, p2_trend, med_trend
      
    print, 'Start time ', systime()
    
  endfor
  
  print, 'End of Program'
  
end