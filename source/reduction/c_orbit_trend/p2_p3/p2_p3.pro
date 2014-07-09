pro p2_p3

; remove a systematic trend in the data by subtracting a median value of the pixels which do not have the target flux
; save the median value removed as a variable in the output file as p2_trend - one for each data point
; 
; Program uses measure_systend.pro 

Compile_opt idl2
;!p.background=cgcolor('white')
;devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
;imagelib    ; adds system variable !IMAGE needed for plotting
  
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
   
  ; define new variables for saving in output file - normalized time and header-board temps
  nfrm=n_elements(jd)
  
  hdname=fname[bb]
  
  ; define p2_trend array - save this in output file
  p2_trend=fltarr(nfrm)
      
  for rr=0, nfrm-1 do begin
    
    dat=reform(data1[*,*,rr])
    rbin=1
    dat2=rebin_data(dat, rbin)
    
  ;  plot_image, dat
    
  ;  stop
      
    ; determine index locations of outside pixels
      
    ; get roi x and y center coords in binned pixels
    roi_xc=float((size(dat2, /dim))[0])/2.
    roi_yc=float((size(dat2, /dim))[1])/2.
        
    ; find approx center of PSF
    mthr=1.5*robust_mean(dat2,3)
    spr=2*rbin
    pks=brite_findpeaks(dat2, minthresh=mthr, spread=spr)
        
    if pks[0] eq -999.9 then goto, next_frame
        
    npks=n_elements(pks)/3.
        
    if npks gt 1 then pks=pks[*,where(pks[2,*] eq max(pks[2,*]))]
        
    pks1=round(pks)
        
    ; determine which quadrent target is in on the CCD or if it is in the center
    xd=pks[0]-roi_xc  ; difference between peak x-loc and center of ROI
    yd=pks[1]-roi_yc  ; difference between peak y-loc and center of ROI
        
    if abs(xd) lt (4.*rbin) AND abs(yd) lt (4.*rbin) then begin ; target is roughly in center
          
      temp_dat=dat2
      temp_dat[roi_xc-(11.*rbin):roi_xc+(10.*rbin),roi_yc-(11.*rbin):roi_yc+(10.*rbin)]=-9999
          
      xx=where(temp_dat ne -9999, nxx)
          
      x_loc=(array_indices(dat2, xx))[0,*]
      y_loc=(array_indices(dat2, xx))[1,*]
          
      ;plot_image, bytscl(dat2, 50, 500)
      ;oplot, x_loc, y_loc, color=cgcolor('orange'), psym=2
      ;stop
      ;wait, 1
          
    endif else begin  ; target is within a quadrent
                        ; define a region based in the other 3 quads
          
      temp_dat=dat2
          
      ; define 4 quads - lower-left, right and upper-left, right
      if xd lt 0. and yd lt 0. then temp_dat[0:roi_xc+(7.*rbin),0:roi_yc+(7.*rbin)]=-9999 ; lower-left
      if xd gt 0. and yd gt 0. then temp_dat[roi_xc-(7.*rbin):(roi_xc*2)-1, 0:roi_yc+(7.*rbin)]=-9999  ; lower-right
      if xd lt 0. and yd gt 0. then temp_dat[0:roi_xc+(7.*rbin),roi_yc-(7*rbin):(roi_yc*2)-1]=-9999  ; upper-left
      if xd gt 0. and yd gt 0. then temp_dat[roi_xc-(7.*rbin):(roi_xc*2)-1,roi_yc-(7*rbin):(roi_yc*2)-1]=-9999  ; upper-right
          
      xx=where(temp_dat ne -9999, nxx)
          
      x_loc=(array_indices(dat2, xx))[0,*]
      y_loc=(array_indices(dat2, xx))[1,*]
          
     plot_image, bytscl(dat2, 50, 500)
      oplot, x_loc, y_loc, color=cgcolor('orange'), psym=2   
      stop 
      ;wait, 1
          
      endelse
        
      ; calculate the median in pixels away from the target
      p2_trend[rr]=median(dat2[xx])
       
      next_frame:
    endfor  ; end loop over individual observations
    
    ; subtract the p2_trend from each image/frame and add a median value
    med_trend=robust_mean(p2_trend,3)
    for rr=0, nfrm-1 do data1[*,*,rr]=data1[*,*,rr]-p2_trend[rr];+med_trend
    
    ; save p3 file
    outfile=outdir+fname[bb]+'_p3.sav'
    save, filename=outfile, roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
  simbad_mags, parlax, otype, sptype, p2_trend, med_trend
  
  print, 'Start time ', systime()
    
  endfor
  
  print, 'End of Program'
  
end