pro compare_p2_p3

  ; program to check the output images from remove_orbit_trend.pro
  ;
  ; check for any remaining residuals
  ;
  ; check orbit does not have an increasing trend?
  ;
  Compile_opt idl2
  
  !p.background=cgcolor('white')
  
  ; input directory
  indir2='~/BRITE/data/UB/p2/'
  indir3='~/BRITE/data/UB/p3/'
  
  filesin2=file_search(indir2+'*.sav', count=nsav2)
  filesin3=file_search(indir3+'*.sav', count=nsav3)
  
  for i=0, nsav2-1 do begin
  
    restore, filesin2[i] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp
    print, filesin2[i]
    
    window, 0
    plot_image, bytscl(data1[*,*,0], 50, 500)
    print, robust_sigma(data1)
    
    med_frm2=fltarr(n_elements(jd))
    tot_frm2=fltarr(n_elements(jd))
    for jj=0, n_elements(jd)-1 do med_frm2[jj]=median(data1[*,*,jj])
    for jj=0, n_elements(jd)-1 do tot_frm2[jj]=total(data1[*,*,jj])
    
    restore, filesin3[i] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, sm_flx, avg_flx, expn
    print, filesin3[i]
    
    window, 1
    plot_image, bytscl(data1[*,*,0], 50, 500)
    print, robust_sigma(data1)

    
    med_frm3=fltarr(n_elements(jd))
    tot_frm3=fltarr(n_elements(jd))
    for jj=0, n_elements(jd)-1 do med_frm3[jj]=median(data1[*,*,jj])
    for jj=0, n_elements(jd)-1 do tot_frm3[jj]=total(data1[*,*,jj])
    
    stop
 
  endfor
  
  stop
end