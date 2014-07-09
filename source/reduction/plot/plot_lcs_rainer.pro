pro plot_lcs_rainer

; double check the data for Rainer - sent December 15th

  compile_opt idl2
  
  !p.background=cgcolor('white')
  
  indir='/Users/gemmawhittaker/BRITE/data/UB/p4_rainer/'
  
  filesin=file_search(indir+'*_gem2.sav', count=nsav)
  fname=file_basename(filesin, '_p4_gem2.sav')
  
  for ii=0, nsav-1 do begin
  
    restore, filesin[ii]  ; flux, jd, bkgd, $
      ;exp_num, roi_dim, ccd_temp, xy_psf
    
    window, 0, xpos=100, ypos=500, xsize=600, ysize=400
    plot, jd, flux, psym=2, color=cgcolor('black'), /ynozero
    
    window, 1, xpos=700, ypos=500, xsize=600, ysize=400
    plot, jd, bkgd, psym=2, color=cgcolor('black'), /ynozero
   
    
    window, 3, xpos=100, ypos=10, xsize=600, ysize=400
    plot, xy_psf[0,*], xy_psf[1,*], psym=2, color=cgcolor('black'), /ynozero
    
    stop
    
   
    
  endfor
  
  print, 'End of program'
end