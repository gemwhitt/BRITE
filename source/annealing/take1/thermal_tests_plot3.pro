pro thermal_tests_plot3

  ; plot results from thermal_testing2
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  
  indir='~/BRITE/results/ThermalTest_240913/hot_pixels/'
  
  outdir='~/BRITE/results/ThermalTest_240913/plots/'
  
  toplot='y'
  
  fileout1=outdir+'hot_pixels_021013.ps'
  fileout2=outdir+'hot_pixels2_021013.ps'
  
  filesin=file_search(indir+'hotpix2_*.txt', count=nfiles)
  
  mean_hps=fltarr(nfiles,4)
  
  fname=file_basename(filesin, '.txt')
  times=strmid(fname, 7)
  
  hps=fltarr(2,4,6)
  
  for i=0, nfiles-1 do begin
  
    ; read in general info
    readcol, filesin[i], skipline=1, numline=4, temp, mean_data, sigma_data, thresh, min_data, max_data, pc_hp
    
    ; get hp locations
    readcol, filesin[i], skipline=6, numline=10, hp_loc, format='l'
    
    ; get hot pixel intensities
    readcol, filesin[i], skipline=13, numline=4, hp1, hp2, hp3, hp4, hp5, hp6
    
    hps[i,*,0]=hp1
    hps[i,*,1]=hp2
    hps[i,*,2]=hp3
    hps[i,*,3]=hp4
    hps[i,*,4]=hp5
    hps[i,*,5]=hp6
    
    endfor
    
    cols1=cgcolor(['blue','red'])
    
    ps_on, '~/Desktop/hps_021013.ps'
    
    ; 1. plot pixel intensities versus temp
    ;window, 0, xsize=600, ysize=600, xpos=2300, ypos=60
    for i=0, 5 do begin
      plot, temp, hps[0,*,i], xtitle='Temperature', title='Hot Pixel '+strtrim(i,2), $
        ytitle='Hot Pixel Intensity', color=cgcolor('black'), charsize=0.8, yrange=[min(hps[*,*,i]), max(hps[*,*,i])], $
        /nodata
        al_legend, ['10s','1s'], color=cols1, linestyle=0
    
      for j=0, 1 do begin
        oplot, temp, hps[j,*,i], color=cols1[j], thick=3
        endfor
        endfor
        
        ps_off
    
    stop
    ;  ; 1. plot pixel intensities versus temp
    ;  window, 0, xsize=600, ysize=600, xpos=2300, ypos=60
    ;  plot, temp, hps_norm, /ynozero, xtitle='Temperature', title=file_basename(filesin[i],'.txt'), $
    ;    ytitle='Hot Pixel Intensity', color=cgcolor('black'), charsize=0.8, yrange=[0.95,1.05], $
    ;    /nodata
    ;  for j=0, 9 do oplot, temp, hps_norm[*,j], color=cgcolor('blue') ;cols1[j]
    
    ;  ; 1. plot pixel intensities versus temp
    ;  window, 0, xsize=600, ysize=600, xpos=2300, ypos=60
    ;  plot, temp, mean_hps, /ynozero, xtitle='Temperature', title=file_basename(filesin[i],'.txt'), $
    ;    ytitle='Hot Pixel Intensity', color=cgcolor('black'), charsize=0.8
    ;
    
  
  
  for i=0, 4 do mean_hps[i,*]=mean_hps[i,*]/robust_mean(mean_hps[i,*],2)
  
  cols1=cgcolor(['red','orange', 'blue', 'green yellow', 'olive'])
  
 ; if toplot eq 'n' then window, 0, xsize=600, ysize=600, xpos=2300, ypos=60 else $
 ;   ps_on, fileout1, xsize=17, ysize=17
 ; plot, temp, mean_hps, xtitle='Temperature', title='Mean Hot Pixels', $
 ;   ytitle='Hot Pixel Intensity', color=cgcolor('black'), charsize=0.8, /nodata, yrange=[0.2,2], xrange=[-10,40]
 ; for i=0, 4 do oplot, temp, mean_hps[i,*], color=cols1[i], thick=3
 ; al_legend, times, color=cols1, linestyle=0, textcolors=cgcolor('black'), /left, box=1, $
 ;   outline_color=cgcolor('black'), charsize=0.7
 ; if toplot eq 'y' then ps_off
  
 
    ps_on, fileout, xsize=17, ysize=17
  plot, temp, mean_hps, xtitle='Temperature', title='Mean Hot Pixels', $
    ytitle='Hot Pixel Intensity', color=cgcolor('black'), charsize=0.8, /nodata, yrange=[0.7,1.1], xrange=[-10,40]
  for i=0,4 do oplot, temp, mean_hps[i,*], color=cols1[i], thick=3
  al_legend, times, color=cols1, linestyle=0, textcolors=cgcolor('black'), /left, box=1, $
    outline_color=cgcolor('black'), charsize=0.7
  if toplot eq 'y' then ps_off
  
  stop
end
