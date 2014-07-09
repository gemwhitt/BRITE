pro thermal_tests_plot2

; plot results from thermal_testing2
; 
Compile_opt idl2
!p.background=cgcolor('white')

indir='~/BRITE/results/ThermalTest_240913/hot_pixels/'

outdir='~/BRITE/results/ThermalTest_240913/plots/'

toplot='y'

fileout1=outdir+'hot_pixels.ps'
fileout2=outdir+'hot_pixels2.ps'

filesin=file_search(indir+'*.txt', count=nfiles)

mean_hps=fltarr(nfiles,4)

fname=file_basename(filesin, '.txt')
times=strmid(fname, 7)

for i=0, nfiles-1 do begin
  
  ; read in general info
  readcol, filesin[i], skipline=1, numline=4, temp, mean_data, sigma_data, thresh, min_data, max_data, pc_hp
  
  ; get hp locations
  readcol, filesin[i], skipline=6, numline=10, hp_loc, format='l'
  
  ; get hot pixel intensities 
  readcol, filesin[i], skipline=17, numline=4, hp1, hp2, hp3, hp4, hp5, hp6, hp7, hp8, hp9, hp10
  
  ; all hps
  hps=[[hp1],[hp2],[hp3],[hp4],[hp5],[hp6],[hp7],[hp8],[hp9],[hp10]]
  
  hps_norm=[[hp1/robust_mean(hp1,2)],[hp2/robust_mean(hp2,2)],[hp3/robust_mean(hp3,2)],[hp4/robust_mean(hp4,2)],$
    [hp5/robust_mean(hp5,2)],[hp6/robust_mean(hp6,2)],[hp7/robust_mean(hp7,2)],[hp8/robust_mean(hp8,2)],$
    [hp9/robust_mean(hp9,2)],[hp10/robust_mean(hp10,2)]]
    
    ; get mean of all hps at each temp
  for j=0, 3 do mean_hps[i,j]=robust_mean(reform(hps[j,*]),2)
  
  ; get the range of values for the y-axis
  min_hp=min(hps)
  max_hp=max(hps)
  spr_hp=(max_hp-min_hp)
  
  ; MAKE PLOTS>>>>>>>>
  
  ;cols1=cgcolor(['])
  
  ; 1. plot pixel intensities versus temp
  ;window, 0, xsize=600, ysize=600, xpos=2300, ypos=60 
  ;plot, temp, hps, /ynozero, xtitle='Temperature', title=file_basename(filesin[i],'.txt'), $
  ;  ytitle='Hot Pixel Intensity', color=cgcolor('black'), charsize=0.8, yrange=[min_hp-spr_hp/2.,max_hp+spr_hp/2.], $
  ;  /nodata
  ;  for j=0, 9 do oplot, temp, hps[*,j], color=cgcolor('blue') ;cols1[j]
    
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
       
endfor

for i=0, 4 do mean_hps[i,*]=mean_hps[i,*]/robust_mean(mean_hps[i,*],2)

cols1=cgcolor(['red','orange', 'blue', 'green yellow', 'olive'])

if toplot eq 'n' then window, 0, xsize=600, ysize=600, xpos=2300, ypos=60 else $
  ps_on, fileout1, xsize=17, ysize=17
plot, temp, mean_hps, xtitle='Temperature', title='Mean Hot Pixels', $
ytitle='Hot Pixel Intensity', color=cgcolor('black'), charsize=0.8, /nodata, yrange=[0.2,2], xrange=[-10,40]
for i=0, 4 do oplot, temp, mean_hps[i,*], color=cols1[i], thick=3
al_legend, times, color=cols1, linestyle=0, textcolors=cgcolor('black'), /left, box=1, $
  outline_color=cgcolor('black'), charsize=0.7
  if toplot eq 'y' then ps_off
  
  if toplot eq 'n' then window, 0, xsize=600, ysize=600, xpos=2300, ypos=60 else $
    ps_on, fileout2, xsize=17, ysize=17
  plot, temp, mean_hps, xtitle='Temperature', title='Mean Hot Pixels', $
    ytitle='Hot Pixel Intensity', color=cgcolor('black'), charsize=0.8, /nodata, yrange=[0.7,1.1], xrange=[-10,40]
  for i=0,4 do oplot, temp, mean_hps[i,*], color=cols1[i], thick=3
  al_legend, times, color=cols1, linestyle=0, textcolors=cgcolor('black'), /left, box=1, $
    outline_color=cgcolor('black'), charsize=0.7
  if toplot eq 'y' then ps_off

stop
end 
