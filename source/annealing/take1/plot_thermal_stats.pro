pro plot_thermal_stats

Compile_opt idl2

!p.background=cgcolor('white')

; plot stats from thermal_testing2.pro
; 
indir='~/BRITE/results/ThermalTest_240913/hot_pixels/'

outdir='~/BRITE/results/ThermalTest_240913/plots/hot_pixels/'

fileout1=outdir+'mean_temp.ps'
fileout2=outdir+'sig_temp.ps'
fileout3=outdir+'pc_hp_temp.ps'

toplot='y'

filesin=file_search(indir+'*.txt', count=nfiles)

mean_int=fltarr(nfiles,4)
sig_int=fltarr(nfiles,4)
pc_hp=fltarr(nfiles,4)
times=strarr(nfiles)

temp=[0.,10.,20.,30.]

for i=0, nfiles-1 do begin
  
  readcol, filesin[i], temp1, mean1, sigma1, thresh, min1, max1, pc_hp1, skipline=1, numline=4
  
  mean_int[i,*]=mean1
  sig_int[i,*]=sigma1
  pc_hp[i,*]=pc_hp1
  times[i]=strmid(file_basename(filesin[i],'.txt'),7)

endfor

; need 5 symbols and colours
plotsym, 0, 0.3, /fill
cols1=cgcolor(['blue','green','orange','purple','red'])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot medians
; get range of medians for plotting
min_med=min(mean_int)
max_med=max(mean_int)

if toplot eq 'n' then window, 0, xsize=500, ysize=500 else $
  ps_on, fileout1, xsize=17, ysize=17
  plot, temp, mean_int[0,*], psym=-8, color=cgcolor('black'), $
    yrange=[min_med,max_med], xrange=[-5., 35.], xstyle=1, charsize=0.7, $
    xtitle='Temperature', ytitle='Data median', title='Data median vs temperature - all times', /nodata
for i=0, nfiles-1 do oplot, temp, mean_int[i,*], psym=-8, color=cols1[i], thick=3
  al_legend, times, psym=[8,8,8,8,8], color=cols1, textcolors=cgcolor('black'), /left, box=1, $
  outline_color=cgcolor('black'), charsize=0.7
   if toplot eq 'y' then ps_off

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
   ; plot sigmas
   ; get range of sigs for plotting
   min_sig=min(sig_int)
   max_sig=max(sig_int)
   
   if toplot eq 'n' then window, 0, xsize=500, ysize=500 else $
     ps_on, fileout2, xsize=17, ysize=17
   plot, temp, sig_int[0,*], psym=-8, color=cgcolor('black'), $
     yrange=[min_sig,max_sig], xrange=[-5., 35.], xstyle=1, charsize=0.7, $
     xtitle='Temperature', ytitle='Data scatter', title='Data scatter vs temperature - all times', /nodata
   for i=0, nfiles-1 do oplot, temp, sig_int[i,*], psym=-8, color=cols1[i], thick=3
   al_legend, times, psym=[8,8,8,8,8], color=cols1, textcolors=cgcolor('black'), /left, box=1, $
     outline_color=cgcolor('black'), charsize=0.7
      if toplot eq 'y' then ps_off
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     ; plot pc hot pixels
   ; get range of medians for plotting
   min_hp=min(pc_hp)
   max_hp=max(pc_hp)


   if toplot eq 'n' then window, 0, xsize=500, ysize=500 else $
     ps_on, fileout3, xsize=17, ysize=17
   plot, temp, pc_hp[0,*], psym=-8, color=cgcolor('black'), $
     yrange=[min_hp,max_hp], xrange=[-5., 35.], xstyle=1, charsize=0.7, $
     xtitle='Temperature', ytitle='Percentage of Hot pixels (>180)', title='Percentage of hot pixels vs temperature - all times', /nodata
   for i=0, nfiles-1 do oplot, temp, pc_hp[i,*], psym=-8, color=cols1[i], thick=3
   al_legend, times, psym=[8,8,8,8,8], color=cols1, textcolors=cgcolor('black'), /left, box=1, $
     outline_color=cgcolor('black'), charsize=0.7
      if toplot eq 'y' then ps_off
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  
print, 'end of program'

end
