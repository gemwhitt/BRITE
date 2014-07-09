pro plot_modelpsf

; plot stats from orion_modelpsf.pro
; 
Compile_opt idl2
plotsym, 0, /fill, 1.2

statfile='~/BRITE/data/UB/p2/ORION/class/stats/psf_time_temp.txt'

readcol, statfile, hdname, vmag, fl1, fl2, fl3, np1, np2, np3, $
  format='(a,f,d,d,d,f,f,f)'
  
  
diff13=fl3-fl1
diff12=fl2-fl1

plot, np1, fl1, color=cgcolor('black'), psym=8
oplot, np2, fl2, color=cgcolor('blue'), psym=8
oplot, np3, fl3, color=cgcolor('green'), psym=8


stop

end