pro plot_rms

Compile_opt idl2

!p.background=cgcolor('white')
plotsym, 0, /fill, 1.2

filein='~/Desktop/phot_gem9.txt'

readcol, filein, fname, vmag, er1, er2, er3, dt1, dt2, dt3,format='(a,f,f,f,f,f,f,f)'

plot, vmag, er1, color=cgcolor('black'), psym=8, /nodata, xtitle='Vmag', ytitle='RMS error', charsize=1.2
oplot, vmag, er1, color=cgcolor('blue'), psym=8
oplot, vmag, er2, color=cgcolor('red'), psym=8
;oplot, vmag, er3, color=cgcolor('green'), psym=8
oplot, [0,5], [1,1], thick=3, color=cgcolor('purple')
al_legend, ['20', '30'], textcolor=['blue','red'], /left, charsize=1.5, color=cgcolor('black')


stop
end