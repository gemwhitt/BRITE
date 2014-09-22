pro medcol_plot, fileout, medcol2, medcol3, im1, i2d, target

totdif=fix(total(medcol2)-total(medcol3))

ps_on, fileout, xsize=27, ysize=12, /landscape
!p.multi=[0,2,1,0,0]
plotsym, 0, /fill, 1.1
plot_image, bytscl(im1, 20, 200), title=target, color=cgcolor('black'), charsize=0.7
oplot, i2d[0,*], i2d[1,*], color=cgcolor('purple'), psym=2
plot, medcol2, color=cgcolor('black'), psym=-8, title=target+' total difference = '+strtrim(totdif,2)+' DN',$
    charsize=0.7, xtitle='Column number', ytitle='Median column DN'
oplot, medcol3, color=cgcolor('purple'), psym=-8
al_legend, ['all pixels', 'avoiding PSF'], color=[cgcolor('black'),cgcolor('purple')], psym=8, /box, /right, $
  charsize=0.8
ps_off
!p.multi=0

spawn, 'open '+fileout+' &'

end