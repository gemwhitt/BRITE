pro ccdloc_vs_medimg0

; Make a plot of CCD readout distance versus medimg0, plot start on left/right hand-side in different colours
;
Compile_opt idl2

sat='BA'
field='ORION'

indir='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/stats/p0_images/'

filein=indir+'ccd_loc_vs_medimg0.txt'

readcol, filein, fname, read_dist, xloc, yloc, med_img0, avg_img0, format='a,i,i,i,i,f)'

a=where(yloc lt 1330)
b=where(yloc ge 1330)

plotsym, 0, /fill, 1.0
plot, read_dist, med_img0, color=cgcolor('black'), /ynozero, /nodata
oplot, read_dist[a], med_img0[a], color=cgcolor('blue'), psym=8
oplot, read_dist[b], med_img0[b], color=cgcolor('green'), psym=8

stop

end