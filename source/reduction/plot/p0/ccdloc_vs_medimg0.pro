pro ccdloc_vs_medimg0

; Make a plot of CCD readout distance versus medimg0, plot start on left/right hand-side in different colours
;
Compile_opt idl2

sat='BA'
field='CENTAURUS'

indir='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/stats/p0/'

outdir='/Users/gemmawhittaker/BRITE/'+sat+'/'+field+'/plots/p0/'

filein=indir+'ccd_loc_vs_medimg0.txt'

plotout=outdir+'medimg0_v_xyccd.ps'

readcol, filein, fname, xloc, yloc, medimg, format='a,i,i,i'

a=where(yloc lt 1330)
b=where(yloc ge 1330)

; plot medimg vs xloc and yloc - save outputs


ps_on, plotout, xsize=16, ysize=26
!p.multi=[0,1,2,0,0]
plotsym, 0, /fill, 0.7
plot, xloc, medimg, color=cgcolor('black'), /ynozero, psym=8, title=sat+' '+field, $
  xtitle='CCD X-location', ytitle='Median raw images', charsize=0.7
plot, yloc, medimg, color=cgcolor('black'), /ynozero, psym=8, title=sat+' '+field, $
  xtitle='CCD Y-location', ytitle='Median raw images', charsize=0.7
ps_off
!p.multi=0

spawn, 'open '+plotout

print, 'End of program'
end