pro warsaw_medcol

  ; plot med col results for presentation

Compile_opt idl2

sat='UB'
field='CENTAURUS'

plotsym, 0, /fill, 1.

indir='~/BRITE/'+sat+'/'+field+'/reduction/medcols/'
outdir='~/BRITE/reports/presentations/warsaw/'

filein=indir+'medcol_results.txt'
fileout=outdir+'medcols.ps'

readcol, filein, fname, medw, medr, format='a,f,f,x,x,x', skipline=1

fname2=strmid(fname, 9)

keep=where(fname2 eq '28_28_0', nkeep)
fname1=strmid(fname[keep], 0, 8)

ps_on, fileout, xsize=15, ysize=10
plot, medw[keep], color=cgcolor('black'), psym=8, ytitle='Median of residual image', charsize=0.7, $
  xrange=[-2,32], xstyle=1, /nodata, xtitle='Target number'
oplot, medw[keep], color=cgcolor('dark grey'), psym=8
oplot, medr[keep], color=cgcolor('purple'), psym=8
al_legend, ['whole', 'raster'], psym=8, colors=[cgcolor('dark grey'),cgcolor('purple')], $
  /top, /right, textcolors=cgcolor('black'), charsize=0.7
ps_off

spawn, 'open '+fileout+' &'
stop


  ; 
  ;
print, 'end of program'
end