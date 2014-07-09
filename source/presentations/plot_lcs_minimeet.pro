pro plot_lcs_minimeet

Compile_opt idl2

indir='~/BRITE/data/UB/testing/nostack/p4_cust/'

outdir='/Users/gemmawhittaker/BRITE/mini_meeting/plots/'

!p.background=cgcolor('white')

filesin=file_search(indir+'HD31237_p4_test.sav', count=nsav)

restore, filesin[0]
plotsym, 0, /fill, 0.5

flux=flux/robust_mean(flux,2)

fileout1=outdir+'gem_lc_2o.ps'
ps_on, fileout1, xsize=20, ysize=7
plot, jd-jd[0], flux, thick=2, xtitle='Time (days)', ytitle='Normalized Flux', charsize=0.7, title='HD31237', $
  color=cgcolor('black'), xrange=[0,0.09], yrange=[0.9,1.1], xstyle=1
  ps_off
  stop
  fileout2=outdir+'gem_lc_1o.ps'
  ps_on, fileout2, xsize=20, ysize=7
  plot, jd-jd[0], flux, psym=8, xtitle='Time (days)', ytitle='Normalized Flux', charsize=0.7, title='HD31237', $
    color=cgcolor('black'), xrange=[0,0.012], yrange=[0.9,1.1], xstyle=1
    ps_off
    
    xx=where(jd-jd[0] lt 0.012)
    f1=flux[xx]
    print, robust_sigma(f1)

stop


end