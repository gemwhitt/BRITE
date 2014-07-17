pro plot_lcs

compile_opt idl2

!p.background=cgcolor('white')

indir='/Users/gemmawhittaker/BRITE/data/UB/p4/'
outdir='/Users/gemmawhittaker/BRITE/results/gemphot/plots/lcs/'

filesin=file_search(indir+'*CF1-7*_gem2.sav', count=nsav)
fname=file_basename(filesin, '_p4_gem2.sav')

for ii=0, nsav-1 do begin

  restore, filesin[ii]  ;flux, jd, bkgd, bk_err, mean_resid, sig_resid, roi_name, $
  ;exp_num, ra_dec1, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
  ;  simbad_mags, parlax, otype, sptype, p2_trend, med_trend, nsubpix
  
  xx=where(flux ne 0.0 AND jd-jd[0] ge 0.0, nxx)
  flux=flux[xx]
  jd=jd[xx]
  counts=nsubpix[xx]
  
;  vpos=strpos(simbad_mags, 'V=')
;  vmag=strmid(simbad_mags, vpos)
  
  nflux=flux/robust_mean(flux,3)
  time1=jd-jd[0]
  ntime=n_elements(time1)
  
  xlim=min([time1[ntime-1],6])
;  stop
  
  ;window, 0, xpos=1500, ypos=200
  ;plot, time1, nflux, xrange=[0,7], psym=2, color=cgcolor('black')
  
  ncounts=counts/robust_mean(counts,3)
  sig_counts=robust_sigma(ncounts)
  
  width=float(n_elements(counts))/3.
  
  sm_trend=smooth(ncounts,width, /edge_truncate)
  
  ;window, 1, xpos=2300, ypos=200
  ;plot, time1, ncounts, xrange=[0,xlim], psym=2, color=cgcolor('black')
  ;oplot, time1, sm_trend, color=cgcolor('orange'), thick=2
  ;oplot, time1, sm_trend+(3*sig_counts), color=cgcolor('green'), thick=2
  ;oplot, time1, sm_trend-(3*sig_counts), color=cgcolor('green'), thick=2
  
  ;stop
  
  yy=where(ncounts gt sm_trend+(3*sig_counts) OR ncounts lt sm_trend-(3*sig_counts), nyy, complement=keep)
  
  if nyy gt 0 then begin
    time1=time1[keep]
    nflux=nflux[keep]
  endif
  
  ; smooth the flux over 1 day or less than this
  npts_tot=n_elements(time1)
  time2=time1[1:n_elements(time1)-1]
  tdiff=time2-time1
  gap=where(tdiff gt 0.015, ngap)
  
  sig_flux=robust_sigma(nflux[0:gap[0]])
  scat=strmid(strtrim(sig_flux,2), 0, 5)
  
  fileout=outdir+fname[ii]+'_lc2.ps'
  plotsym, 0, /fill, 0.4
  ;window, 0, xpos=1500, ypos=200
  ps_on, fileout, xsize=16, ysize=11
  plot, time1, nflux, xrange=[0,xlim], psym=8, color=cgcolor('black'), yrange=[0.8,1.2], $
    ;xtitle='Time (days)', ytitle='Normalized Flux', charsize=0.7, title=fname[ii]+', '+vmag+', scatter='+scat+'%'
    xtitle='Time (days)', ytitle='Normalized Flux', charsize=0.7, title=fname[ii]+', scatter='+scat+'%'
    ps_off
  
 ; print, sig_flux
  
  
  
;  stop
  
  

endfor

print, 'End of program'
end