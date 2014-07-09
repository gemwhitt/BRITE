pro compare_p4results

  ; program written to check the results from p3_p4_custom_4 - using no binning
  ;
  ; check mags, number of counts, flux, residual, tracking...
  ;
  Compile_opt idl2
  
  !p.background=cgcolor('white')
  
  indir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot1/'
 
  filesin=file_search(indir+'*_gem.sav', count=nf)
    
  ; set up arrays
  lc_sig=fltarr(4,nf)
  npixels=fltarr(4,nf)
  vmags=fltarr(nf)
  
  for i=0, nf-1 do begin
  
  restore, filesin[i]  ; vmag, nsubpix, max_snr, max_dn, xy_psf, flux1, flux2, jd1, nsat, $
;    roi_name, exp_num, ra_dec1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
;      simbad_mags, parlax, otype, sptype, p2_trend, med_trend

vmags[i]=vmag

mags=-2.5*alog10(flux)

for j=0, 3 do begin
  npixels[j,i]=robust_mean(nsubpix[j,*],2)
  lc_sig[j,i]=robust_sigma(mags[j,*]/robust_mean(mags[j,*],2))
  
endfor

  endfor
  
  cols1=cgcolor(['black', 'purple', 'red', 'orange'])
  
  window, 1, xsize=550, ysize=500
  plotsym, 0, /fill, 0.9
  plot, vmags, lc_sig, /nodata, xtitle='V-magnitude', ytitle='Scatter', charsize=1., color=cgcolor('black'), $
    yrange=[0,0.004]
  for j=0, 3 do oplot, vmags, lc_sig[j,*], psym=8, color=cols1[j]
 
 
  stop
  print, 'End of program'
end

