pro stk_vs_non

; compare different elements of stacked versus non-stacked images - at various stages in the reduction pipeline
; 
Compile_opt idl2
!p.background=cgcolor('white')

indir='~/BRITE/data/UB/p3_2/'  ; HP cleaned + p2_trend removed

stkfiles=file_search(indir+'Orion-CF1-3*_p3.sav', count=nstk)
sinfiles=file_search(indir+'HD*', count=nsin)

window, 1, xsize=700, ysize=500
plotsym, 0, /fill, 1.1

vmags=fltarr(15)
sdiff=fltarr(15)

for i=0, 14 do begin
  
  restore, stkfiles[i]  ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
                        ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend, sig_trend
                        
  vpos=strpos(simbad_mags, 'V')
                        
  vmags[i]=float(strmid(simbad_mags, vpos+2))
                       
  ; calculate average scatter in background from ALL observations
  avg_sig1=median(sig_trend)  
  
  restore, sinfiles[i]  ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
  ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend, sig_trend
  
  ; calculate average scatter in background from ALL observations
  avg_sig2=median(sig_trend)
  
  sdiff[i]=(avg_sig2-avg_sig1)*100.
  
endfor

!p.multi=0
plot, vmags, sdiff, color=cgcolor('black'), psym=8



stop
print, 'end of program'
end