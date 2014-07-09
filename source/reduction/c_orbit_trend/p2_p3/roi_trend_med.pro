pro roi_trend_med

; Read in result from calculate_trend_orbit.pro :
; (outfile='/Users/gemmawhittaker/BRITE/reduction/orbit_trends/roi_med_all/roi_med_all.sav')
; disregard results with roi_out=0 - which should take care of "off-target" images
; split the rest into observation sequences. - cross-correlate trends between orbits
; 
filein='/Users/gemmawhittaker/BRITE/reduction/orbit_trends/roi_med_all/roi_med_all_p3.sav'

restore, filein ;roi_out, jd1, header_temps, roi_loc

nroi=(size(roi_out, /dim))[0]

nfrm=n_elements(jd1)

jd2=jd1[1:nfrm-1]
jdiff=jd2-jd1
gap=where(jdiff gt 0.015, ngap) ; there are ngap+1 observation sequences - but not all have valid data points
gap=[-1,gap,nfrm-1]
count=roi_loc[0,*]+roi_loc[1,*]

for ii=0, 0 do begin  ;ngap do begin
  
  ; watch out for 0.0 values in roi_out - iterate over these?
  time1=jd1[gap[ii]+1:gap[ii+1]]
  roi_med=roi_out[*,gap[ii]+1:gap[ii+1]]
  header_temp=header_temps[gap[ii]+1:gap[ii+1]]
  
  npts=n_elements(time1)
  
  avg_med=fltarr(npts)
  for jj=0, npts-1 do avg_med[jj]=robust_mean(roi_med[*,jj],4)
  width=3.
  sm_med=smooth(avg_med, width, /edge_truncate)
  
  avg_roi=fltarr(nroi)
  for jj=0, nroi-1 do avg_roi[jj]=robust_mean(roi_med[jj,*],2)
  count=roi_loc[0,*]+roi_loc[1,*]
  
  plot, count, avg_roi, color=cgcolor('black'), psym=2
  stop
  
  ; record the trend
  outfile='/Users/gemmawhittaker/BRITE/reduction/orbit_trends/roi_med_all/roi_trend_all.sav'
  ;save, filename=outfile, time

  window, 1, xpos=1500, ypos=200, xsize=700, ysize=500
  plotsym, 0, /fill, 0.9
  plot, time1, roi_med[0,*], color=cgcolor('black'), psym=2, yrange=[80,120], /nodata
  for jj=0, nroi-1 do begin
    oplot, time1, roi_med[jj,*], color=cgcolor('orchid'), psym=8
;    stop
  endfor
  oplot, time1, avg_med, color=cgcolor('purple'), thick=3

  stop

endfor


; do x-correlations
lag=indgen(100)

for ii=0, 13 do begin

series1=reform(roi_out[ii,0:700])
series2=reform(roi_out[ii+1,0:700])

result2=a_correlate(series1, lag, /double)

window, 3, xpos=3100, ypos=-400, xsize=700, ysize=500
plot, lag, result2, color=cgcolor('black'), thick=2, xtitle='Lag', ytitle='Auto-Correlation Function', charsize=0.7

result=c_correlate(series1, series2, lag, /double)

;window, 2, xpos=3100, ypos=200, xsize=700, ysize=500
;plot, lag, result, color=cgcolor('black'), thick=2, xtitle='Lag', ytitle='Cross-Correlation Function', charsize=0.7

;result2=a_correlate(series1, lag, /double)

;window, 3, xpos=3100, ypos=-400, xsize=700, ysize=500
;plot, lag, result2, color=cgcolor('black'), thick=2, xtitle='Lag', ytitle='Auto-Correlation Function', charsize=0.7

stop
endfor







stop

end 