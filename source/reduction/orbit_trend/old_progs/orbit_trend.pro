pro orbit_trend, expn,flx,temp1,nbins1,width, avg_flx,sm_flx,avg_temp

; Using the "phase-folded" time series, folded on a period of i-orbit, calculate the average trend accross the orbit

; bin the points by exposure number
result=histogram(expn, nbins=43, locations=loc, reverse_indices=ri)
  
avg_flx=fltarr(43)
avg_temp=fltarr(43)
  
for i=0, 42 do begin
  
  avg_flx[i]=robust_mean(flx[ri[ri[i]:ri[i+1]-1]],2)
  avg_temp[i]=robust_mean(temp1[ri[ri[i]:ri[i+1]-1]],2)
    
endfor

sm_flx=smooth(avg_flx,width, /edge_truncate)

end