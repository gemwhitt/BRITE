pro brite_autocor1

; Investigate correlations in the BRITE datasets
; Output: Autocorrelation graph

Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 1.3

sat='TOR'
field='CENTAURUS'
target='HD127973'

indir='~/BRITE/'+sat+'/'+field+'/data/aper_lc_final/' 
filesin=file_search(indir+target+'*', count=nf)

for f=0, nf-1 do begin
  
  ; write out new file with correction applied
  readcol, filesin[f], time, frame, flux, xcom, $
      ycom, temperature, vmag, bmag, resid,$
      format='(d,i,f,f,f,f,f,f,f)'

;restore, filesin[f]
;      stop
  ;time=jd
  time1=time-time[0]
  plotsym, 0, /fill, 0.8
  
  nfrm=n_elements(time)
  
 ; flux=fltarr(nfrm)
 ; for im=0, nfrm-1 do flux[im]=total(data1[*,*,im])
  
  ;temp=fltarr(nfrm)
  ;for im=0, nfrm-1 do temp[im]=total(ccd_temp[*,im])
  temp=temperature
  
  plot, time1, flux, psym=8, color=cgcolor('black')
  
; stop
  plot, time1, temp, psym=8, color=cgcolor('black')
  
  ; find where data becomes consistent
;  xx=where(time1 ge 70 AND time1 le 80, nxx)
;  
;  time2=time1[xx]
  flux2=flux;[xx]
;  temp=temperature[xx]
;  
  
  ; DO CROSS-CORRELATIONS....
  print, correlate(temp, flux2)
  lag=indgen(100)
  
  result=c_correlate(temp, flux2, lag)
  
  plot, result, color=cgcolor('black'), psym=8
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  stop
  ; DO AUTO-CORRELATIONS.......
  ac=fltarr(n_elements(lag))  
  for i=0, n_elements(lag)-1 do begin 
    nitt=n1-lag[i]
    ac1=fltarr(nitt)
    
    for k=0, nitt-1 do ac1[k]=(totdn1[k] - mean(totdn1)) * (totdn1[k+lag[i]] - mean(totdn1))
        
    ac[i]=total(ac1)*(1./(n1-lag[i]))
  endfor
  
  result=ac/ac[0]  
  
  plot, result, color=cgcolor('black'), psym=8
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;Now what??
  
  stop
  
endfor

end