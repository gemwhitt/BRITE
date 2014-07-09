pro plot_p0

; plot p0 info 
; 
; 
Compile_opt idl2

!p.background=cgcolor('white')

indir='/Users/gemmawhittaker/BRITE/data/UB/roi_raw_sav/ORION/'

outdir='/Users/gemmawhittaker/BRITE/results/meeting_jun/'

filesin=file_search(indir+'*.sav', count=nf)

for i=0, nf-1 do begin
  
  restore, filesin[i] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl
  
  nfrm=n_elements(jd)
  
  totdn=fltarr(nfrm)
  medimg0=fltarr(nfrm)
  avg_temp=fltarr(nfrm)

  for j=0, nfrm-1 do begin
    totdn[j]=total(data1[*,*,j])
    medimg0[j]=median(data1[*,*,j])
    avg_temp[j]=average(ccd_temp[*,j])
  endfor
  
  jd1=jd-jd[0]
  
  plotsym, 0, /fill, 0.3
  ps_on, outdir+'medimg_all.ps', xsize=15, ysize=11
  plot, jd1, medimg0, color=cgcolor('black'), psym=8, xtitle='Days since start of observations', $
    ytitle='Median images - before warm column removal', charsize=0.7
  oplot, [0,200], [16400, 16400], color=cgcolor('blue'), thick=5
  oplot, [0,200], [5000,5000], color=cgcolor('red'), thick=5
  ps_off
      
  spawn, 'open '+outdir+'medimg_all.ps'
  
  stop
  
endfor




end

